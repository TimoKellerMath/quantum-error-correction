/*
 * exactB: exhaustive search for light nontrivial logical operators of a CSS code with checks of
 * any (mixed) weight, any number of qubits and any number k of logical qubits.
 *
 * Implementation B of Track T-E (new_codes), generalized from src/cluster2.c of
 * quantum-error-correction/weight6-exact-distances.  "ED" below is exact_distances.tex there;
 * NOTES.md maps every part of this file to a statement of ED and lists what changes for general
 * weights.
 *
 * Side Z: A = H_X (the logicals are z in ker H_X), B = H_Z (z is trivial iff z in rowsp H_Z).
 * Side X: A = H_Z, B = H_X.
 *
 *   search      ED, Algorithm alg:search; completeness ED, Proposition prop:exact
 *   P_c         all even subsets of c, so the first rule of PROPAGATE is XOR unit propagation
 *               (NOTES.md, Lemma N1)
 *   S           rows of B and sums of two overlapping rows of B (ED, Lemma lem:minimal(b))
 *   beta        max of the degree-sequence bound (NOTES.md, Lemma N2; it dominates ceil(u/h) of
 *               ED) and a greedy packing of unsatisfied checks with disjoint free sets (ED)
 *   anchors     ED, Lemma lem:anchors and Corollary cor:runs: orbits of a group whose generators
 *               are verified to preserve rowsp A and rowsp B; without a group, single qubits
 *   parity      ED, Lemma lem:parity in the form "1 in rowsp A" (NOTES.md, Lemma N3)
 *   triviality  z in ker A is trivial iff <z, L_j> = 0 for a basis L_j of ker B modulo rowsp A
 *               (ED, Section 2, first paragraph)
 */
#define _GNU_SOURCE
#include <errno.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXW 64      /* largest check weight of A accepted (candidate array) */
#define INF (1 << 28)

static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "exactB: error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(2);
}

static void *xmalloc(size_t s) { void *p = malloc(s ? s : 1); if (!p) die("out of memory"); return p; }
static void *xcalloc(size_t a, size_t b) { void *p = calloc(a ? a : 1, b ? b : 1); if (!p) die("out of memory"); return p; }

static double wall(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static double cputime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int cmp_int(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

/* ------------------------------------------------------------------ sparse matrices (CSR) */

typedef struct { int rows; int *off, *idx; } Csr;

/* "n m" on the first line, then m lines with the 0-based qubits of one check each */
static Csr read_checks(const char *fn, int *n_out) {
    FILE *f = fopen(fn, "r");
    if (!f) die("cannot open %s", fn);
    char *line = NULL;
    size_t cap = 0;
    long nn, mm;
    if (getline(&line, &cap, f) < 0) die("%s: empty file", fn);
    char tail[8];
    if (sscanf(line, "%ld %ld %7s", &nn, &mm, tail) != 2 || nn <= 0 || mm < 0 || nn > 100000000L || mm > 100000000L)
        die("%s: first line must be \"n m\"", fn);
    Csr M;
    M.rows = (int)mm;
    M.off = xmalloc((mm + 1) * sizeof(int));
    size_t icap = 1024, ilen = 0;
    M.idx = xmalloc(icap * sizeof(int));
    char *seen = xcalloc(nn, 1);
    for (long r = 0; r < mm; r++) {
        if (getline(&line, &cap, f) < 0) die("%s: %ld check lines expected, found %ld", fn, mm, r);
        M.off[r] = (int)ilen;
        char *p = line, *e;
        for (;;) {
            while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
            if (!*p) break;
            errno = 0;
            long v = strtol(p, &e, 10);
            if (e == p || errno) die("%s: line %ld: not an integer", fn, r + 2);
            if (v < 0 || v >= nn) die("%s: line %ld: qubit %ld out of range", fn, r + 2, v);
            if (seen[v]) die("%s: line %ld: qubit %ld repeated", fn, r + 2, v);
            seen[v] = 1;
            if (ilen == icap) { icap *= 2; M.idx = realloc(M.idx, icap * sizeof(int)); if (!M.idx) die("out of memory"); }
            M.idx[ilen++] = (int)v;
            p = e;
        }
        for (size_t i = M.off[r]; i < ilen; i++) seen[M.idx[i]] = 0;
        qsort(M.idx + M.off[r], ilen - M.off[r], sizeof(int), cmp_int);
    }
    M.off[mm] = (int)ilen;
    while (getline(&line, &cap, f) >= 0)
        for (char *p = line; *p; p++)
            if (*p != ' ' && *p != '\t' && *p != '\r' && *p != '\n') die("%s: more than %ld check lines", fn, mm);
    free(line);
    free(seen);
    fclose(f);
    *n_out = (int)nn;
    return M;
}

static void transpose(const Csr *M, int n, int **toff, int **tval) {
    int *cnt = xcalloc(n + 1, sizeof(int));
    for (int i = 0; i < M->off[M->rows]; i++) cnt[M->idx[i] + 1]++;
    for (int q = 0; q < n; q++) cnt[q + 1] += cnt[q];
    int *tv = xmalloc((M->off[M->rows] + 1) * sizeof(int));
    int *pos = xmalloc((n + 1) * sizeof(int));
    memcpy(pos, cnt, (n + 1) * sizeof(int));
    for (int r = 0; r < M->rows; r++)
        for (int i = M->off[r]; i < M->off[r + 1]; i++) tv[pos[M->idx[i]]++] = r;
    free(pos);
    *toff = cnt;
    *tval = tv;
}

/* ------------------------------------------------------------------ dense GF(2) */

typedef struct { int r, c, nw; uint64_t *a; } Bm;
#define ROW(M, i) ((M)->a + (size_t)(i) * (M)->nw)
#define GETB(v, j) (((v)[(j) >> 6] >> ((j) & 63)) & 1)
#define SETB(v, j) ((v)[(j) >> 6] |= 1ULL << ((j) & 63))

static Bm bm_new(int r, int c) {
    Bm M;
    M.r = r; M.c = c; M.nw = (c + 63) / 64;
    M.a = xcalloc((size_t)(r > 0 ? r : 1) * M.nw, sizeof(uint64_t));
    return M;
}
static Bm bm_from_csr(const Csr *S, int n) {
    Bm M = bm_new(S->rows, n);
    for (int r = 0; r < S->rows; r++)
        for (int i = S->off[r]; i < S->off[r + 1]; i++) SETB(ROW(&M, r), S->idx[i]);
    return M;
}
static inline void vxor(uint64_t *a, const uint64_t *b, int nw) { for (int w = 0; w < nw; w++) a[w] ^= b[w]; }
static inline int viszero(const uint64_t *a, int nw) { for (int w = 0; w < nw; w++) if (a[w]) return 0; return 1; }

/* reduced row echelon form in place; piv[i] is the pivot column of row i < rank */
static int bm_rref(Bm *M, int *piv) {
    int rank = 0;
    uint64_t *tmp = xmalloc(M->nw * sizeof(uint64_t));
    for (int col = 0; col < M->c && rank < M->r; col++) {
        int w = col >> 6;
        uint64_t bit = 1ULL << (col & 63);
        int p = -1;
        for (int i = rank; i < M->r; i++) if (ROW(M, i)[w] & bit) { p = i; break; }
        if (p < 0) continue;
        if (p != rank) {
            memcpy(tmp, ROW(M, p), M->nw * 8);
            memcpy(ROW(M, p), ROW(M, rank), M->nw * 8);
            memcpy(ROW(M, rank), tmp, M->nw * 8);
        }
        uint64_t *pr = ROW(M, rank);
        for (int i = 0; i < M->r; i++) {
            if (i == rank) continue;
            uint64_t *ri = ROW(M, i);
            if (ri[w] & bit) for (int k = w; k < M->nw; k++) ri[k] ^= pr[k];
        }
        piv[rank++] = col;
    }
    free(tmp);
    return rank;
}

/* v <- v + sum of the RREF rows i with v[piv[i]] = 1; the result is 0 iff v was in the row space */
static void bm_reduce(const Bm *R, int rank, const int *piv, uint64_t *v) {
    for (int i = 0; i < rank; i++) if (GETB(v, piv[i])) vxor(v, ROW(R, i), R->nw);
}

/* ------------------------------------------------------------------ problem data (shared) */

static int n, mA, mB, K, clsw, cmax, wAmax, wBmax, rkA, rkB, parity;
static Csr A, B;
static int *qaoff, *qa;                  /* qubit -> A-checks */
static int ns, *soff, *sq, *sthr;        /* majority family S */
static int *qsoff, *qs;                  /* qubit -> members of S */
static uint64_t *cls;                    /* n x clsw, bit j of row q is L_j[q] */
static Bm RA, RB;                        /* RREF of A and B */
static int *pivA, *pivB;
static int W, Wsearch;
static int opt_bound = 3, opt_majo = 2, opt_stabw = 0, nthreads = 1, splitD = -1, verbose = 0;
static long long node_limit = -1;
static int nruns, *run_start, *ord, only_run = -1;
static atomic_int stopflag;
static atomic_llong next_job;
static atomic_llong total_nodes;
static pthread_mutex_t found_mtx = PTHREAD_MUTEX_INITIALIZER;
static int found_w = -1, *found_set;

/* ------------------------------------------------------------------ majority family */

static void push_stab(const int *v, int w, int *cap_s, int *cap_q) {
    if (ns + 1 >= *cap_s) {
        *cap_s *= 2;
        soff = realloc(soff, (*cap_s + 1) * sizeof(int));
        sthr = realloc(sthr, *cap_s * sizeof(int));
        if (!soff || !sthr) die("out of memory");
    }
    while (soff[ns] + w > *cap_q) { *cap_q *= 2; sq = realloc(sq, *cap_q * sizeof(int)); if (!sq) die("out of memory"); }
    memcpy(sq + soff[ns], v, w * sizeof(int));
    sthr[ns] = w / 2;
    soff[ns + 1] = soff[ns] + w;
    ns++;
}

/* S = rows of B (majo >= 1) and sums of two rows of B with a common qubit (majo >= 2), weight <= wmax;
 * every member lies in rowsp B, which is all that ED, Lemma lem:minimal(b) needs */
static void build_stabs(int majo, int wmax) {
    int cap_s = 1024, cap_q = 8192;
    soff = xmalloc((cap_s + 1) * sizeof(int));
    sthr = xmalloc(cap_s * sizeof(int));
    sq = xmalloc(cap_q * sizeof(int));
    ns = 0;
    soff[0] = 0;
    if (majo >= 1)
        for (int r = 0; r < mB; r++) {
            int w = B.off[r + 1] - B.off[r];
            if (w > 0 && w <= wmax) push_stab(B.idx + B.off[r], w, &cap_s, &cap_q);
        }
    if (majo >= 2) {
        int *qbo, *qbv;
        transpose(&B, n, &qbo, &qbv);
        int *mk = xmalloc(mB * sizeof(int));
        for (int r = 0; r < mB; r++) mk[r] = -1;
        int *buf = xmalloc((2 * wBmax + 2) * sizeof(int));
        for (int r = 0; r < mB; r++)
            for (int i = B.off[r]; i < B.off[r + 1]; i++)
                for (int j = qbo[B.idx[i]]; j < qbo[B.idx[i] + 1]; j++) {
                    int t = qbv[j];
                    if (t <= r || mk[t] == r) continue;
                    mk[t] = r;
                    int a = B.off[r], ae = B.off[r + 1], b = B.off[t], be = B.off[t + 1], w = 0;
                    while (a < ae || b < be) {
                        if (b == be || (a < ae && B.idx[a] < B.idx[b])) buf[w++] = B.idx[a++];
                        else if (a == ae || B.idx[b] < B.idx[a]) buf[w++] = B.idx[b++];
                        else { a++; b++; }
                    }
                    if (w > 0 && w <= wmax) push_stab(buf, w, &cap_s, &cap_q);
                }
        free(buf); free(mk); free(qbo); free(qbv);
    }
    Csr S = { ns, soff, sq };
    transpose(&S, n, &qsoff, &qs);
}

/* ------------------------------------------------------------------ search state (per thread) */

typedef struct {
    int tid;
    signed char *st;                     /* 0 free, 1 in T, 2 in X (excluded) */
    unsigned char *syn;
    int *avail, *ulist, *upos, ucount;
    uint64_t *ubits;
    int *hq, *cntH;                      /* hq[q]: unsatisfied A-checks through q; cntH[h]: free q with hq = h */
    int *scnt;
    int *trail, tlen, *tin, Tsize;
    int *cstack, ctop, *sstack, stop_;
    int *mark, markgen;
    int *pk1, *pk2, *pkc;
    uint64_t *cbuf;
    long long jobseq, myclaim;
    int cur_run;
    long long nodes, generated, leaves;
    long long *runnodes;
} Ctx;

static void ctx_init(Ctx *x, int tid) {
    memset(x, 0, sizeof *x);
    x->tid = tid;
    x->st = xcalloc(n, 1);
    x->syn = xcalloc(mA, 1);
    x->avail = xmalloc(mA * sizeof(int));
    for (int c = 0; c < mA; c++) x->avail[c] = A.off[c + 1] - A.off[c];
    x->ulist = xmalloc(mA * sizeof(int));
    x->upos = xmalloc(mA * sizeof(int));
    x->ubits = xcalloc((mA + 63) / 64 + 1, sizeof(uint64_t));
    x->hq = xcalloc(n, sizeof(int));
    x->cntH = xcalloc(cmax + 2, sizeof(int));
    x->cntH[0] = n;
    x->scnt = xcalloc(ns + 1, sizeof(int));
    x->trail = xmalloc((n + 1) * sizeof(int));
    x->tin = xmalloc((n + 2) * sizeof(int));
    x->cstack = xmalloc(((size_t)A.off[mA] + MAXW + 64) * sizeof(int));
    x->sstack = xmalloc(((size_t)soff[ns] + 64) * sizeof(int));
    x->mark = xcalloc(n, sizeof(int));
    x->pk1 = xmalloc((mA + 1) * sizeof(int));
    x->pk2 = xmalloc((mA + 1) * sizeof(int));
    x->pkc = xmalloc((wAmax + 2) * sizeof(int));
    x->cbuf = xmalloc(clsw * sizeof(uint64_t));
    x->runnodes = xcalloc(nruns, sizeof(long long));
}

static inline void flip(Ctx *x, int c) {
    int d;
    x->syn[c] ^= 1;
    if (x->syn[c]) {
        x->upos[c] = x->ucount;
        x->ulist[x->ucount++] = c;
        x->ubits[c >> 6] |= 1ULL << (c & 63);
        d = 1;
    } else {
        int p = x->upos[c], l = x->ulist[--x->ucount];
        x->ulist[p] = l;
        x->upos[l] = p;
        x->ubits[c >> 6] &= ~(1ULL << (c & 63));
        d = -1;
    }
    for (int i = A.off[c]; i < A.off[c + 1]; i++) {
        int q = A.idx[i];
        if (x->st[q] == 0) { x->cntH[x->hq[q]]--; x->hq[q] += d; x->cntH[x->hq[q]]++; }
        else x->hq[q] += d;
    }
}

/* put q into T; returns 0 on an immediate contradiction (the assignment is on the trail anyway) */
static int set_in(Ctx *x, int q) {
    int ok = 1;
    x->cntH[x->hq[q]]--;
    x->st[q] = 1;
    x->trail[x->tlen++] = q;
    x->tin[x->Tsize++] = q;
    for (int i = qaoff[q]; i < qaoff[q + 1]; i++) {
        int c = qa[i];
        flip(x, c);
        x->avail[c]--;
        x->cstack[x->ctop++] = c;
    }
    for (int i = qsoff[q]; i < qsoff[q + 1]; i++) {
        int s = qs[i];
        x->scnt[s]++;
        if (x->scnt[s] > sthr[s]) ok = 0;
        else if (x->scnt[s] == sthr[s]) x->sstack[x->stop_++] = s;
    }
    if (x->Tsize > Wsearch) ok = 0;
    return ok;
}

static void set_out(Ctx *x, int q) {
    x->cntH[x->hq[q]]--;
    x->st[q] = 2;
    x->trail[x->tlen++] = q;
    for (int i = qaoff[q]; i < qaoff[q + 1]; i++) {
        int c = qa[i];
        x->avail[c]--;
        x->cstack[x->ctop++] = c;
    }
}

static void undo(Ctx *x, int m) {
    while (x->tlen > m) {
        int q = x->trail[--x->tlen];
        if (x->st[q] == 1) {
            x->Tsize--;
            for (int i = qaoff[q + 1] - 1; i >= qaoff[q]; i--) { int c = qa[i]; flip(x, c); x->avail[c]++; }
            for (int i = qsoff[q]; i < qsoff[q + 1]; i++) x->scnt[qs[i]]--;
        } else {
            for (int i = qaoff[q]; i < qaoff[q + 1]; i++) x->avail[qa[i]]++;
        }
        x->st[q] = 0;
        x->cntH[x->hq[q]]++;
    }
}

/* Lower bound beta for |E|, E = z \ T, which consists of free qubits and meets every unsatisfied
 * check (ED, Lemma lem:minimal(c)).  All quantities depend only on the sets T, X (not on the
 * history), which the thread distribution below needs. */
static int lower_bound(Ctx *x) {
    int u = x->ucount;
    if (u == 0) return 0;
    int lb;
    if (opt_bound == 1) {
        int hmax = 0;
        for (int h = cmax; h >= 1; h--) if (x->cntH[h] > 0) { hmax = h; break; }
        if (hmax == 0) return INF;
        lb = (u + hmax - 1) / hmax;
    } else {
        /* NOTES.md, Lemma N2: sum over E of hq >= u; take the largest hq first */
        int need = u, t = 0;
        for (int h = cmax; h >= 1 && need > 0; h--) {
            int c = x->cntH[h];
            if (!c) continue;
            int m = (need + h - 1) / h;
            if (m > c) m = c;
            t += m;
            need -= m * h;
        }
        if (need > 0) return INF;
        lb = t;
    }
    if (opt_bound >= 2 && u > lb && lb <= Wsearch - x->Tsize) {
        /* unsatisfied checks with pairwise disjoint free sets need distinct elements of E; greedy in
         * the order (free count, index), which depends only on the state */
        int m = 0;
        for (int wd = 0; wd <= (mA - 1) >> 6; wd++) {
            uint64_t b = x->ubits[wd];
            while (b) { x->pk1[m++] = (wd << 6) + __builtin_ctzll(b); b &= b - 1; }
        }
        for (int a = 0; a <= wAmax + 1; a++) x->pkc[a] = 0;
        for (int i = 0; i < m; i++) x->pkc[x->avail[x->pk1[i]] + 1]++;
        for (int a = 1; a <= wAmax + 1; a++) x->pkc[a] += x->pkc[a - 1];
        for (int i = 0; i < m; i++) x->pk2[x->pkc[x->avail[x->pk1[i]]]++] = x->pk1[i];
        x->markgen++;
        int pk = 0;
        for (int i = 0; i < m; i++) {
            int c = x->pk2[i], ok = 1;
            for (int j = A.off[c]; j < A.off[c + 1]; j++) { int q = A.idx[j]; if (!x->st[q] && x->mark[q] == x->markgen) { ok = 0; break; } }
            if (!ok) continue;
            pk++;
            for (int j = A.off[c]; j < A.off[c + 1]; j++) { int q = A.idx[j]; if (!x->st[q]) x->mark[q] = x->markgen; }
        }
        if (pk > lb) lb = pk;
    }
    return lb;
}

/* PROPAGATE of ED, Algorithm alg:search, with P_c = all even subsets of c (NOTES.md, Lemma N1) */
static int propagate(Ctx *x, int ok) {
    if (!ok) goto fail;
    while (x->ctop > 0 || x->stop_ > 0) {
        if (x->stop_ > 0) {
            int s = x->sstack[--x->stop_];
            if (x->scnt[s] > sthr[s]) goto fail;
            if (x->scnt[s] == sthr[s])
                for (int i = soff[s]; i < soff[s + 1]; i++) if (x->st[sq[i]] == 0) set_out(x, sq[i]);
            continue;
        }
        int c = x->cstack[--x->ctop];
        int a = x->avail[c];
        if (a == 0) { if (x->syn[c]) goto fail; continue; }
        if (a != 1) continue;
        int q = -1;
        for (int i = A.off[c]; i < A.off[c + 1]; i++) if (x->st[A.idx[i]] == 0) { q = A.idx[i]; break; }
        if (x->syn[c]) { if (!set_in(x, q)) goto fail; }
        else set_out(x, q);
    }
    if (x->Tsize + lower_bound(x) > Wsearch) return 0;
    return 1;
fail:
    x->ctop = 0;
    x->stop_ = 0;
    return 0;
}

static int nontrivial(Ctx *x) {
    memset(x->cbuf, 0, clsw * sizeof(uint64_t));
    for (int i = 0; i < x->Tsize; i++) vxor(x->cbuf, cls + (size_t)x->tin[i] * clsw, clsw);
    return !viszero(x->cbuf, clsw);
}

static void report(Ctx *x) {
    pthread_mutex_lock(&found_mtx);
    if (found_w < 0) {
        found_w = x->Tsize;
        memcpy(found_set, x->tin, x->Tsize * sizeof(int));
        atomic_store(&stopflag, 1);
    }
    pthread_mutex_unlock(&found_mtx);
}

static int choose_check(Ctx *x) {
    int best = -1, ba = INF;
    for (int k = 0; k < x->ucount; k++) {
        int c = x->ulist[k], a = x->avail[c];
        if (a < ba || (a == ba && c < best)) { ba = a; best = c; }
    }
    return best;
}

static void search(Ctx *x, int depth);

/* Thread distribution: every thread walks the same tree above depth splitD; the depth-splitD
 * subtrees are numbered in walk order, and each number is executed by the one thread that
 * claimed it from the shared counter (NOTES.md, Lemma N4). */
static void subtree(Ctx *x, int depth) {
    if (nthreads > 1 && depth == splitD) {
        long long j = x->jobseq++;
        if (j != x->myclaim) return;
        search(x, depth);
        x->myclaim = atomic_fetch_add(&next_job, 1);
        return;
    }
    search(x, depth);
}

/* GROW of ED, Algorithm alg:search */
static void search(Ctx *x, int depth) {
    if (atomic_load_explicit(&stopflag, memory_order_relaxed)) return;
    int counted = (nthreads == 1 || depth >= splitD || x->tid == 0);
    if (counted) {
        x->nodes++;
        x->runnodes[x->cur_run]++;
        if (node_limit > 0 && (x->nodes & 1023) == 0) {
            long long tot = atomic_fetch_add(&total_nodes, 1024) + 1024;
            if (tot >= node_limit) { atomic_store(&stopflag, 2); return; }
        }
    }
    if (x->ucount == 0) {
        if (counted) x->leaves++;
        if (nontrivial(x)) report(x);
        return;
    }
    int c = choose_check(x);
    int cand[MAXW], r = 0;
    for (int i = A.off[c]; i < A.off[c + 1]; i++) if (x->st[A.idx[i]] == 0) cand[r++] = A.idx[i];
    int m0 = x->tlen;
    for (int i = 0; i < r; i++) {
        if (atomic_load_explicit(&stopflag, memory_order_relaxed)) break;
        int q = cand[i];
        if (x->st[q] == 2) continue;
        int already = (x->st[q] == 1);
        int m2 = x->tlen;
        if (counted) x->generated++;
        if (already || propagate(x, set_in(x, q))) subtree(x, depth + 1);
        undo(x, m2);
        if (already) break;
        if (i < r - 1) {
            set_out(x, q);
            if (!propagate(x, 1)) break;
        }
    }
    undo(x, m0);
}

static double *run_time;   /* wall time per run, single thread only */
static atomic_char *run_done; /* run completed by every thread that walked it */

static void *worker(void *arg) {
    Ctx *x = arg;
    x->myclaim = atomic_fetch_add(&next_job, 1);
    for (int r = 0; r < nruns; r++) {
        if (only_run >= 0 && r != only_run) continue;
        if (atomic_load(&stopflag)) break;
        double t0 = wall();
        x->cur_run = r;
        int s0 = run_start[r];
        for (int i = 0; i < s0; i++) if (x->st[ord[i]] == 0) set_out(x, ord[i]);
        int ok = propagate(x, 1);
        int a = ord[s0];
        if (ok) {
            if (x->st[a] == 0) ok = propagate(x, set_in(x, a));
            else ok = (x->st[a] == 1);
        }
        if (ok) subtree(x, 0);
        undo(x, 0);
        if (nthreads == 1) run_time[r] = wall() - t0;
        if (!atomic_load(&stopflag)) atomic_store(&run_done[r], 1);
    }
    return NULL;
}

/* ------------------------------------------------------------------ setup */

static void usage(void) {
    fprintf(stderr,
            "usage: exactB NAME.hx NAME.hz --side Z|X --W W [--group GENS] [--threads t]\n"
            "       [--order id|bfs] [--orbit-order min|size] [--split D] [--bound 1|2|3]\n"
            "       [--majority 0|1|2] [--stab-wmax w] [--run r] [--node-limit N] [--out FILE] [-v]\n"
            "  searches for a nontrivial Z-logical (side Z: in ker H_X, not in rowsp H_Z) or X-logical\n"
            "  of weight <= W; anchors: orbits of the group generated by GENS (one permutation of\n"
            "  0..n-1 per line, image form), else every qubit; earlier orbits are excluded.\n");
    exit(2);
}

/* generators: one line per permutation, n integers p(0) .. p(n-1); '#' starts a comment line */
static int read_gens(const char *fn, int **gens_out) {
    FILE *f = fopen(fn, "r");
    if (!f) die("cannot open %s", fn);
    int cap = 4, ng = 0;
    int *g = xmalloc((size_t)cap * n * sizeof(int));
    char *line = NULL;
    size_t lcap = 0;
    long ln = 0;
    while (getline(&line, &lcap, f) >= 0) {
        ln++;
        char *p = line, *e;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == 0) continue;
        if (ng == cap) { cap *= 2; g = realloc(g, (size_t)cap * n * sizeof(int)); if (!g) die("out of memory"); }
        int *perm = g + (size_t)ng * n, len = 0;
        for (;;) {
            while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
            if (!*p) break;
            errno = 0;
            long v = strtol(p, &e, 10);
            if (e == p || errno) die("%s: line %ld: not an integer", fn, ln);
            if (len >= n) die("%s: line %ld: more than n = %d entries", fn, ln, n);
            if (v < 0 || v >= n) die("%s: line %ld: entry %ld out of range", fn, ln, v);
            perm[len++] = (int)v;
            p = e;
        }
        if (len != n) die("%s: line %ld: %d entries, expected n = %d", fn, ln, len, n);
        char *seen = xcalloc(n, 1);
        for (int i = 0; i < n; i++) { if (seen[perm[i]]) die("%s: line %ld: not a permutation", fn, ln); seen[perm[i]] = 1; }
        free(seen);
        ng++;
    }
    free(line);
    fclose(f);
    *gens_out = g;
    return ng;
}

/* gamma(row) in rowsp M for every row of S (M = RREF of S), i.e. gamma(rowsp S) = rowsp S */
static int preserves(const int *perm, const Csr *S, const Bm *R, int rank, const int *piv) {
    uint64_t *v = xmalloc(R->nw * sizeof(uint64_t));
    int ok = 1;
    for (int r = 0; r < S->rows && ok; r++) {
        memset(v, 0, R->nw * sizeof(uint64_t));
        for (int i = S->off[r]; i < S->off[r + 1]; i++) SETB(v, perm[S->idx[i]]);
        bm_reduce(R, rank, piv, v);
        if (!viszero(v, R->nw)) ok = 0;
    }
    free(v);
    return ok;
}

static int uf_find(int *p, int a) { while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; } return a; }

int main(int argc, char **argv) {
    const char *fhx = NULL, *fhz = NULL, *fgens = NULL, *forder = "id", *foorder = "min", *fout = NULL;
    char side = 0;
    W = -1;
    for (int i = 1; i < argc; i++) {
        const char *o = argv[i];
#define ARG (i + 1 < argc ? argv[++i] : (usage(), (char *)0))
        if (!strcmp(o, "--side")) { const char *s = ARG; if (strcmp(s, "Z") && strcmp(s, "X")) usage(); side = s[0]; }
        else if (!strcmp(o, "--W")) W = atoi(ARG);
        else if (!strcmp(o, "--group")) fgens = ARG;
        else if (!strcmp(o, "--threads")) nthreads = atoi(ARG);
        else if (!strcmp(o, "--order")) forder = ARG;
        else if (!strcmp(o, "--orbit-order")) foorder = ARG;
        else if (!strcmp(o, "--split")) splitD = atoi(ARG);
        else if (!strcmp(o, "--bound")) opt_bound = atoi(ARG);
        else if (!strcmp(o, "--majority")) opt_majo = atoi(ARG);
        else if (!strcmp(o, "--stab-wmax")) opt_stabw = atoi(ARG);
        else if (!strcmp(o, "--run")) only_run = atoi(ARG);
        else if (!strcmp(o, "--node-limit")) node_limit = atoll(ARG);
        else if (!strcmp(o, "--out")) fout = ARG;
        else if (!strcmp(o, "-v")) verbose = 1;
        else if (o[0] == '-') usage();
        else if (!fhx) fhx = o;
        else if (!fhz) fhz = o;
        else usage();
    }
    if (!fhx || !fhz || !side || W < 0 || nthreads < 1 || opt_bound < 1 || opt_bound > 3 || opt_majo < 0 || opt_majo > 2) usage();
    if (nthreads > 1 && splitD < 0) splitD = 4;
    if (nthreads == 1) splitD = -1;

    double t0 = wall(), c0 = cputime();
    printf("cmd:");
    for (int i = 0; i < argc; i++) printf(" %s", argv[i]);
    printf("\n");
    int nx, nz;
    Csr HX = read_checks(fhx, &nx), HZ = read_checks(fhz, &nz);
    if (nx != nz) die("n differs: %d in %s, %d in %s", nx, fhx, nz, fhz);
    n = nx;
    if (side == 'Z') { A = HX; B = HZ; } else { A = HZ; B = HX; }
    mA = A.rows;
    mB = B.rows;
    wAmax = wBmax = 0;
    int wAmin = INF, wBmin = INF;
    for (int r = 0; r < mA; r++) { int w = A.off[r + 1] - A.off[r]; if (w > wAmax) wAmax = w; if (w < wAmin) wAmin = w; }
    for (int r = 0; r < mB; r++) { int w = B.off[r + 1] - B.off[r]; if (w > wBmax) wBmax = w; if (w < wBmin) wBmin = w; }
    if (wAmax > MAXW) die("check weight %d > %d", wAmax, MAXW);
    transpose(&A, n, &qaoff, &qa);
    cmax = 0;
    int cmin = INF;
    for (int q = 0; q < n; q++) { int c = qaoff[q + 1] - qaoff[q]; if (c > cmax) cmax = c; if (c < cmin) cmin = c; }

    /* linear algebra: CSS condition, ranks, logical basis, parity */
    RA = bm_from_csr(&A, n);
    RB = bm_from_csr(&B, n);
    {
        Bm DA = bm_from_csr(&A, n);
        for (int r = 0; r < mA; r++)
            for (int s = 0; s < mB; s++) {
                int par = 0;
                for (int i = B.off[s]; i < B.off[s + 1]; i++) par ^= GETB(ROW(&DA, r), B.idx[i]);
                if (par) die("H_X H_Z^T != 0 (rows %d, %d)", r, s);
            }
        free(DA.a);
    }
    pivA = xmalloc((mA + 1) * sizeof(int));
    pivB = xmalloc((mB + 1) * sizeof(int));
    rkA = bm_rref(&RA, pivA);
    rkB = bm_rref(&RB, pivB);
    K = n - rkA - rkB;
    clsw = K > 0 ? (K + 63) / 64 : 1;
    cls = xcalloc((size_t)n * clsw, sizeof(uint64_t));
    {
        /* L_1..L_K: kernel vectors of B (from its RREF) that extend the RREF rows of A to a basis of
         * ker B; the rows of A lie in ker B by the CSS condition */
        int nw = RA.nw;
        Bm bas = bm_new(rkA + (n - rkB) + 1, n);
        int *bpiv = xmalloc((rkA + (n - rkB) + 1) * sizeof(int));
        int nb = 0;
        for (int i = 0; i < rkA; i++) { memcpy(ROW(&bas, nb), ROW(&RA, i), nw * 8); bpiv[nb++] = pivA[i]; }
        char *isp = xcalloc(n, 1);
        for (int i = 0; i < rkB; i++) isp[pivB[i]] = 1;
        uint64_t *v = xmalloc(nw * 8);
        int nl = 0;
        for (int f = 0; f < n; f++) {
            if (isp[f]) continue;
            memset(v, 0, nw * 8);
            SETB(v, f);
            for (int i = 0; i < rkB; i++) if (GETB(ROW(&RB, i), f)) SETB(v, pivB[i]);
            for (int s = 0; s < mB; s++) {   /* v in ker B */
                int par = 0;
                for (int j = B.off[s]; j < B.off[s + 1]; j++) par ^= GETB(v, B.idx[j]);
                if (par) die("internal: kernel vector not in ker B");
            }
            uint64_t *u = ROW(&bas, nb);
            memcpy(u, v, nw * 8);
            for (int i = 0; i < nb; i++) if (GETB(u, bpiv[i])) vxor(u, ROW(&bas, i), nw);
            if (viszero(u, nw)) continue;
            int p = 0;
            while (!GETB(u, p)) p++;
            bpiv[nb++] = p;
            if (nl >= K) die("internal: more than k logicals");
            for (int q = 0; q < n; q++) if (GETB(u, q)) cls[(size_t)q * clsw + (nl >> 6)] |= 1ULL << (nl & 63);
            nl++;
        }
        if (nl != K) die("internal: %d logicals, k = %d", nl, K);
        free(v); free(isp); free(bpiv); free(bas.a);
    }
    {
        uint64_t *one = xcalloc(RA.nw, 8);
        for (int q = 0; q < n; q++) SETB(one, q);
        bm_reduce(&RA, rkA, pivA, one);
        parity = viszero(one, RA.nw);
        free(one);
    }
    Wsearch = parity ? W - (W & 1) : W;

    printf("exactB: %s %s side %c: n=%d k=%d | A (%s): m=%d rank=%d weights %d..%d, qubit degrees %d..%d | B (%s): m=%d rank=%d weights %d..%d\n",
           fhx, fhz, side, n, K, side == 'Z' ? "H_X" : "H_Z", mA, rkA, wAmin, wAmax, cmin, cmax,
           side == 'Z' ? "H_Z" : "H_X", mB, rkB, wBmin, wBmax);
    printf("parity: all-ones vector %s rowsp %s%s\n", parity ? "in" : "not in", side == 'Z' ? "H_X" : "H_Z",
           parity ? ", so every logical of this side has even weight" : "");
    if (K == 0) { printf("k = 0: no nontrivial logical exists\n"); return 0; }

    /* anchors: an ordering of the qubits by orbits; run r anchors ord[run_start[r]] and excludes ord[0 .. run_start[r]-1] */
    ord = xmalloc(n * sizeof(int));
    run_start = xmalloc((n + 1) * sizeof(int));
    if (fgens) {
        int *gens, ng = read_gens(fgens, &gens);
        for (int g = 0; g < ng; g++) {
            if (!preserves(gens + (size_t)g * n, &A, &RA, rkA, pivA) || !preserves(gens + (size_t)g * n, &B, &RB, rkB, pivB))
                die("generator %d of %s does not preserve rowsp H_X and rowsp H_Z", g + 1, fgens);
        }
        int *par = xmalloc(n * sizeof(int));
        for (int q = 0; q < n; q++) par[q] = q;
        for (int g = 0; g < ng; g++)
            for (int q = 0; q < n; q++) {
                int a = uf_find(par, q), b = uf_find(par, gens[(size_t)g * n + q]);
                if (a != b) { if (a < b) par[b] = a; else par[a] = b; }
            }
        /* orbits, each listed in increasing order; root = least element */
        int *osz = xcalloc(n, sizeof(int)), *roots = xmalloc(n * sizeof(int)), nr = 0;
        for (int q = 0; q < n; q++) { int r = uf_find(par, q); if (r == q) roots[nr++] = q; osz[r]++; }
        if (!strcmp(foorder, "size")) {
            for (int i = 1; i < nr; i++) {
                int v = roots[i], j = i - 1;
                while (j >= 0 && (osz[roots[j]] < osz[v] || (osz[roots[j]] == osz[v] && roots[j] > v))) { roots[j + 1] = roots[j]; j--; }
                roots[j + 1] = v;
            }
        } else if (strcmp(foorder, "min")) usage();
        int pos = 0;
        for (int i = 0; i < nr; i++) {
            run_start[i] = pos;
            for (int q = 0; q < n; q++) if (uf_find(par, q) == roots[i]) ord[pos++] = q;
        }
        nruns = nr;
        printf("group: %d generators from %s verified (each maps rowsp H_X and rowsp H_Z into themselves); %d orbits, sizes", ng, fgens, nr);
        for (int i = 0; i < nr && i < 20; i++) printf(" %d", osz[roots[i]]);
        if (nr > 20) printf(" ...");
        printf("\n");
        free(par); free(osz); free(roots); free(gens);
    } else {
        if (!strcmp(forder, "id")) {
            for (int q = 0; q < n; q++) ord[q] = q;
        } else if (!strcmp(forder, "bfs")) {
            char *vis = xcalloc(n, 1);
            int h = 0, t = 0;
            for (int s = 0; s < n; s++) {
                if (vis[s]) continue;
                vis[s] = 1; ord[t++] = s;
                while (h < t) {
                    int q = ord[h++];
                    for (int i = qaoff[q]; i < qaoff[q + 1]; i++) {
                        int c = qa[i];
                        for (int j = A.off[c]; j < A.off[c + 1]; j++) { int r = A.idx[j]; if (!vis[r]) { vis[r] = 1; ord[t++] = r; } }
                    }
                }
            }
            free(vis);
        } else usage();
        for (int i = 0; i < n; i++) run_start[i] = i;
        nruns = n;
        printf("anchors: every qubit (trivial group), order %s\n", forder);
    }
    if (only_run >= nruns) die("--run %d: only %d runs", only_run, nruns);

    if (opt_stabw <= 0) opt_stabw = 2 * wBmax;
    build_stabs(opt_majo, opt_stabw);
    found_set = xmalloc((n + 1) * sizeof(int));
    run_time = xcalloc(nruns, sizeof(double));
    run_done = xcalloc(nruns, sizeof(atomic_char));
    printf("search: W=%d (searched %d)%s, bound %d, majority family %d (%d members, weight <= %d), threads %d, split depth %d, runs %d%s\n",
           W, Wsearch, parity && Wsearch != W ? " by parity" : "", opt_bound, opt_majo, ns, opt_stabw, nthreads, splitD,
           nruns, only_run >= 0 ? " (one run only)" : "");
    fflush(stdout);

    double ts = wall(), cs = cputime();
    Ctx *ctx = xmalloc(nthreads * sizeof(Ctx));
    for (int t = 0; t < nthreads; t++) ctx_init(&ctx[t], t);
    atomic_store(&stopflag, 0);
    atomic_store(&next_job, 0);
    atomic_store(&total_nodes, 0);
    if (nthreads == 1) worker(&ctx[0]);
    else {
        pthread_t *th = xmalloc(nthreads * sizeof(pthread_t));
        for (int t = 0; t < nthreads; t++) pthread_create(&th[t], NULL, worker, &ctx[t]);
        for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
        free(th);
    }
    double dt = wall() - ts, dc = cputime() - cs;
    long long nodes = 0, gen = 0, leaves = 0;
    long long *rn = xcalloc(nruns, sizeof(long long));
    for (int t = 0; t < nthreads; t++) {
        nodes += ctx[t].nodes; gen += ctx[t].generated; leaves += ctx[t].leaves;
        for (int r = 0; r < nruns; r++) rn[r] += ctx[t].runnodes[r];
    }
    if (verbose || nruns <= 16) {
        for (int r = 0; r < nruns; r++) {
            if (only_run >= 0 && r != only_run) continue;
            printf("  run %d: anchor %d, excluded %d, nodes %lld", r, ord[run_start[r]], run_start[r], rn[r]);
            if (nthreads == 1) printf(", %.3f s", run_time[r]);
            if (!atomic_load(&run_done[r])) printf(" (not completed: search stopped)");
            printf("\n");
        }
    } else {
        /* the five runs with most nodes, and the share of the first tenth of the runs */
        int top[5], nt = 0;
        for (int r = 0; r < nruns; r++) {
            int j;
            if (nt < 5) j = nt++;
            else if (rn[r] > rn[top[4]]) j = 4;
            else continue;
            top[j] = r;
            while (j > 0 && rn[top[j]] > rn[top[j - 1]]) { int t = top[j]; top[j] = top[j - 1]; top[j - 1] = t; j--; }
        }
        long long s10 = 0, sall = 0;
        for (int r = 0; r < nruns; r++) { sall += rn[r]; if (r < (nruns + 9) / 10) s10 += rn[r]; }
        printf("  runs with most nodes:");
        for (int i = 0; i < nt; i++) printf(" run %d (anchor %d) %lld;", top[i], ord[run_start[top[i]]], rn[top[i]]);
        printf(" first tenth of the runs: %.1f%% of the nodes\n", sall ? 100.0 * s10 / sall : 0.0);
    }
    int st = atomic_load(&stopflag);
    printf("RESULT side=%c W=%d searched=%d found=%s nodes=%lld generated=%lld leaves=%lld wall=%.2f cpu=%.2f setup=%.2f threads=%d%s\n",
           side, W, Wsearch, found_w >= 0 ? "yes" : "none", nodes, gen, leaves, dt, dc, ts - t0, nthreads,
           st == 2 ? " NODE-LIMIT (incomplete)" : "");
    (void)c0;
    if (found_w >= 0) {
        /* recheck without the logical basis: A z = 0 and z not in rowsp B */
        int w = found_w;
        int *z = found_set;
        qsort(z, w, sizeof(int), cmp_int);
        char *inz = xcalloc(n, 1);
        for (int i = 0; i < w; i++) inz[z[i]] = 1;
        int synok = 1;
        for (int c = 0; c < mA; c++) {
            int p = 0;
            for (int i = A.off[c]; i < A.off[c + 1]; i++) p ^= inz[A.idx[i]];
            if (p) synok = 0;
        }
        uint64_t *v = xcalloc(RB.nw, 8);
        for (int i = 0; i < w; i++) SETB(v, z[i]);
        bm_reduce(&RB, rkB, pivB, v);
        int ntriv = !viszero(v, RB.nw);
        printf("WORD weight=%d syndrome=%s rowsp-test=%s :", w, synok ? "zero" : "NONZERO", ntriv ? "nontrivial" : "TRIVIAL");
        for (int i = 0; i < w; i++) printf(" %d", z[i]);
        printf("\n");
        if (!synok || !ntriv) { printf("RECHECK FAILED\n"); return 3; }
        printf("CONCLUSION: d_%c <= %d\n", side, w);
        if (fout) {
            FILE *f = fopen(fout, "w");
            if (f) { for (int i = 0; i < w; i++) fprintf(f, "%d%c", z[i], i + 1 < w ? ' ' : '\n'); fclose(f); }
        }
        free(inz); free(v);
    } else if (st == 0) {
        if (only_run >= 0) printf("CONCLUSION: run %d of %d reports nothing (not a bound by itself)\n", only_run, nruns);
        else {
            int lbd = W + 1;
            if (parity && (lbd & 1)) lbd++;
            printf("CONCLUSION: no nontrivial %c-logical of weight <= %d, d_%c >= %d%s\n", side, W, side, lbd,
                   parity ? " (parity)" : "");
        }
    }
    fflush(stdout);
    return 0;
}
