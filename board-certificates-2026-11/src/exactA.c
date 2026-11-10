/*
 * exactA: exhaustive search for nontrivial logical operators of weight <= W of a CSS code.
 *
 * Generalization of vsearch_orb.c (implementation A of exact_distances.tex, Algorithm 1 there) to
 * arbitrary and mixed check weights <= 16, any number of qubits and of logical qubits (multiword
 * classes), per-qubit or orbit anchors (Lemma 2.4 of the tex) and pattern classes for a group given
 * by generators (prompt Section 8, C1 item 2). NOTES.md maps every part of the code to the statement
 * of exact_distances.tex it implements.
 *
 *   exactA NAME.hx NAME.hz --side Z|X --W W [--group GENS] [--mode anchor|pattern] [--threads t]
 *
 * Instance files: first line "n m", then m lines of 0-based qubit indices (one check per line).
 * Side Z searches z in ker H_X outside rowsp H_Z (checks: rows of H_X, majority family and
 * nontriviality from H_Z); side X is the mirror image.
 * GENS: first the number g of generators, then g permutations of {0..n-1} (n images each); every
 * generator must map the rows of H_X onto the rows of H_X and the rows of H_Z onto the rows of H_Z.
 * Other options: --order index|bfs, --maj 0|1|2, --pack 0|1, --tie 0|1, --noparity, --split 0|1,
 * --jobs a:b, --tabmax w, --enumerate FILE, -v.
 */
#define _GNU_SOURCE
#include <pthread.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXW 16
#define TABW 12 /* largest weight with a 3^w table for restricted checks */

static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "exactA: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(2);
}
static void *xmalloc(size_t s) { void *p = malloc(s ? s : 1); if (!p) die("out of memory"); return p; }
static void *xcalloc(size_t a, size_t b) { void *p = calloc(a ? a : 1, b ? b : 1); if (!p) die("out of memory"); return p; }
static void *xrealloc(void *p, size_t s) { p = realloc(p, s ? s : 1); if (!p) die("out of memory"); return p; }
static double wall(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double cputime(void) { struct timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static int cmpint(const void *a, const void *b) { int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }

/* ------------------------------------------------------------------ sparse rows, sets of rows */
typedef struct { int n, m; int *off, *v, *raw; } Csr; /* v sorted rows; raw: input order (read_csr only) */

static Csr read_csr(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) die("cannot open %s", fn);
    char *line = NULL;
    size_t cap = 0;
    Csr M = {0};
    if (getline(&line, &cap, f) < 0 || sscanf(line, "%d %d", &M.n, &M.m) != 2 || M.n <= 0 || M.m < 0)
        die("%s: bad header", fn);
    M.off = xmalloc(sizeof(int) * (M.m + 1));
    int vcap = 1024, nv = 0, rcap = 1024;
    M.v = xmalloc(sizeof(int) * vcap);
    M.raw = xmalloc(sizeof(int) * rcap);
    M.off[0] = 0;
    for (int r = 0; r < M.m; r++) {
        if (getline(&line, &cap, f) < 0) die("%s: %d rows announced, %d found", fn, M.m, r);
        char *p = line, *e;
        int start = nv;
        for (;;) {
            long x = strtol(p, &e, 10);
            if (e == p) break;
            if (x < 0 || x >= M.n) die("%s: row %d: index %ld out of range", fn, r, x);
            if (nv == vcap) { vcap *= 2; M.v = xrealloc(M.v, sizeof(int) * vcap); }
            M.v[nv++] = (int)x;
            p = e;
        }
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p) die("%s: row %d: unexpected character '%c'", fn, r, *p);
        if (nv + 0 > rcap) { while (nv > rcap) rcap *= 2; M.raw = xrealloc(M.raw, sizeof(int) * rcap); }
        memcpy(M.raw + start, M.v + start, sizeof(int) * (nv - start));
        qsort(M.v + start, nv - start, sizeof(int), cmpint);
        for (int i = start + 1; i < nv; i++)
            if (M.v[i] == M.v[i - 1]) die("%s: row %d: repeated index %d", fn, r, M.v[i]);
        M.off[r + 1] = nv;
    }
    free(line);
    fclose(f);
    return M;
}

static uint64_t hvec(const int *a, int len) {
    uint64_t h = 0x9E3779B97F4A7C15ull ^ (uint64_t)len;
    for (int i = 0; i < len; i++) {
        h ^= (uint64_t)(uint32_t)a[i];
        h *= 0xff51afd7ed558ccdull;
        h ^= h >> 32;
    }
    return h;
}

/* growable set of sorted int vectors with lookup */
typedef struct { int m, mcap, nv, vcap, hcap; int *off, *v, *slot; } VecSet;

static void vs_init(VecSet *S) {
    memset(S, 0, sizeof *S);
    S->mcap = 64; S->vcap = 1024; S->hcap = 256;
    S->off = xmalloc(sizeof(int) * (S->mcap + 1)); S->off[0] = 0;
    S->v = xmalloc(sizeof(int) * S->vcap);
    S->slot = xmalloc(sizeof(int) * S->hcap);
    memset(S->slot, -1, sizeof(int) * S->hcap);
}
static int vs_find(const VecSet *S, const int *a, int len) {
    uint64_t h = hvec(a, len);
    for (int i = (int)(h & (S->hcap - 1));; i = (i + 1) & (S->hcap - 1)) {
        int s = S->slot[i];
        if (s < 0) return -1;
        if (S->off[s + 1] - S->off[s] == len && !memcmp(S->v + S->off[s], a, sizeof(int) * len)) return s;
    }
}
static void vs_place(VecSet *S, int s) {
    uint64_t h = hvec(S->v + S->off[s], S->off[s + 1] - S->off[s]);
    int i = (int)(h & (S->hcap - 1));
    while (S->slot[i] >= 0) i = (i + 1) & (S->hcap - 1);
    S->slot[i] = s;
}
/* returns the index; *isnew set to 1 if inserted */
static int vs_insert(VecSet *S, const int *a, int len, int *isnew) {
    int s = vs_find(S, a, len);
    if (s >= 0) { if (isnew) *isnew = 0; return s; }
    if (S->m + 1 > S->mcap) { S->mcap *= 2; S->off = xrealloc(S->off, sizeof(int) * (S->mcap + 1)); }
    while (S->nv + len > S->vcap) { S->vcap *= 2; S->v = xrealloc(S->v, sizeof(int) * S->vcap); }
    memcpy(S->v + S->nv, a, sizeof(int) * len);
    S->nv += len;
    s = S->m++;
    S->off[S->m] = S->nv;
    if (2 * S->m > S->hcap) {
        S->hcap *= 2;
        S->slot = xrealloc(S->slot, sizeof(int) * S->hcap);
        memset(S->slot, -1, sizeof(int) * S->hcap);
        for (int t = 0; t < S->m; t++) vs_place(S, t);
    } else
        vs_place(S, s);
    if (isnew) *isnew = 1;
    return s;
}

/* rows of M without empty rows and later duplicates, in the order of first occurrence */
static Csr dedupe(const Csr *M, VecSet *S, int *ndrop) {
    vs_init(S);
    *ndrop = 0;
    int *rawv = xmalloc(sizeof(int) * (M->off[M->m] + 1)), nr = 0;
    for (int r = 0; r < M->m; r++) {
        int len = M->off[r + 1] - M->off[r], isnew;
        if (len == 0) { (*ndrop)++; continue; }
        vs_insert(S, M->v + M->off[r], len, &isnew);
        if (!isnew) { (*ndrop)++; continue; }
        memcpy(rawv + nr, M->raw + M->off[r], sizeof(int) * len);
        nr += len;
    }
    Csr R = {M->n, S->m, S->off, S->v, rawv};
    return R;
}

/* ------------------------------------------------------------------ dense GF(2) echelon forms */
typedef struct { int n, nw, r, cap; uint64_t *rows; int *piv; } Ech; /* fully reduced, piv = lowest bit */

static void ech_init(Ech *E, int n, int cap) {
    E->n = n; E->nw = (n + 63) / 64; E->r = 0; E->cap = cap;
    E->rows = xcalloc((size_t)cap * E->nw, sizeof(uint64_t));
    E->piv = xmalloc(sizeof(int) * (cap + 1));
}
static inline int getbit(const uint64_t *v, int i) { return (int)(v[i >> 6] >> (i & 63) & 1); }
static void ech_reduce(const Ech *E, uint64_t *v) {
    for (int i = 0; i < E->r; i++)
        if (getbit(v, E->piv[i])) {
            const uint64_t *row = E->rows + (size_t)i * E->nw;
            for (int w = 0; w < E->nw; w++) v[w] ^= row[w];
        }
}
static int is_zero(const uint64_t *v, int nw) { for (int w = 0; w < nw; w++) if (v[w]) return 0; return 1; }
/* reduces v (destroyed); if nonzero, adds it and returns 1 */
static int ech_insert(Ech *E, uint64_t *v) {
    ech_reduce(E, v);
    int p = -1;
    for (int w = 0; w < E->nw; w++) if (v[w]) { p = 64 * w + __builtin_ctzll(v[w]); break; }
    if (p < 0) return 0;
    if (E->r >= E->cap) die("echelon capacity");
    for (int i = 0; i < E->r; i++) {
        uint64_t *row = E->rows + (size_t)i * E->nw;
        if (getbit(row, p)) for (int w = 0; w < E->nw; w++) row[w] ^= v[w];
    }
    memcpy(E->rows + (size_t)E->r * E->nw, v, sizeof(uint64_t) * E->nw);
    E->piv[E->r++] = p;
    return 1;
}
static void ech_from_csr(Ech *E, const Csr *M, int cap) {
    ech_init(E, M->n, cap);
    uint64_t *v = xmalloc(sizeof(uint64_t) * E->nw);
    for (int r = 0; r < M->m; r++) {
        memset(v, 0, sizeof(uint64_t) * E->nw);
        for (int i = M->off[r]; i < M->off[r + 1]; i++) v[M->v[i] >> 6] |= 1ull << (M->v[i] & 63);
        ech_insert(E, v);
    }
    free(v);
}

/* ------------------------------------------------------------------ global instance (read-only in the search) */
static int n, m;                  /* qubits; checks of the search side */
static int *Coff, *Cv;            /* checks, qubits in role order */
static uint16_t *cfull;           /* (1 << weight) - 1 */
static int *qoff, *qchk, qdegmax; /* qubit -> checks; largest qubit degree */
static uint8_t *qrole;            /* role of the qubit in that check */
static int ns, *soff, *sv, *sthr, *qsoff, *qsv; /* majority family S in rowsp of the dual side */
static int K, KW;
static uint64_t *lxm;             /* n x KW: bit j of qubit q = entry q of the j-th dual logical */
static int *corb;                 /* orbit index of a check (pattern mode), else -2 */
static int use_maj = 2, use_pack = 1, tie = 0, enumerate = 0, verbose = 0, keeporder = 0;
static atomic_int stopflag;
static pthread_mutex_t repmx = PTHREAD_MUTEX_INITIALIZER;
static int *word, wordlen = -1;   /* first reported word */
static long long nfound_total = 0;
static FILE *enumf = NULL;

/* anchor mode */
static int *orank, *anchorq, nanch;
/* pattern mode */
typedef struct { int w, ncls; uint16_t *rank; int *rep; long *pref; } ClassTab; /* pref[j]: masks in classes < j */
static int norb, *crep, *orbtab, *qmin;
static ClassTab *ctabs;
static int nctabs;
static uint8_t *checkfree;
static int t3[1 << TABW], p3[TABW + 1];
static uint8_t *low2;             /* lowest ternary digit equal to 2, or 255 */

typedef struct { int kind, a, b; } Job; /* kind 0: anchor a; kind 1: orbit a, class b */
typedef struct { int job, br; } Item;
static Job *jobs;
static int njobs;
static Item *items;
static long nitems;
static atomic_long nextitem;
static long long *itnodes;
static double *ittime;
static int tabmax = TABW;

/* ------------------------------------------------------------------ per-thread search state */
typedef struct {
    uint8_t *st;                 /* 0 free, 1 in T, 2 in X */
    uint16_t *inm, *outm;        /* role masks per check */
    uint8_t *unsat;
    int *ulist, *upos, nu;       /* unsatisfied checks */
    int *ucnt;                   /* number of unsatisfied checks through each qubit */
    int *scnt;
    int Tsz, *tstack;
    int *trail, tlen;
    int *cq, ncq, *sq, nsq;
    int conflict;
    uint32_t *mk, mkgen;
    int rorb, rthr;              /* restricted orbit (pattern mode) and class threshold */
    const uint16_t *rrank;
    long nforb;
    uint32_t *tab, *tabbuf;      /* 3^w table for restricted checks */
    int tab_orb, tab_thr;
    int W;
    long long nodes, leaves;
    uint64_t *cls;
} Ctx;

static void ctx_init(Ctx *x, int W) {
    memset(x, 0, sizeof *x);
    x->st = xcalloc(n, 1);
    x->inm = xcalloc(m, sizeof(uint16_t));
    x->outm = xcalloc(m, sizeof(uint16_t));
    x->unsat = xcalloc(m, 1);
    x->ucnt = xcalloc(n, sizeof(int));
    x->ulist = xmalloc(sizeof(int) * (m + 1));
    x->upos = xmalloc(sizeof(int) * (m + 1));
    x->scnt = xcalloc(ns + 1, sizeof(int));
    x->tstack = xmalloc(sizeof(int) * (n + 1));
    x->trail = xmalloc(sizeof(int) * (n + 1));
    x->cq = xmalloc(sizeof(int) * (qoff[n] + 64));
    x->sq = xmalloc(sizeof(int) * (ns + 64));
    x->mk = xcalloc(n, sizeof(uint32_t));
    x->rorb = -1;
    x->tab_orb = -1;
    x->W = W;
    x->cls = xmalloc(sizeof(uint64_t) * (KW + 1));
}

static inline void toggle_unsat(Ctx *x, int c) {
    if ((x->unsat[c] ^= 1)) {
        x->upos[c] = x->nu; x->ulist[x->nu++] = c;
        for (int i = Coff[c]; i < Coff[c + 1]; i++) x->ucnt[Cv[i]]++;
    } else {
        int p = x->upos[c], l = x->ulist[--x->nu]; x->ulist[p] = l; x->upos[l] = p;
        for (int i = Coff[c]; i < Coff[c + 1]; i++) x->ucnt[Cv[i]]--;
    }
}

static inline void put_in(Ctx *x, int q) {
    x->st[q] = 1;
    x->trail[x->tlen++] = q;
    x->tstack[x->Tsz++] = q;
    if (x->Tsz > x->W) x->conflict = 1;
    for (int i = qoff[q]; i < qoff[q + 1]; i++) {
        int c = qchk[i];
        x->inm[c] |= (uint16_t)(1u << qrole[i]);
        toggle_unsat(x, c);
        x->cq[x->ncq++] = c;
    }
    if (use_maj)
        for (int i = qsoff[q]; i < qsoff[q + 1]; i++) {
            int s = qsv[i], v = ++x->scnt[s];
            if (v > sthr[s]) x->conflict = 1;
            else if (v == sthr[s]) x->sq[x->nsq++] = s;
        }
}

static inline void put_out(Ctx *x, int q) {
    x->st[q] = 2;
    x->trail[x->tlen++] = q;
    for (int i = qoff[q]; i < qoff[q + 1]; i++) {
        int c = qchk[i];
        x->outm[c] |= (uint16_t)(1u << qrole[i]);
        x->cq[x->ncq++] = c;
    }
}

static void undo(Ctx *x, int to) {
    while (x->tlen > to) {
        int q = x->trail[--x->tlen];
        if (x->st[q] == 1) {
            x->Tsz--;
            for (int i = qoff[q]; i < qoff[q + 1]; i++) {
                int c = qchk[i];
                x->inm[c] &= (uint16_t)~(1u << qrole[i]);
                toggle_unsat(x, c);
            }
            if (use_maj) for (int i = qsoff[q]; i < qsoff[q + 1]; i++) x->scnt[qsv[i]]--;
        } else {
            for (int i = qoff[q]; i < qoff[q + 1]; i++) x->outm[qchk[i]] &= (uint16_t)~(1u << qrole[i]);
        }
        x->st[q] = 0;
    }
}

/* beta(T, X): lower bound for |E| over sets E of free qubits meeting every unsatisfied check. Propagate only
 * tests |T| + beta <= W, so both bounds stop as soon as that test is decided; the returned value then
 * exceeds the budget W - |T| exactly when beta does. */
static int lbound(Ctx *x) {
    if (x->nu == 0) return 0;
    int budget = x->W - x->Tsz;
    if (budget < 1) return budget + 1;
    int pk = 0;
    if (use_pack) {
        if (++x->mkgen == 0) { memset(x->mk, 0, sizeof(uint32_t) * n); x->mkgen = 1; }
        uint32_t g = x->mkgen;
        for (int a = 0; a < x->nu; a++) {
            int c = x->ulist[a];
            unsigned fr = cfull[c] & ~(unsigned)(x->inm[c] | x->outm[c]);
            const int *cv = Cv + Coff[c];
            int clash = 0;
            for (unsigned b = fr; b; b &= b - 1) if (x->mk[cv[__builtin_ctz(b)]] == g) { clash = 1; break; }
            if (clash) continue;
            if (++pk > budget) return pk;
            for (unsigned b = fr; b; b &= b - 1) x->mk[cv[__builtin_ctz(b)]] = g;
        }
    }
    /* ceil(nu / hmax) <= budget iff hmax >= ceil(nu / budget) */
    int hneed = (x->nu + budget - 1) / budget, hmax = 0;
    for (int a = 0; a < x->nu && hmax < hneed; a++) {
        int c = x->ulist[a];
        unsigned fr = cfull[c] & ~(unsigned)(x->inm[c] | x->outm[c]);
        const int *cv = Cv + Coff[c];
        for (unsigned b = fr; b; b &= b - 1) {
            int h = x->ucnt[cv[__builtin_ctz(b)]];
            if (h > hmax) hmax = h;
        }
    }
    if (hmax == 0) return 1 << 20;
    if (hmax < hneed) return budget + 1;
    return pk;
}

/* the rule of Propagate for a check of the restricted orbit; masks in role order of the check */
static inline int rrule(Ctx *x, unsigned in, unsigned fr, unsigned *fin, unsigned *fout) {
    if (x->tab) {
        uint32_t e = x->tab[t3[in] + 2 * t3[fr]];
        if (!(e >> 31)) return 0;
        *fin = (e & 0xFFFu) & ~in;
        *fout = fr & ~((e >> 12) & 0xFFFu);
        return 1;
    }
    int f = __builtin_popcount(fr);
    /* 2^(f-2) compatible even sets avoid (contain) any given free role; fewer are forbidden */
    if (f >= 2 && (f - 2 >= 40 || (1L << (f - 2)) > x->nforb)) { *fin = *fout = 0; return 1; }
    int par = __builtin_popcount(in) & 1, any = 0;
    unsigned inter = 0xFFFFu, uni = 0, D = fr;
    for (;;) {
        if ((__builtin_popcount(D) & 1) == par) {
            unsigned P = in | D;
            if (P == 0 || x->rrank[P] >= x->rthr) { any = 1; inter &= P; uni |= P; }
        }
        if (!D) break;
        D = (D - 1) & fr;
    }
    if (!any) return 0;
    *fin = inter & ~in;
    *fout = fr & ~uni;
    return 1;
}

static int propagate(Ctx *x) {
    while (!x->conflict && (x->nsq || x->ncq)) {
        if (x->nsq) {
            int s = x->sq[--x->nsq];
            if (x->scnt[s] > sthr[s]) { x->conflict = 1; break; }
            if (x->scnt[s] == sthr[s])
                for (int i = soff[s]; i < soff[s + 1]; i++) if (!x->st[sv[i]]) put_out(x, sv[i]);
            continue;
        }
        int c = x->cq[--x->ncq];
        unsigned in = x->inm[c], fr = cfull[c] & ~(in | x->outm[c]);
        const int *cv = Cv + Coff[c];
        if (corb[c] == x->rorb) {
            unsigned fin, fout;
            if (!rrule(x, in, fr, &fin, &fout)) { x->conflict = 1; break; }
            for (unsigned b = fin | fout; b && !x->conflict; b &= b - 1) {
                int r = __builtin_ctz(b), q = cv[r];
                if (x->st[q]) continue;
                if (fin >> r & 1) put_in(x, q); else put_out(x, q);
            }
        } else {
            int par = __builtin_popcount(in) & 1;
            if (fr == 0) { if (par) { x->conflict = 1; break; } continue; }
            if (fr & (fr - 1)) continue;
            int q = cv[__builtin_ctz(fr)];
            if (par) put_in(x, q); else put_out(x, q);
        }
    }
    int ok = !x->conflict && x->Tsz + lbound(x) <= x->W;
    x->conflict = 0;
    x->ncq = 0;
    x->nsq = 0;
    return ok;
}

static void leaf(Ctx *x) {
    x->leaves++;
    memset(x->cls, 0, sizeof(uint64_t) * KW);
    for (int i = 0; i < x->Tsz; i++) {
        const uint64_t *l = lxm + (size_t)x->tstack[i] * KW;
        for (int w = 0; w < KW; w++) x->cls[w] ^= l[w];
    }
    if (is_zero(x->cls, KW)) return;
    pthread_mutex_lock(&repmx);
    if (!enumerate && atomic_load(&stopflag)) { pthread_mutex_unlock(&repmx); return; }
    nfound_total++;
    if (wordlen < 0) {
        wordlen = x->Tsz;
        word = xmalloc(sizeof(int) * (x->Tsz + 1));
        memcpy(word, x->tstack, sizeof(int) * x->Tsz);
    }
    if (enumf) {
        int *t = xmalloc(sizeof(int) * (x->Tsz + 1));
        memcpy(t, x->tstack, sizeof(int) * x->Tsz);
        qsort(t, x->Tsz, sizeof(int), cmpint);
        fprintf(enumf, "%d :", x->Tsz);
        for (int i = 0; i < x->Tsz; i++) fprintf(enumf, " %d", t[i]);
        fprintf(enumf, "\n");
        free(t);
    }
    if (!enumerate) atomic_store(&stopflag, 1);
    pthread_mutex_unlock(&repmx);
}

/* Grow of Algorithm 1 */
static void grow(Ctx *x) {
    if (atomic_load_explicit(&stopflag, memory_order_relaxed)) return;
    x->nodes++;
    if (x->nu == 0) { leaf(x); return; }
    int best = -1, bf = 99;
    for (int a = 0; a < x->nu; a++) {
        int c = x->ulist[a];
        int f = __builtin_popcount(cfull[c] & ~(unsigned)(x->inm[c] | x->outm[c]));
        if (f < bf || (f == bf && (tie ? c > best : c < best))) { bf = f; best = c; }
    }
    int cand[MAXW], nc = 0;
    unsigned fr = cfull[best] & ~(unsigned)(x->inm[best] | x->outm[best]);
    for (unsigned b = fr; b; b &= b - 1) cand[nc++] = Cv[Coff[best] + __builtin_ctz(b)];
    int m0 = x->tlen;
    for (int i = 0; i < nc; i++) {
        if (atomic_load_explicit(&stopflag, memory_order_relaxed)) break;
        int q = cand[i];
        if (x->st[q] == 2) continue;
        if (x->st[q] == 1) { grow(x); break; }
        int m1 = x->tlen;
        put_in(x, q);
        if (propagate(x)) grow(x);
        undo(x, m1);
        if (i < nc - 1) { put_out(x, q); if (!propagate(x)) break; }
    }
    undo(x, m0);
}

/* the root call of Grow restricted to its branch br (the loop iteration i = br) */
static void root_branch(Ctx *x, int br) {
    if (br == 0) x->nodes++;
    if (x->nu == 0) { if (br == 0) leaf(x); return; }
    int best = -1, bf = 99;
    for (int a = 0; a < x->nu; a++) {
        int c = x->ulist[a];
        int f = __builtin_popcount(cfull[c] & ~(unsigned)(x->inm[c] | x->outm[c]));
        if (f < bf || (f == bf && (tie ? c > best : c < best))) { bf = f; best = c; }
    }
    int cand[MAXW], nc = 0;
    unsigned fr = cfull[best] & ~(unsigned)(x->inm[best] | x->outm[best]);
    for (unsigned b = fr; b; b &= b - 1) cand[nc++] = Cv[Coff[best] + __builtin_ctz(b)];
    if (br >= nc) return;
    for (int i = 0; i <= br; i++) {
        int q = cand[i];
        if (x->st[q] == 2) { if (i == br) return; continue; }
        if (x->st[q] == 1) { if (i == br) grow(x); return; }
        if (i == br) {
            int m1 = x->tlen;
            put_in(x, q);
            if (propagate(x)) grow(x);
            undo(x, m1);
            return;
        }
        put_out(x, q);
        if (!propagate(x)) return;
    }
}

/* ------------------------------------------------------------------ jobs */
static void build_table(Ctx *x, int orb, int thr) {
    const ClassTab *T = &ctabs[orbtab[orb]];
    int w = T->w, N3 = p3[w];
    uint32_t *tab = x->tabbuf;
    for (unsigned P = 0; P < (1u << w); P++) {
        int ok = !(__builtin_popcount(P) & 1) && (P == 0 || T->rank[P] >= thr);
        tab[t3[P]] = ok ? (0x80000000u | P | (P << 12)) : 0;
    }
    for (int t = 0; t < N3; t++) {
        int j = low2[t];
        if (j == 255) continue;
        uint32_t a = tab[t - 2 * p3[j]], b = tab[t - p3[j]];
        if (!(a >> 31)) tab[t] = b;
        else if (!(b >> 31)) tab[t] = a;
        else tab[t] = 0x80000000u | ((a & b) & 0xFFFu) | ((a | b) & (0xFFFu << 12));
    }
    x->tab = tab;
    x->tab_orb = orb;
    x->tab_thr = thr;
}

/* applies T_0 and X_0 of the job; returns the result of Propagate */
static int setup_job(Ctx *x, const Job *J) {
    if (J->kind == 0) {
        for (int q = 0; q < n; q++) if (orank[q] < J->a) put_out(x, q);
        put_in(x, anchorq[J->a]);
        return propagate(x);
    }
    int i = J->a, j = J->b;
    const ClassTab *T = &ctabs[orbtab[i]];
    x->rorb = i;
    x->rthr = j;
    x->rrank = T->rank;
    x->nforb = T->pref[j];
    x->tab = NULL;
    if (T->w <= tabmax) {
        if (!(x->tab_orb == i && x->tab_thr == j)) build_table(x, i, j);
        else x->tab = x->tabbuf;
    }
    for (int q = 0; q < n; q++) if (qmin[q] < i || checkfree[q]) put_out(x, q);
    int c = crep[i];
    unsigned pi = (unsigned)T->rep[j];
    for (int r = 0; r < T->w; r++) {
        int q = Cv[Coff[c] + r];
        if (pi >> r & 1) {
            if (x->st[q]) { x->ncq = x->nsq = 0; x->conflict = 0; return 0; }
            put_in(x, q);
        } else if (!x->st[q]) put_out(x, q);
    }
    return propagate(x);
}

static void *worker(void *arg) {
    Ctx *x = arg;
    for (;;) {
        if (atomic_load(&stopflag)) break;
        long it = atomic_fetch_add(&nextitem, 1);
        if (it >= nitems) break;
        const Job *J = &jobs[items[it].job];
        long long n0 = x->nodes;
        double t0 = wall();
        if (setup_job(x, J)) {
            if (items[it].br < 0) grow(x);
            else root_branch(x, items[it].br);
        }
        undo(x, 0);
        x->rorb = -1;
        x->tab = NULL;
        itnodes[it] = x->nodes - n0;
        ittime[it] = wall() - t0;
        if (verbose)
            fprintf(stderr, "item %ld job %d (%s %d %d) branch %d: nodes %lld %.2f s\n", it, items[it].job,
                    J->kind ? "pattern" : "anchor", J->a, J->b, items[it].br, itnodes[it], ittime[it]);
    }
    return NULL;
}

/* ------------------------------------------------------------------ group, orbits, classes */
static int ngen;
static int *gq;   /* ngen x n */
static int *gcc;  /* ngen x m: action on the checks of the search side */

static void read_group(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) die("cannot open %s", fn);
    if (fscanf(f, "%d", &ngen) != 1 || ngen < 0) die("%s: bad header", fn);
    gq = xmalloc(sizeof(int) * ((size_t)ngen * n + 1));
    uint8_t *seen = xmalloc(n);
    for (int g = 0; g < ngen; g++) {
        memset(seen, 0, n);
        for (int q = 0; q < n; q++) {
            int v;
            if (fscanf(f, "%d", &v) != 1 || v < 0 || v >= n || seen[v]) die("%s: generator %d is not a permutation of 0..%d", fn, g, n - 1);
            seen[v] = 1;
            gq[(size_t)g * n + q] = v;
        }
    }
    int extra;
    if (fscanf(f, "%d", &extra) == 1) die("%s: more numbers than %d permutations of %d points", fn, ngen, n);
    free(seen);
    fclose(f);
}

/* image of every row of R under generator g, as an index of R; dies if some image is not a row */
static void act_on_rows(int g, const Csr *R, const VecSet *S, int *img, const char *what) {
    int *buf = xmalloc(sizeof(int) * (n + 1));
    for (int r = 0; r < R->m; r++) {
        int len = R->off[r + 1] - R->off[r];
        for (int i = 0; i < len; i++) buf[i] = gq[(size_t)g * n + R->v[R->off[r] + i]];
        qsort(buf, len, sizeof(int), cmpint);
        int s = vs_find(S, buf, len);
        if (s < 0) die("generator %d does not map the rows of %s onto themselves (row %d)", g, what, r);
        if (img) img[r] = s;
    }
    free(buf);
}

static int cmp_heavy_first(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b, px = __builtin_popcount(x), py = __builtin_popcount(y);
    if (px != py) return py - px;
    return (x > y) - (x < y);
}
static const int *sortkey;
static int cmp_by_key(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    if (sortkey[x] != sortkey[y]) return (sortkey[x] > sortkey[y]) - (sortkey[x] < sortkey[y]);
    return (x > y) - (x < y);
}

static int find_(int *p, int a) { while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; } return a; }

/* classes of nonempty even role masks under the role permutations perms (np x w) */
static void make_classes(ClassTab *T, int w, const uint8_t *perms, int np) {
    int M = 1 << w;
    int *par = xmalloc(sizeof(int) * M);
    for (int i = 0; i < M; i++) par[i] = i;
    unsigned lo[256], hi[256];
    for (int p = 0; p < np; p++) {
        const uint8_t *pm = perms + (size_t)p * w;
        for (int b = 0; b < 256; b++) {
            unsigned a = 0, h = 0;
            for (int r = 0; r < 8; r++) {
                if (b >> r & 1) {
                    if (r < w) a |= 1u << pm[r];
                    if (r + 8 < w) h |= 1u << pm[r + 8];
                }
            }
            lo[b] = a; hi[b] = h;
        }
        for (int P = 1; P < M; P++) {
            if (__builtin_popcount(P) & 1) continue;
            int Q = (int)(lo[P & 255] | hi[P >> 8]);
            int a = find_(par, P), b = find_(par, Q);
            if (a != b) par[a < b ? b : a] = a < b ? a : b;
        }
    }
    /* roots are minimal elements; classes sorted by decreasing weight, then by representative */
    int *reps = xmalloc(sizeof(int) * M), nr = 0;
    for (int P = 1; P < M; P++) if (!(__builtin_popcount(P) & 1) && find_(par, P) == P) reps[nr++] = P;
    qsort(reps, nr, sizeof(int), cmp_heavy_first);
    int *idx = xmalloc(sizeof(int) * M);
    for (int a = 0; a < nr; a++) idx[reps[a]] = a;
    T->w = w;
    T->ncls = nr;
    T->rep = reps;
    T->rank = xmalloc(sizeof(uint16_t) * M);
    T->pref = xcalloc(nr + 1, sizeof(long));
    long *cnt = xcalloc(nr + 1, sizeof(long));
    for (int P = 0; P < M; P++) {
        if (P == 0 || (__builtin_popcount(P) & 1)) { T->rank[P] = 0xFFFF; continue; }
        int a = idx[find_(par, P)];
        T->rank[P] = (uint16_t)a;
        cnt[a]++;
    }
    for (int a = 0; a < nr; a++) T->pref[a + 1] = T->pref[a] + cnt[a];
    free(cnt);
    free(idx);
    free(par);
}

/* ------------------------------------------------------------------ main */
static void usage(void) {
    die("usage: exactA NAME.hx NAME.hz --side Z|X --W W [--group GENS] [--mode anchor|pattern] [--threads t]\n"
        "  [--order index|bfs] [--maj 0|1|2] [--pack 0|1] [--tie 0|1] [--noparity] [--split 0|1] [--jobs a:b]\n"
        "  [--tabmax w] [--enumerate FILE] [--keeporder] [-v]");
}

int main(int argc, char **argv) {
    if (argc < 3) usage();
    const char *fx = argv[1], *fz = argv[2], *gfile = NULL, *efile = NULL;
    char side = 0;
    int W = -1, threads = 1, mode = -1, order = 0, noparity = 0, split = -1, ja = 0, jb = -1;
    for (int i = 3; i < argc; i++) {
        const char *a = argv[i];
        const char *v = i + 1 < argc ? argv[i + 1] : NULL;
        if (!strcmp(a, "-v")) { verbose = 1; continue; }
        if (!strcmp(a, "--noparity")) { noparity = 1; continue; }
        if (!strcmp(a, "--keeporder")) { keeporder = 1; continue; }
        if (!v) usage();
        i++;
        if (!strcmp(a, "--side")) { if (strcmp(v, "Z") && strcmp(v, "X")) usage(); side = v[0]; }
        else if (!strcmp(a, "--W")) W = atoi(v);
        else if (!strcmp(a, "--group")) gfile = v;
        else if (!strcmp(a, "--mode")) { if (!strcmp(v, "anchor")) mode = 0; else if (!strcmp(v, "pattern")) mode = 1; else usage(); }
        else if (!strcmp(a, "--threads")) threads = atoi(v);
        else if (!strcmp(a, "--order")) { if (!strcmp(v, "index")) order = 0; else if (!strcmp(v, "bfs")) order = 1; else usage(); }
        else if (!strcmp(a, "--maj")) use_maj = atoi(v);
        else if (!strcmp(a, "--pack")) use_pack = atoi(v);
        else if (!strcmp(a, "--tie")) tie = atoi(v);
        else if (!strcmp(a, "--split")) split = atoi(v);
        else if (!strcmp(a, "--tabmax")) tabmax = atoi(v);
        else if (!strcmp(a, "--enumerate")) { efile = v; enumerate = 1; }
        else if (!strcmp(a, "--jobs")) { if (sscanf(v, "%d:%d", &ja, &jb) != 2) usage(); }
        else usage();
    }
    if (!side || W < 0 || threads < 1 || use_maj < 0 || use_maj > 2 || tabmax > TABW) usage();
    if (mode < 0) mode = 0;
    double tw0 = wall(), tc0 = cputime();

    /* instance */
    Csr HXr = read_csr(fx), HZr = read_csr(fz);
    if (HXr.n != HZr.n) die("H_X has %d columns, H_Z %d", HXr.n, HZr.n);
    n = HXr.n;
    VecSet SX, SZ;
    int dropx, dropz;
    Csr HX = dedupe(&HXr, &SX, &dropx), HZ = dedupe(&HZr, &SZ, &dropz);
    /* CSS condition: every row of H_Z meets every row of H_X evenly */
    {
        int *xo = xcalloc(n + 1, sizeof(int)), *xl;
        for (int r = 0; r < HX.m; r++) for (int i = HX.off[r]; i < HX.off[r + 1]; i++) xo[HX.v[i] + 1]++;
        for (int q = 0; q < n; q++) xo[q + 1] += xo[q];
        xl = xmalloc(sizeof(int) * (HX.off[HX.m] + 1));
        int *pos = xmalloc(sizeof(int) * (n + 1));
        memcpy(pos, xo, sizeof(int) * (n + 1));
        for (int r = 0; r < HX.m; r++) for (int i = HX.off[r]; i < HX.off[r + 1]; i++) xl[pos[HX.v[i]]++] = r;
        int *cnt = xcalloc(HX.m + 1, sizeof(int)), *touched = xmalloc(sizeof(int) * (HX.off[HX.m] + 1));
        for (int r = 0; r < HZ.m; r++) {
            int nt = 0;
            for (int i = HZ.off[r]; i < HZ.off[r + 1]; i++)
                for (int k = xo[HZ.v[i]]; k < xo[HZ.v[i] + 1]; k++) { if (!cnt[xl[k]]) touched[nt++] = xl[k]; cnt[xl[k]]++; }
            for (int t = 0; t < nt; t++) { if (cnt[touched[t]] & 1) die("not CSS: row %d of H_Z meets row %d of H_X oddly", r, touched[t]); cnt[touched[t]] = 0; }
        }
        free(xo); free(xl); free(pos); free(cnt); free(touched);
    }
    Csr *C = side == 'Z' ? &HX : &HZ, *D = side == 'Z' ? &HZ : &HX;
    VecSet *SC = side == 'Z' ? &SX : &SZ;
    const char *cname = side == 'Z' ? "H_X" : "H_Z", *dname = side == 'Z' ? "H_Z" : "H_X";
    m = C->m;

    /* linear algebra: ranks, k, dual logical basis (ker D modulo rowsp C), parity */
    Ech EC, ED;
    ech_from_csr(&EC, C, n);
    ech_from_csr(&ED, D, n);
    int rC = EC.r, rD = ED.r, nw = EC.nw;
    int kexp = n - rC - rD;
    if (kexp < 0) die("negative k");
    int *ispiv = xcalloc(n, sizeof(int));
    for (int i = 0; i < rD; i++) ispiv[ED.piv[i]] = 1;
    Ech EX;
    ech_init(&EX, n, n);
    memcpy(EX.rows, EC.rows, sizeof(uint64_t) * (size_t)rC * nw);
    memcpy(EX.piv, EC.piv, sizeof(int) * rC);
    EX.r = rC;
    K = 0;
    KW = kexp > 0 ? (kexp + 63) / 64 : 1;
    lxm = xcalloc((size_t)n * KW, sizeof(uint64_t));
    uint64_t *v = xmalloc(sizeof(uint64_t) * nw), *tmp = xmalloc(sizeof(uint64_t) * nw);
    for (int f = 0; f < n && K < kexp; f++) {
        if (ispiv[f]) continue;
        memset(v, 0, sizeof(uint64_t) * nw);
        v[f >> 6] |= 1ull << (f & 63);
        for (int i = 0; i < rD; i++) if (getbit(ED.rows + (size_t)i * nw, f)) v[ED.piv[i] >> 6] |= 1ull << (ED.piv[i] & 63);
        memcpy(tmp, v, sizeof(uint64_t) * nw);
        if (ech_insert(&EX, tmp)) {
            for (int q = 0; q < n; q++) if (getbit(v, q)) lxm[(size_t)q * KW + K / 64] |= 1ull << (K % 64);
            K++;
        }
    }
    if (K != kexp) die("logical basis: found %d classes, expected %d", K, kexp);
    free(EX.rows); free(EX.piv); free(ispiv);
    for (int w = 0; w < nw; w++) v[w] = ~0ull;
    if (n % 64) v[nw - 1] = (1ull << (n % 64)) - 1;
    ech_reduce(&EC, v);
    int parity = is_zero(v, nw);
    free(EC.rows); free(EC.piv);

    /* checks, weights */
    int wmax = 0, wcount[MAXW + 1] = {0}, ncheckfree = 0;
    for (int c = 0; c < m; c++) {
        int w = C->off[c + 1] - C->off[c];
        if (w > MAXW) die("check %d of %s has weight %d > %d", c, cname, w, MAXW);
        wcount[w]++;
        if (w > wmax) wmax = w;
    }
    Coff = xmalloc(sizeof(int) * (m + 1));
    Cv = xmalloc(sizeof(int) * (C->off[m] + 1));
    memcpy(Coff, C->off, sizeof(int) * (m + 1));
    memcpy(Cv, keeporder ? C->raw : C->v, sizeof(int) * C->off[m]);
    cfull = xmalloc(sizeof(uint16_t) * (m + 1));
    for (int c = 0; c < m; c++) cfull[c] = (uint16_t)((1u << (Coff[c + 1] - Coff[c])) - 1);
    checkfree = xcalloc(n, 1);
    {
        int *deg = xcalloc(n, sizeof(int));
        for (int i = 0; i < Coff[m]; i++) deg[Cv[i]]++;
        for (int q = 0; q < n; q++) if (!deg[q]) { checkfree[q] = 1; ncheckfree++; }
        free(deg);
    }

    /* majority family: rows of D and sums of two overlapping rows of D */
    VecSet S;
    vs_init(&S);
    if (use_maj >= 1) for (int r = 0; r < D->m; r++) vs_insert(&S, D->v + D->off[r], D->off[r + 1] - D->off[r], NULL);
    if (use_maj >= 2) {
        int *dofs = xcalloc(n + 1, sizeof(int));
        for (int i = 0; i < D->off[D->m]; i++) dofs[D->v[i] + 1]++;
        for (int q = 0; q < n; q++) dofs[q + 1] += dofs[q];
        int *dl = xmalloc(sizeof(int) * (D->off[D->m] + 1)), *pos = xmalloc(sizeof(int) * (n + 1));
        memcpy(pos, dofs, sizeof(int) * (n + 1));
        for (int r = 0; r < D->m; r++) for (int i = D->off[r]; i < D->off[r + 1]; i++) dl[pos[D->v[i]]++] = r;
        int *mark = xmalloc(sizeof(int) * (D->m + 1)), *buf = xmalloc(sizeof(int) * (2 * n + 2));
        int *nb = xmalloc(sizeof(int) * (D->m + 1));
        for (int r = 0; r < D->m; r++) mark[r] = -1;
        for (int r = 0; r < D->m; r++) {
            int nnb = 0;
            for (int i = D->off[r]; i < D->off[r + 1]; i++)
                for (int k = dofs[D->v[i]]; k < dofs[D->v[i] + 1]; k++) {
                    int s = dl[k];
                    if (s <= r || mark[s] == r) continue;
                    mark[s] = r;
                    nb[nnb++] = s;
                }
            /* pairs (r, s) in increasing s, as in codes.py */
            qsort(nb, nnb, sizeof(int), cmpint);
            for (int t = 0; t < nnb; t++) {
                int s = nb[t];
                {
                    int a = D->off[r], ae = D->off[r + 1], b = D->off[s], be = D->off[s + 1], len = 0;
                    while (a < ae || b < be) {
                        if (b >= be || (a < ae && D->v[a] < D->v[b])) buf[len++] = D->v[a++];
                        else if (a >= ae || D->v[b] < D->v[a]) buf[len++] = D->v[b++];
                        else { a++; b++; }
                    }
                    if (len) vs_insert(&S, buf, len, NULL);
                }
            }
        }
        free(dofs); free(dl); free(pos); free(mark); free(buf); free(nb);
    }
    ns = S.m;
    soff = S.off;
    sv = S.v;
    sthr = xmalloc(sizeof(int) * (ns + 1));
    for (int s = 0; s < ns; s++) sthr[s] = (soff[s + 1] - soff[s]) / 2;
    qsoff = xcalloc(n + 1, sizeof(int));
    for (int i = 0; i < soff[ns]; i++) qsoff[sv[i] + 1]++;
    for (int q = 0; q < n; q++) qsoff[q + 1] += qsoff[q];
    qsv = xmalloc(sizeof(int) * (soff[ns] + 1));
    {
        int *pos = xmalloc(sizeof(int) * (n + 1));
        memcpy(pos, qsoff, sizeof(int) * (n + 1));
        for (int s = 0; s < ns; s++) for (int i = soff[s]; i < soff[s + 1]; i++) qsv[pos[sv[i]]++] = s;
        free(pos);
    }
    if (!use_maj) ns = 0;

    /* group */
    if (gfile) {
        read_group(gfile);
        gcc = xmalloc(sizeof(int) * ((size_t)ngen * m + 1));
        for (int g = 0; g < ngen; g++) {
            act_on_rows(g, C, SC, gcc + (size_t)g * m, cname);
            act_on_rows(g, D, side == 'Z' ? &SZ : &SX, NULL, dname);
        }
    } else
        ngen = 0;

    /* BFS order of the qubits along the checks, from qubit 0 (used by --order bfs) */
    int *bfsrank = xmalloc(sizeof(int) * n);
    {
        int *cof = xcalloc(n + 1, sizeof(int));
        for (int i = 0; i < Coff[m]; i++) cof[Cv[i] + 1]++;
        for (int q = 0; q < n; q++) cof[q + 1] += cof[q];
        int *cl = xmalloc(sizeof(int) * (Coff[m] + 1)), *pos = xmalloc(sizeof(int) * (n + 1));
        memcpy(pos, cof, sizeof(int) * (n + 1));
        for (int c = 0; c < m; c++) for (int i = Coff[c]; i < Coff[c + 1]; i++) cl[pos[Cv[i]]++] = c;
        int *queue = xmalloc(sizeof(int) * n), head = 0, tail = 0, r = 0;
        uint8_t *cseen = xcalloc(m + 1, 1);
        for (int q = 0; q < n; q++) bfsrank[q] = -1;
        for (int s = 0; s < n; s++) {
            if (bfsrank[s] >= 0) continue;
            bfsrank[s] = r++;
            queue[tail++] = s;
            while (head < tail) {
                int q = queue[head++];
                for (int k = cof[q]; k < cof[q + 1]; k++) {
                    int c = cl[k];
                    if (cseen[c]) continue;
                    cseen[c] = 1;
                    for (int i = Coff[c]; i < Coff[c + 1]; i++) if (bfsrank[Cv[i]] < 0) { bfsrank[Cv[i]] = r++; queue[tail++] = Cv[i]; }
                }
            }
        }
        free(cof); free(cl); free(pos); free(queue); free(cseen);
    }

    corb = xmalloc(sizeof(int) * (m + 1));
    for (int c = 0; c < m; c++) corb[c] = -2;
    int nqorb = 0;
    if (mode == 0) {
        /* Lemma 2.4: orbits O_1, ..., O_r of the group on the qubits, anchor q_i in O_i, X_i = O_1 u ... u O_{i-1} */
        int *par = xmalloc(sizeof(int) * n);
        for (int q = 0; q < n; q++) par[q] = q;
        for (int g = 0; g < ngen; g++)
            for (int q = 0; q < n; q++) {
                int a = find_(par, q), b = find_(par, gq[(size_t)g * n + q]);
                if (a != b) par[a < b ? b : a] = a < b ? a : b;
            }
        int *key = xmalloc(sizeof(int) * n); /* orbit key: least index (order index) or least BFS rank */
        for (int q = 0; q < n; q++) key[q] = -1;
        for (int q = 0; q < n; q++) {
            int rt = find_(par, q), kq = order ? bfsrank[q] : q;
            if (key[rt] < 0 || kq < key[rt]) key[rt] = kq;
        }
        int *roots = xmalloc(sizeof(int) * n);
        for (int q = 0; q < n; q++) if (find_(par, q) == q) roots[nqorb++] = q;
        sortkey = key;
        qsort(roots, nqorb, sizeof(int), cmp_by_key);
        int *rk = xmalloc(sizeof(int) * n);
        for (int a = 0; a < nqorb; a++) rk[roots[a]] = a;
        orank = xmalloc(sizeof(int) * n);
        anchorq = xmalloc(sizeof(int) * nqorb);
        for (int a = 0; a < nqorb; a++) anchorq[a] = -1;
        for (int q = 0; q < n; q++) {
            int a = rk[find_(par, q)];
            orank[q] = a;
            int kq = order ? bfsrank[q] : q;
            if (anchorq[a] < 0 || kq < (order ? bfsrank[anchorq[a]] : anchorq[a])) anchorq[a] = q;
        }
        nanch = nqorb;
        njobs = nqorb;
        jobs = xmalloc(sizeof(Job) * njobs);
        for (int a = 0; a < njobs; a++) { jobs[a].kind = 0; jobs[a].a = a; jobs[a].b = 0; }
        free(par); free(key); free(roots); free(rk);
    } else {
        /* C1 item 2: orbits c_1, ..., c_r on the checks, role orders rho_c = gamma_c o rho_{c_i}, classes of S_i */
        int *rho = xmalloc(sizeof(int) * (Coff[m] + 1));
        int *queue = xmalloc(sizeof(int) * (m + 1));
        crep = xmalloc(sizeof(int) * (m + 1));
        norb = 0;
        int *corder = xmalloc(sizeof(int) * (m + 1));
        for (int c = 0; c < m; c++) corder[c] = c;
        if (order) {
            /* checks by least BFS rank of their qubits */
            int *ck = xmalloc(sizeof(int) * (m + 1));
            for (int c = 0; c < m; c++) { ck[c] = n; for (int i = Coff[c]; i < Coff[c + 1]; i++) if (bfsrank[Cv[i]] < ck[c]) ck[c] = bfsrank[Cv[i]]; }
            sortkey = ck;
            qsort(corder, m, sizeof(int), cmp_by_key);
            free(ck);
        }
        int *pos = xmalloc(sizeof(int) * n);
        for (int q = 0; q < n; q++) pos[q] = -1;
        uint8_t *perms = NULL;
        size_t pcap = 0;
        orbtab = xmalloc(sizeof(int) * (m + 1));
        ctabs = xmalloc(sizeof(ClassTab) * (m + 1));
        nctabs = 0;
        for (int a = 0; a < m; a++) {
            int c0 = corder[a];
            if (corb[c0] >= 0) continue;
            int i = norb++, w = Coff[c0 + 1] - Coff[c0];
            crep[i] = c0;
            corb[c0] = i;
            memcpy(rho + Coff[c0], Cv + Coff[c0], sizeof(int) * w);
            int head = 0, tail = 0;
            queue[tail++] = c0;
            while (head < tail) {
                int c = queue[head++];
                for (int g = 0; g < ngen; g++) {
                    int c2 = gcc[(size_t)g * m + c];
                    if (corb[c2] >= 0) continue;
                    corb[c2] = i;
                    for (int r = 0; r < w; r++) rho[Coff[c2] + r] = gq[(size_t)g * n + rho[Coff[c] + r]];
                    queue[tail++] = c2;
                }
            }
            /* Schreier generators gamma_{gc}^{-1} g gamma_c restricted to c_i, as role permutations */
            int np = 0;
            for (int h = 0; h < tail; h++) {
                int c = queue[h];
                for (int g = 0; g < ngen; g++) {
                    int c2 = gcc[(size_t)g * m + c];
                    for (int r = 0; r < w; r++) pos[rho[Coff[c2] + r]] = r;
                    if ((size_t)(np + 1) * w > pcap) { pcap = 2 * (size_t)(np + 1) * w + 64; perms = xrealloc(perms, pcap); }
                    int ident = 1;
                    for (int r = 0; r < w; r++) {
                        int r2 = pos[gq[(size_t)g * n + rho[Coff[c] + r]]];
                        if (r2 < 0) die("internal: Schreier generator");
                        perms[(size_t)np * w + r] = (uint8_t)r2;
                        if (r2 != r) ident = 0;
                    }
                    for (int r = 0; r < w; r++) pos[rho[Coff[c2] + r]] = -1;
                    if (ident) continue;
                    int dup = 0;
                    for (int p = 0; p < np && !dup; p++) if (!memcmp(perms + (size_t)p * w, perms + (size_t)np * w, w)) dup = 1;
                    if (!dup) np++;
                }
            }
            /* reuse a class table of the same weight when the stabilizer acts trivially */
            int t = -1;
            if (np == 0)
                for (int s = 0; s < nctabs; s++) if (ctabs[s].w == w && ctabs[s].ncls == (1 << (w - 1)) - 1) { t = s; break; }
            if (t < 0) { make_classes(&ctabs[nctabs], w, perms, np); t = nctabs++; }
            orbtab[i] = t;
        }
        memcpy(Cv, rho, sizeof(int) * Coff[m]);
        free(rho); free(queue); free(pos); free(perms); free(corder);
        qmin = xmalloc(sizeof(int) * n);
        for (int q = 0; q < n; q++) qmin[q] = 1 << 30;
        for (int c = 0; c < m; c++) for (int i = Coff[c]; i < Coff[c + 1]; i++) if (corb[c] < qmin[Cv[i]]) qmin[Cv[i]] = corb[c];
        njobs = 0;
        for (int i = 0; i < norb; i++) njobs += ctabs[orbtab[i]].ncls;
        jobs = xmalloc(sizeof(Job) * (njobs + 1));
        int jj = 0;
        for (int i = 0; i < norb; i++)
            for (int j = 0; j < ctabs[orbtab[i]].ncls; j++) { jobs[jj].kind = 1; jobs[jj].a = i; jobs[jj].b = j; jj++; }
        /* ternary indexing for the tables */
        p3[0] = 1;
        for (int j = 1; j <= TABW; j++) p3[j] = 3 * p3[j - 1];
        for (int P = 0; P < (1 << TABW); P++) { t3[P] = 0; for (int r = 0; r < TABW; r++) if (P >> r & 1) t3[P] += p3[r]; }
        int wt = 0;
        for (int s = 0; s < nctabs; s++) if (ctabs[s].w <= tabmax && ctabs[s].w > wt) wt = ctabs[s].w;
        if (wt > 0) {
            low2 = xmalloc(p3[wt]);
            for (int t = 0; t < p3[wt]; t++) {
                int u = t, j = 0;
                low2[t] = 255;
                while (u) { if (u % 3 == 2) { low2[t] = (uint8_t)j; break; } u /= 3; j++; }
            }
        }
    }

    /* qubit -> (check, role), after the role orders are fixed */
    qoff = xcalloc(n + 1, sizeof(int));
    for (int i = 0; i < Coff[m]; i++) qoff[Cv[i] + 1]++;
    for (int q = 0; q < n; q++) qoff[q + 1] += qoff[q];
    qchk = xmalloc(sizeof(int) * (Coff[m] + 1));
    qrole = xmalloc(Coff[m] + 1);
    qdegmax = 0;
    for (int q = 0; q < n; q++) if (qoff[q + 1] - qoff[q] > qdegmax) qdegmax = qoff[q + 1] - qoff[q];
    {
        int *pos = xmalloc(sizeof(int) * (n + 1));
        memcpy(pos, qoff, sizeof(int) * (n + 1));
        for (int c = 0; c < m; c++)
            for (int i = Coff[c]; i < Coff[c + 1]; i++) { int q = Cv[i]; qchk[pos[q]] = c; qrole[pos[q]] = (uint8_t)(i - Coff[c]); pos[q]++; }
        free(pos);
    }

    int Weff = W;
    if (parity && !noparity && (W & 1)) Weff = W - 1;

    printf("exactA: side %c, n=%d, checks %s: %d rows (%d dropped: empty or repeated), dual %s: %d rows (%d dropped)\n",
           side, n, cname, m, side == 'Z' ? dropx : dropz, dname, D->m, side == 'Z' ? dropz : dropx);
    printf("  check weights:");
    for (int w = 0; w <= MAXW; w++) if (wcount[w]) printf(" %d^%d", w, wcount[w]);
    printf("; qubits in no check: %d\n", ncheckfree);
    printf("  rank %s=%d, rank %s=%d, k=%d (%d words per class), majority family %d sets (level %d)\n", cname, rC, dname, rD, K, KW, ns, use_maj);
    printf("  parity: all-ones vector %s rowsp %s, so %s\n", parity ? "in" : "not in", cname,
           parity ? "every logical of this side has even weight" : "no parity reduction");
    if (gfile) printf("  group: %d generators from %s, preserving the rows of H_X and H_Z\n", ngen, gfile);
    if (mode == 0) printf("  mode anchor (Lemma 2.4): %d qubit orbits, order %s\n", nqorb, order ? "bfs" : "index");
    else {
        printf("  mode pattern (C1 item 2): %d check orbits, %d jobs, order %s; classes per orbit:", norb, njobs, order ? "bfs" : "index");
        for (int i = 0; i < norb && i < 8; i++) printf(" %d", ctabs[orbtab[i]].ncls);
        if (norb > 8) printf(" ...");
        printf("\n");
    }
    printf("  W=%d, searching weights <= %d\n", W, Weff);
    fflush(stdout);

    /* weight-one words on qubits in no check (Lemma 2.1(a): a minimal word through such a qubit is that qubit) */
    if (K > 0 && W >= 1)
        for (int q = 0; q < n && wordlen < 0; q++)
            if (checkfree[q] && !is_zero(lxm + (size_t)q * KW, KW)) { word = xmalloc(sizeof(int)); word[0] = q; wordlen = 1; nfound_total = 1; atomic_store(&stopflag, 1); }

    if (jb < 0 || jb > njobs) jb = njobs;
    if (ja < 0) ja = 0;
    if (ja > jb) ja = jb;
    if (split < 0) split = (threads > 1 && (jb - ja) < 4 * threads) ? 1 : 0;
    nitems = (long)(jb - ja) * (split ? wmax : 1);
    items = xmalloc(sizeof(Item) * (nitems + 1));
    {
        long t = 0;
        for (int j = ja; j < jb; j++) {
            if (split) for (int b = 0; b < wmax; b++) { items[t].job = j; items[t].br = b; t++; }
            else { items[t].job = j; items[t].br = -1; t++; }
        }
    }
    itnodes = xcalloc(nitems + 1, sizeof(long long));
    ittime = xcalloc(nitems + 1, sizeof(double));
    if (efile) { enumf = fopen(efile, "w"); if (!enumf) die("cannot open %s", efile); }
    atomic_store(&nextitem, 0);

    long long nodes = 0, leaves = 0;
    double tws = wall();
    if (K > 0 && wordlen < 0) {
        Ctx *ctx = xmalloc(sizeof(Ctx) * threads);
        pthread_t *th = xmalloc(sizeof(pthread_t) * threads);
        for (int t = 0; t < threads; t++) {
            ctx_init(&ctx[t], Weff);
            if (mode == 1) ctx[t].tabbuf = xmalloc(sizeof(uint32_t) * (size_t)p3[TABW < wmax ? TABW : wmax]);
        }
        if (threads == 1) worker(&ctx[0]);
        else {
            for (int t = 0; t < threads; t++) pthread_create(&th[t], NULL, worker, &ctx[t]);
            for (int t = 0; t < threads; t++) pthread_join(th[t], NULL);
        }
        for (int t = 0; t < threads; t++) { nodes += ctx[t].nodes; leaves += ctx[t].leaves; }
    }
    double tw = wall() - tw0, tc = cputime() - tc0, tsearch = wall() - tws;
    if (enumf) fclose(enumf);

    int complete = (ja == 0 && jb == njobs);
    if (verbose && nitems > 0) {
        /* the most expensive items */
        long best[5] = {-1, -1, -1, -1, -1};
        for (long it = 0; it < nitems; it++)
            for (int a = 0; a < 5; a++)
                if (best[a] < 0 || itnodes[it] > itnodes[best[a]]) { for (int b = 4; b > a; b--) best[b] = best[b - 1]; best[a] = it; break; }
        for (int a = 0; a < 5 && best[a] >= 0; a++)
            fprintf(stderr, "top item %ld (job %d branch %d): nodes %lld, %.2f s\n", best[a], items[best[a]].job, items[best[a]].br, itnodes[best[a]], ittime[best[a]]);
    }

    if (wordlen >= 0) {
        qsort(word, wordlen, sizeof(int), cmpint);
        /* recheck: in ker C, outside rowsp D (by the echelon form of D, independent of the class bits) */
        uint8_t *in = xcalloc(n, 1);
        for (int i = 0; i < wordlen; i++) in[word[i]] = 1;
        int ker = 1;
        for (int c = 0; c < m && ker; c++) { int s = 0; for (int i = Coff[c]; i < Coff[c + 1]; i++) s ^= in[Cv[i]]; if (s) ker = 0; }
        memset(v, 0, sizeof(uint64_t) * nw);
        for (int i = 0; i < wordlen; i++) v[word[i] >> 6] |= 1ull << (word[i] & 63);
        ech_reduce(&ED, v);
        int nontriv = !is_zero(v, nw);
        printf("RESULT side=%c W=%d Weff=%d verdict=found weight=%d nodes=%lld leaves=%lld wall=%.2f cpu=%.2f search=%.2f jobs=%d items=%ld threads=%d mode=%s found=%lld recheck_ker=%s recheck_nontrivial=%s\n",
               side, W, Weff, wordlen, nodes, leaves, tw, tc, tsearch, jb - ja, nitems, threads, mode ? "pattern" : "anchor",
               nfound_total, ker ? "ok" : "FAIL", nontriv ? "ok" : "FAIL");
        printf("WORD %d :", wordlen);
        for (int i = 0; i < wordlen; i++) printf(" %d", word[i]);
        printf("\n");
        if (!ker || !nontriv || wordlen > Weff) { printf("ERROR: reported word fails the recheck\n"); return 3; }
        printf("CONCLUSION d_%c <= %d\n", side, wordlen);
    } else {
        printf("RESULT side=%c W=%d Weff=%d verdict=%s nodes=%lld leaves=%lld wall=%.2f cpu=%.2f search=%.2f jobs=%d items=%ld threads=%d mode=%s\n",
               side, W, Weff, complete ? "none" : "none_in_subset", nodes, leaves, tw, tc, tsearch, jb - ja, nitems, threads, mode ? "pattern" : "anchor");
        if (K == 0) printf("CONCLUSION k=0, no logical operators\n");
        else if (complete) {
            int lb = W + 1;
            if (parity && !noparity && (lb & 1)) lb++;
            printf("CONCLUSION none of weight <= %d; d_%c >= %d\n", W, side, lb);
        } else
            printf("CONCLUSION partial run (jobs %d:%d of %d), no bound\n", ja, jb, njobs);
    }
    return 0;
}
