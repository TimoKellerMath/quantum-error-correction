/*
 * Exact search for light Z-logical operators of a CSS code (anchored cluster search).
 *
 * Base algorithm: Rowshan-Devitt, arXiv:2609.39395v2, Algorithm 1 (grow a qubit set T from an
 * anchor, always resolving one unsatisfied X-check, prune by |T| + ceil(u/c) > W, stop a branch at
 * zero syndrome).  Additions (each valid for MINIMUM-weight logicals z containing the anchor):
 *   -x 1  exclusion branching: branch i puts q_i in z and q_1..q_{i-1} out of z, so every z is
 *         reached exactly once (the branches partition the candidates);
 *   -p 1  XOR unit propagation: a check with one free qubit forces it (in iff the check is
 *         unsatisfied), a check with no free qubit must be satisfied;
 *   -m 1  stabilizer majority: for every stabilizer s and minimum-weight z, |z+s| >= |z|, i.e.
 *         |z cap s| <= |s|/2; s ranges over the Z-checks (-m 2: also sums of two overlapping
 *         Z-checks); when |T cap s| = floor(|s|/2) the other qubits of s are forced out;
 *   -b 1,2 stronger lower bounds for the weight still needed (see lower_bound()).
 * A zero-syndrome T is never extended: if T is a boundary it would be a proper subset of a minimum
 * weight z with zero syndrome, and z+T would be lighter in the same class (Rowshan-Devitt, Lemma 2).
 *
 * Instance format (written by codes.py): header "n mX mZ k N a b c", then the supports of the X-checks,
 * of the Z-checks and of a basis of ker H_Z modulo rowsp H_X, one per line as "w q_1 .. q_w".
 * Implementation B of exact_distances.tex. Output: node counts and every logical found.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, mx, mz, K, NN;
static int *xoff, *xq, *qxoff, *qx;      /* X checks and qubit -> X checks */
static int *zoff, *zq;                   /* Z checks */
static int ns, *soff, *sq, *sthr, *qsoff, *qs; /* majority stabilizers */
static uint64_t *cls;                    /* bit j set iff qubit lies in X-logical j */
static int cmax;                         /* max number of X checks per qubit */

/* options */
static int W = 0, anchor = -1, rightonly = 0, rule = 1, excl = 1, prop = 1, majo = 1, bound = 1;
static int enumerate = 0, splitD = -1, P = 1, I = 0, verbose = 0, histo = 0, samp_t = -1, samp_every = 1000;
static long long samp_cnt = 0;
static FILE *sampf = NULL;
static long long hist2[64][64];
static const char *outfn = NULL, *ckfn = NULL;
static long long node_limit = -1;

/* state */
static signed char *st;
static unsigned char *syn, *avail;
static int ucount, *ulist, *upos;
static int *scnt;
static int Tsize;
static uint64_t tcls;
static int *trail, tlen;
static int *cstack, ctop, *sstack, stop_;
static int *mark;
static int markgen = 1;
static unsigned char *hq;          /* number of unsatisfied checks containing the qubit */
static int cntH[8];                /* number of FREE qubits with a given hq */
static uint64_t *ubits;            /* unsatisfied checks as a bitset */
static int La, Lb, Lc;             /* lattice HNF */
static int pairmode = 0, *gadd = NULL, *Flist = NULL, nF = 0, *qdist0 = NULL;

/* statistics */
static long long nodes = 0, generated = 0, nfound = 0, jobid = -1, jobs_done = 0, jobs_skipped = 0;
static long long found_hist[1024];
static int stopflag = 0;
static FILE *outf = NULL, *ckf = NULL;
static char *done = NULL;
static long long done_cap = 0;
static double t_start;

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void die(const char *m) { fprintf(stderr, "error: %s\n", m); exit(1); }

static int *read_rows(FILE *f, int rows, int **off_out) {
    int cap = 64, len = 0;
    int *v = malloc(cap * sizeof(int));
    int *off = malloc((rows + 1) * sizeof(int));
    for (int r = 0; r < rows; r++) {
        int w;
        if (fscanf(f, "%d", &w) != 1) die("bad instance");
        off[r] = len;
        for (int i = 0; i < w; i++) {
            if (len == cap) { cap *= 2; v = realloc(v, cap * sizeof(int)); }
            if (fscanf(f, "%d", &v[len]) != 1) die("bad instance");
            len++;
        }
    }
    off[rows] = len;
    *off_out = off;
    return v;
}

static void transpose(int rows, const int *off, const int *val, int **toff, int **tval) {
    int *cnt = calloc(n + 1, sizeof(int));
    for (int i = 0; i < off[rows]; i++) cnt[val[i] + 1]++;
    for (int q = 0; q < n; q++) cnt[q + 1] += cnt[q];
    int *tv = malloc((off[rows] + 1) * sizeof(int));
    int *pos = malloc((n + 1) * sizeof(int));
    memcpy(pos, cnt, (n + 1) * sizeof(int));
    for (int r = 0; r < rows; r++)
        for (int i = off[r]; i < off[r + 1]; i++) tv[pos[val[i]]++] = r;
    free(pos);
    *toff = cnt;
    *tval = tv;
}

static int cmp_int(const void *a, const void *b) { return (*(const int *)a > *(const int *)b) - (*(const int *)a < *(const int *)b); }

/* majority stabilizers: sums of connected sets of at most `majo` Z checks (connected: the overlap
 * graph of their supports is connected), enumerated once each (ESU algorithm), weight <= wmax */
static int wmax_stab = 0;
static int *ovoff, *ovv;            /* overlap graph of Z checks */
static int st_cap, st_len, st_scap;
static int *esu_cur, esu_ncur, esu_root;
static char *in1;
static int *buf;
static void emit_stab(void) {
    int w = 0;
    for (int k = 0; k < esu_ncur; k++) { int s = esu_cur[k]; for (int i = zoff[s]; i < zoff[s + 1]; i++) in1[zq[i]] ^= 1; }
    for (int k = 0; k < esu_ncur; k++) { int s = esu_cur[k]; for (int i = zoff[s]; i < zoff[s + 1]; i++) if (in1[zq[i]]) { buf[w++] = zq[i]; in1[zq[i]] = 0; } }
    for (int k = 0; k < esu_ncur; k++) { int s = esu_cur[k]; for (int i = zoff[s]; i < zoff[s + 1]; i++) in1[zq[i]] = 0; }
    if (w == 0 || w > wmax_stab) return;
    if (ns == st_scap) { st_scap *= 2; soff = realloc(soff, (st_scap + 1) * sizeof(int)); sthr = realloc(sthr, st_scap * sizeof(int)); }
    while (st_len + w > st_cap) { st_cap *= 2; sq = realloc(sq, st_cap * sizeof(int)); }
    qsort(buf, w, sizeof(int), cmp_int);
    memcpy(sq + st_len, buf, w * sizeof(int));
    st_len += w;
    sthr[ns] = w / 2;
    ns++;
    soff[ns] = st_len;
}
static int adjacent_to_cur(int u) {
    for (int k = 0; k < esu_ncur; k++) {
        int s = esu_cur[k];
        if (s == u) return 1;
        for (int i = ovoff[s]; i < ovoff[s + 1]; i++) if (ovv[i] == u) return 1;
    }
    return 0;
}
static void esu(int *ext, int next) {
    emit_stab();
    if (esu_ncur == majo) return;
    int *ext2 = malloc((next + 64 * (size_t)(majo + 1)) * sizeof(int));
    while (next > 0) {
        int w = ext[--next];
        int n2 = 0;
        for (int i = 0; i < next; i++) ext2[n2++] = ext[i];
        for (int i = ovoff[w]; i < ovoff[w + 1]; i++) {
            int u = ovv[i];
            if (u <= esu_root || adjacent_to_cur(u)) continue;
            int dup = 0;
            for (int j = 0; j < n2; j++) if (ext2[j] == u) { dup = 1; break; }
            if (!dup) ext2[n2++] = u;
        }
        esu_cur[esu_ncur++] = w;
        esu(ext2, n2);
        esu_ncur--;
    }
    free(ext2);
}
static void build_stabs(void) {
    if (wmax_stab <= 0) wmax_stab = 6 * majo;
    st_cap = 1024; st_len = 0; st_scap = 256;
    sq = malloc(st_cap * sizeof(int));
    soff = malloc((st_scap + 1) * sizeof(int));
    sthr = malloc(st_scap * sizeof(int));
    ns = 0;
    soff[0] = 0;
    buf = malloc(2 * n * sizeof(int));
    in1 = calloc(n, 1);
    int *zqo, *zqv;
    transpose(mz, zoff, zq, &zqo, &zqv);
    /* overlap graph */
    ovoff = malloc((mz + 1) * sizeof(int));
    int cap = 16 * mz + 16, len = 0;
    ovv = malloc(cap * sizeof(int));
    int *tmp = malloc(4096 * sizeof(int));
    for (int s = 0; s < mz; s++) {
        ovoff[s] = len;
        int nt = 0;
        for (int i = zoff[s]; i < zoff[s + 1]; i++)
            for (int j = zqo[zq[i]]; j < zqo[zq[i] + 1]; j++) if (zqv[j] != s) tmp[nt++] = zqv[j];
        qsort(tmp, nt, sizeof(int), cmp_int);
        for (int i = 0; i < nt; i++) {
            if (i > 0 && tmp[i] == tmp[i - 1]) continue;
            if (len == cap) { cap *= 2; ovv = realloc(ovv, cap * sizeof(int)); }
            ovv[len++] = tmp[i];
        }
    }
    ovoff[mz] = len;
    free(tmp);
    esu_cur = malloc((majo + 1) * sizeof(int));
    for (int r = 0; r < mz; r++) {
        esu_root = r;
        esu_ncur = 1;
        esu_cur[0] = r;
        int next = 0;
        int *ext = malloc((ovoff[r + 1] - ovoff[r] + 1) * sizeof(int));
        for (int i = ovoff[r]; i < ovoff[r + 1]; i++) if (ovv[i] > r) ext[next++] = ovv[i];
        esu(ext, next);
        free(ext);
    }
    free(in1); free(buf); free(zqo); free(zqv);
    transpose(ns, soff, sq, &qsoff, &qs);
}

static void read_instance(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) die("cannot open instance");
    int a, b, c;
    if (fscanf(f, "%d %d %d %d %d %d %d %d", &n, &mx, &mz, &K, &NN, &a, &b, &c) != 8) die("bad header");
    La = a; Lb = b; Lc = c;
    if (K > 64) die("k > 64 not supported");
    xq = read_rows(f, mx, &xoff);
    zq = read_rows(f, mz, &zoff);
    int *lo, *lq = read_rows(f, K, &lo);
    fclose(f);
    transpose(mx, xoff, xq, &qxoff, &qx);
    cls = calloc(n, sizeof(uint64_t));
    for (int j = 0; j < K; j++)
        for (int i = lo[j]; i < lo[j + 1]; i++) cls[lq[i]] |= 1ULL << j;
    free(lo); free(lq);
    cmax = 0;
    for (int q = 0; q < n; q++) if (qxoff[q + 1] - qxoff[q] > cmax) cmax = qxoff[q + 1] - qxoff[q];
}

static void alloc_state(void) {
    st = calloc(n, 1);
    syn = calloc(mx, 1);
    avail = malloc(mx);
    for (int c = 0; c < mx; c++) avail[c] = (unsigned char)(xoff[c + 1] - xoff[c]);
    ulist = malloc(mx * sizeof(int));
    upos = malloc(mx * sizeof(int));
    scnt = calloc(ns > 0 ? ns : 1, sizeof(int));
    trail = malloc((n + 1) * sizeof(int));
    cstack = malloc((8 * n + 64) * sizeof(int));
    sstack = malloc((8 * n + 64 + 64 * (ns + 1)) * sizeof(int));
    mark = calloc(n, sizeof(int));
    hq = calloc(n, 1);
    memset(cntH, 0, sizeof(cntH));
    cntH[0] = n;
    ubits = calloc((mx + 63) / 64 + 1, sizeof(uint64_t));
    ucount = 0; Tsize = 0; tcls = 0; tlen = 0; ctop = 0; stop_ = 0;
}

static inline void flip(int c) {
    syn[c] ^= 1;
    int d;
    if (syn[c]) { upos[c] = ucount; ulist[ucount++] = c; ubits[c >> 6] |= 1ULL << (c & 63); d = 1; }
    else { int p = upos[c], l = ulist[--ucount]; ulist[p] = l; upos[l] = p; ubits[c >> 6] &= ~(1ULL << (c & 63)); d = -1; }
    for (int i = xoff[c]; i < xoff[c + 1]; i++) {
        int q = xq[i];
        if (st[q] == 0) { cntH[hq[q]]--; hq[q] += d; cntH[hq[q]]++; }
        else hq[q] += d;
    }
}

/* returns 0 on immediate conflict; the assignment is on the trail in any case */
static void set_out(int q);
static int set_in(int q) {
    int ok = 1;
    cntH[hq[q]]--;
    st[q] = 1;
    trail[tlen++] = q;
    Tsize++;
    tcls ^= cls[q];
    for (int i = qxoff[q]; i < qxoff[q + 1]; i++) {
        int c = qx[i];
        flip(c);
        avail[c]--;
        cstack[ctop++] = c;
    }
    for (int i = qsoff ? qsoff[q] : 0; qsoff && i < qsoff[q + 1]; i++) {
        int s = qs[i];
        scnt[s]++;
        if (scnt[s] > sthr[s]) ok = 0;
        else if (scnt[s] == sthr[s]) sstack[stop_++] = s;
    }
    if (Tsize > W) ok = 0;
    if (pairmode && q < NN) {
        /* difference avoidance: no two qubits of block 1 with difference in Flist */
        const int *row = gadd + (size_t)q * NN;
        for (int k = 0; k < nF; k++) {
            int r = row[Flist[k]];
            if (st[r] == 1) ok = 0;
            else if (st[r] == 0) set_out(r);
        }
    }
    return ok;
}

static void set_out(int q) {
    cntH[hq[q]]--;
    st[q] = 2;
    trail[tlen++] = q;
    for (int i = qxoff[q]; i < qxoff[q + 1]; i++) {
        int c = qx[i];
        avail[c]--;
        cstack[ctop++] = c;
    }
}

static void undo(int m) {
    while (tlen > m) {
        int q = trail[--tlen];
        if (st[q] == 1) {
            Tsize--;
            tcls ^= cls[q];
            for (int i = qxoff[q + 1] - 1; i >= qxoff[q]; i--) { int c = qx[i]; flip(c); avail[c]++; }
            for (int i = qsoff ? qsoff[q] : 0; qsoff && i < qsoff[q + 1]; i++) scnt[qs[i]]--;
        } else {
            for (int i = qxoff[q]; i < qxoff[q + 1]; i++) avail[qx[i]]++;
        }
        st[q] = 0;
        cntH[hq[q]]++;
    }
}

/* Lower bound for the number of qubits still to be added.  Every unsatisfied check needs a free
 * qubit of its support in z; a free qubit meets at most h_max unsatisfied checks (bound 1), and
 * unsatisfied checks with pairwise disjoint sets of free qubits need distinct qubits (bound 2). */
/* exact small hitting-set test used by bound 4 */
static int hs_u, hs_chk[64], hs_cov[64];
static long long hs_nodes;
static int hs_rec(int R) {
    hs_nodes++;
    int best = -1, bestn = 99, unc = 0;
    for (int k = 0; k < hs_u; k++) {
        if (hs_cov[k]) continue;
        unc++;
        int c = hs_chk[k], nf = 0;
        for (int i = xoff[c]; i < xoff[c + 1]; i++) if (!st[xq[i]]) nf++;
        if (nf < bestn) { bestn = nf; best = k; }
    }
    if (unc == 0) return 1;
    if (R <= 0 || bestn == 0 || unc > cmax * R) return 0;
    int c = hs_chk[best];
    for (int i = xoff[c]; i < xoff[c + 1]; i++) {
        int q = xq[i];
        if (st[q]) continue;
        int saved[8], ns_ = 0;
        for (int j = qxoff[q]; j < qxoff[q + 1]; j++) {
            int cc = qx[j];
            if (!syn[cc]) continue;
            for (int k = 0; k < hs_u; k++) if (hs_chk[k] == cc && !hs_cov[k]) { hs_cov[k] = 1; saved[ns_++] = k; }
        }
        st[q] = 3;  /* temporarily mark as used so later branches do not reuse it */
        int r = hs_rec(R - 1);
        st[q] = 0;
        for (int m = 0; m < ns_; m++) hs_cov[saved[m]] = 0;
        if (r) return 1;
    }
    return 0;
}
static int hitting_set_le(int R) {
    hs_u = ucount < 64 ? ucount : 64;
    for (int k = 0; k < hs_u; k++) { hs_chk[k] = ulist[k]; hs_cov[k] = 0; }
    return hs_rec(R);
}

static int lower_bound(void) {
    int u = ucount;
    if (u == 0) return 0;
    if (bound == 0) return (u + cmax - 1) / cmax;
    int hmax = 0;
    for (int h = 7; h >= 1; h--) if (cntH[h] > 0) { hmax = h; break; }
    if (hmax == 0) return 1 << 20;
    int lb = (u + hmax - 1) / hmax;
    if (bound >= 2 && u > lb) {
        markgen++;
        int pk = 0;
        for (int wd = 0; wd <= (mx - 1) >> 6; wd++) {
            uint64_t bits = ubits[wd];
            while (bits) {
                int c = (wd << 6) + __builtin_ctzll(bits);
                bits &= bits - 1;
                int free_ok = 1;
                for (int i = xoff[c]; i < xoff[c + 1]; i++) { int q = xq[i]; if (!st[q] && mark[q] == markgen) { free_ok = 0; break; } }
                if (!free_ok) continue;
                pk++;
                for (int i = xoff[c]; i < xoff[c + 1]; i++) { int q = xq[i]; if (!st[q]) mark[q] = markgen; }
            }
        }
        if (pk > lb) lb = pk;
    }
    if (bound >= 3 && lb <= W - Tsize) {
        /* degree-sequence bound: a free qubit covers at most h of the unsatisfied checks, where h is
         * the number of its checks that are unsatisfied; t qubits cover at most the t largest h. */
        int cnt[4] = {0, 0, 0, 0};
        markgen++;
        for (int k = 0; k < u; k++) {
            int c = ulist[k];
            for (int i = xoff[c]; i < xoff[c + 1]; i++) {
                int q = xq[i];
                if (st[q] || mark[q] == markgen) continue;
                mark[q] = markgen;
                int h = 0;
                for (int j = qxoff[q]; j < qxoff[q + 1]; j++) h += syn[qx[j]];
                cnt[h > 3 ? 3 : h]++;
            }
        }
        int t = 0, cov = 0;
        for (int h = 3; h >= 1 && cov < u; h--)
            for (int m = 0; m < cnt[h] && cov < u; m++) { t++; cov += h; }
        if (cov < u) return 1 << 20;
        if (t > lb) lb = t;
    }
    if (bound >= 4 && lb <= W - Tsize && u <= 40) {
        /* exact decision: is there a set of at most R free qubits meeting every unsatisfied check? */
        int R = W - Tsize;
        if (!hitting_set_le(R)) return R + 1;
    }
    return lb;
}

static int propagate(int ok) {
    if (!ok) { ctop = 0; stop_ = 0; return 0; }
    while (ctop > 0 || stop_ > 0) {
        if (stop_ > 0) {
            int s = sstack[--stop_];
            if (scnt[s] > sthr[s]) { ctop = 0; stop_ = 0; return 0; }
            if (scnt[s] == sthr[s] && prop)
                for (int i = soff[s]; i < soff[s + 1]; i++) if (st[sq[i]] == 0) set_out(sq[i]);
            continue;
        }
        int c = cstack[--ctop];
        if (avail[c] == 0) {
            if (syn[c]) { ctop = 0; stop_ = 0; return 0; }
            continue;
        }
        if (!prop || avail[c] != 1) continue;
        int q = -1;
        for (int i = xoff[c]; i < xoff[c + 1]; i++) if (st[xq[i]] == 0) { q = xq[i]; break; }
        if (syn[c]) { if (!set_in(q)) { ctop = 0; stop_ = 0; return 0; } }
        else set_out(q);
    }
    if (Tsize + lower_bound() > W) return 0;
    return 1;
}

static void report(void) {
    nfound++;
    if (Tsize < 1024) found_hist[Tsize]++;
    if (outf) {
        int *s = malloc(Tsize * sizeof(int)), m = 0;
        for (int i = 0; i < tlen; i++) if (st[trail[i]] == 1) s[m++] = trail[i];
        qsort(s, m, sizeof(int), cmp_int);
        fprintf(outf, "%d %llx :", Tsize, (unsigned long long)tcls);
        for (int i = 0; i < m; i++) fprintf(outf, " %d", s[i]);
        fprintf(outf, "\n");
        fflush(outf);
        free(s);
    }
    if (!enumerate) stopflag = 1;
}

static int choose_check(void) {
    int best = -1;
    for (int k = 0; k < ucount; k++) {
        int c = ulist[k];
        if (best < 0) { best = c; continue; }
        if (rule == 0) { if (c < best) best = c; }
        else if (avail[c] < avail[best] || (avail[c] == avail[best] && c < best)) best = c;
    }
    return best;
}

static void search(int depth);

static void subtree(int depth) {
    if (depth == splitD) {
        jobid++;
        if (jobid % P != I) return;
        if (jobid < done_cap && done[jobid]) { jobs_skipped++; return; }
        long long n0 = nodes, f0 = nfound;
        search(depth);
        jobs_done++;
        if (ckf && !stopflag) {
            fprintf(ckf, "%lld %lld %lld\n", jobid, nodes - n0, nfound - f0);
            fflush(ckf);
        }
        if (verbose && jobs_done % 1000 == 0)
            fprintf(stderr, "  jobs done %lld (last id %lld), nodes %lld, found %lld, %.1f s\n", jobs_done, jobid, nodes, nfound, now() - t_start);
        return;
    }
    search(depth);
}

static void search(int depth) {
    if (stopflag) return;
    nodes++;
    if (node_limit > 0 && nodes >= node_limit) { stopflag = 2; return; }
    if (verbose && (nodes & ((1LL << 28) - 1)) == 0)
        fprintf(stderr, "  nodes %lld, found %lld, jobs %lld, %.1f s\n", nodes, nfound, jobs_done, now() - t_start);
    if (histo) hist2[Tsize < 63 ? Tsize : 63][ucount < 63 ? ucount : 63]++;
    if (Tsize == samp_t && ucount > 0 && (samp_cnt++ % samp_every) == 0) {
        if (!sampf) sampf = fopen("samples.txt", "w");
        fprintf(sampf, "T");
        for (int q = 0; q < n; q++) if (st[q] == 1) fprintf(sampf, " %d", q);
        fprintf(sampf, " X");
        for (int q = 0; q < n; q++) if (st[q] == 2) fprintf(sampf, " %d", q);
        fprintf(sampf, " u %d\n", ucount);
    }
    if (ucount == 0) {
        if (tcls) report();
        return;
    }
    int c = choose_check();
    int cand[64], r = 0;
    for (int i = xoff[c]; i < xoff[c + 1]; i++) if (st[xq[i]] == 0) cand[r++] = xq[i];
    int m0 = tlen;
    for (int i = 0; i < r && !stopflag; i++) {
        int q = cand[i];
        if (st[q] == 2) continue;
        int already = (st[q] == 1);
        int m2 = tlen;
        generated++;
        if (already || propagate(set_in(q))) subtree(depth + 1);
        undo(m2);
        if (already) break;
        if (!excl) continue;
        if (i < r - 1) {
            set_out(q);
            if (!propagate(1)) break;
        }
    }
    undo(m0);
}

/* group G = Z^2/Lambda, element index i*Lc + j for the representative (i, j) */
static int gindex(long u, long v) {
    long q = v >= 0 ? v / Lc : -((-v + Lc - 1) / Lc);
    long vv = v - q * Lc;
    long uu = ((u - q * Lb) % La + La) % La;
    return (int)(uu * Lc + vv);
}

/* Pair-anchored search.  For a minimum-weight z with at least two qubits in block 1 let D(z) be the
 * set of differences of its block-1 positions; with the classes u_1, u_2, ... of {g, -g} (g != 0)
 * in a fixed order and j = min{ j : u_j meets D(z) }, a translate of z contains the block-1 qubits
 * 0 and g_j and no two block-1 qubits with difference in u_1 + ... + u_{j-1}.  Job j searches
 * exactly these; job 0 covers z with one block-1 qubit; the run -R covers block 2 only. */
static void pair_driver(void) {
    gadd = malloc((size_t)NN * NN * sizeof(int));
    for (int t1 = 0; t1 < NN; t1++)
        for (int t2 = 0; t2 < NN; t2++)
            gadd[(size_t)t1 * NN + t2] = gindex(t1 / Lc + t2 / Lc, t1 % Lc + t2 % Lc);
    int *neg = malloc(NN * sizeof(int));
    for (int t = 0; t < NN; t++) neg[t] = gindex(-(t / Lc), -(t % Lc));
    for (int t = 0; t < NN; t++) if (gadd[(size_t)t * NN + neg[t]] != 0) die("group table");
    /* BFS distances from qubit 0 in the qubit graph (adjacent: common X check) */
    qdist0 = malloc(n * sizeof(int));
    for (int q = 0; q < n; q++) qdist0[q] = -1;
    int *queue = malloc(n * sizeof(int)), qh = 0, qt = 0;
    qdist0[0] = 0; queue[qt++] = 0;
    while (qh < qt) {
        int q = queue[qh++];
        for (int i = qxoff[q]; i < qxoff[q + 1]; i++) {
            int c = qx[i];
            for (int k = xoff[c]; k < xoff[c + 1]; k++) { int r = xq[k]; if (qdist0[r] < 0) { qdist0[r] = qdist0[q] + 1; queue[qt++] = r; } }
        }
    }
    free(queue);
    /* classes {g, -g}, ordered by (distance of block-1 qubit g from qubit 0, g) */
    int *cls_rep = malloc(NN * sizeof(int)), M = 0;
    char *seen = calloc(NN, 1);
    for (int t = 1; t < NN; t++) if (!seen[t]) { seen[t] = 1; seen[neg[t]] = 1; cls_rep[M++] = t; }
    for (int i = 1; i < M; i++) {      /* insertion sort */
        int x = cls_rep[i], j = i - 1;
        while (j >= 0 && (qdist0[cls_rep[j]] > qdist0[x] || (qdist0[cls_rep[j]] == qdist0[x] && cls_rep[j] > x))) { cls_rep[j + 1] = cls_rep[j]; j--; }
        cls_rep[j + 1] = x;
    }
    Flist = malloc(2 * NN * sizeof(int));
    printf("pair mode: %d difference classes\n", M);
    for (int job = 0; job <= M && !stopflag; job++) {
        if (job % P != I) continue;
        if (job < done_cap && done[job]) { jobs_skipped++; continue; }
        long long n0 = nodes, f0 = nfound;
        double tj = now();
        int ok = 1;
        nF = 0;
        if (job == 0) {
            pairmode = 0;
            for (int q = 1; q < NN; q++) if (st[q] == 0) set_out(q);
            ok = propagate(1);
            if (ok) ok = propagate(set_in(0));
        } else {
            int g = cls_rep[job - 1];
            for (int i = 0; i < job - 1; i++) {
                Flist[nF++] = cls_rep[i];
                if (neg[cls_rep[i]] != cls_rep[i]) Flist[nF++] = neg[cls_rep[i]];
            }
            pairmode = 1;
            if (qdist0[g] >= W) ok = 0;   /* a connected set containing qubits 0 and g has > W elements */
            if (ok) ok = propagate(set_in(0));
            if (ok && st[g] == 0) ok = propagate(set_in(g));
            else if (ok && st[g] != 1) ok = 0;
        }
        generated++;
        if (ok) search(0);
        undo(0);
        jobs_done++;
        if (ckf && !stopflag) { fprintf(ckf, "%d %lld %lld\n", job, nodes - n0, nfound - f0); fflush(ckf); }
        if (verbose) fprintf(stderr, "  job %d/%d g=%d dist=%d: nodes %lld, found %lld, %.2f s\n", job, M,
                             job ? cls_rep[job - 1] : 0, job ? qdist0[cls_rep[job - 1]] : 0, nodes - n0, nfound - f0, now() - tj);
    }
    pairmode = 0;
    free(seen); free(cls_rep); free(neg);
}

static void usage(void) {
    fprintf(stderr,
            "usage: cluster INSTANCE -W w [-a anchor] [-R] [-r rule] [-x 0|1] [-p 0|1] [-m 0|1|2] [-b 0|1|2]\n"
            "       [-e] [-o out] [-D depth -P workers -I id -C checkpoint] [-L nodelimit] [-v] [-H] [-M wmax]\n");
    exit(1);
}

int main(int argc, char **argv) {
    if (argc < 2) usage();
    const char *inst = argv[1];
    for (int i = 2; i < argc; i++) {
        char *o = argv[i];
        #define ARG (i + 1 < argc ? argv[++i] : (usage(), (char *)0))
        if (!strcmp(o, "-W")) W = atoi(ARG);
        else if (!strcmp(o, "-a")) anchor = atoi(ARG);
        else if (!strcmp(o, "-R")) rightonly = 1;
        else if (!strcmp(o, "-r")) rule = atoi(ARG);
        else if (!strcmp(o, "-x")) excl = atoi(ARG);
        else if (!strcmp(o, "-p")) prop = atoi(ARG);
        else if (!strcmp(o, "-m")) majo = atoi(ARG);
        else if (!strcmp(o, "-b")) bound = atoi(ARG);
        else if (!strcmp(o, "-e")) enumerate = 1;
        else if (!strcmp(o, "-o")) outfn = ARG;
        else if (!strcmp(o, "-D")) splitD = atoi(ARG);
        else if (!strcmp(o, "-P")) P = atoi(ARG);
        else if (!strcmp(o, "-I")) I = atoi(ARG);
        else if (!strcmp(o, "-C")) ckfn = ARG;
        else if (!strcmp(o, "-L")) node_limit = atoll(ARG);
        else if (!strcmp(o, "-v")) verbose = 1;
        else if (!strcmp(o, "-H")) histo = 1;
        else if (!strcmp(o, "-S")) { samp_t = atoi(ARG); samp_every = atoi(ARG); }
        else if (!strcmp(o, "-M")) wmax_stab = atoi(ARG);
        else if (!strcmp(o, "-A")) pairmode = 1;
        else usage();
    }
    read_instance(inst);
    if (majo) build_stabs(); else { ns = 0; qsoff = NULL; }
    alloc_state();
    if (anchor < 0) anchor = rightonly ? NN : 0;
    if (ckfn) {
        FILE *f = fopen(ckfn, "r");
        if (f) {
            long long id, a, b;
            while (fscanf(f, "%lld %lld %lld", &id, &a, &b) == 3) {
                if (id >= done_cap) {
                    long long nc = done_cap ? done_cap : 1024;
                    while (nc <= id) nc *= 2;
                    done = realloc(done, nc);
                    memset(done + done_cap, 0, nc - done_cap);
                    done_cap = nc;
                }
                done[id] = 1;
            }
            fclose(f);
        }
        ckf = fopen(ckfn, "a");
    }
    if (outfn) outf = fopen(outfn, enumerate ? "w" : "a");
    t_start = now();
    printf("instance %s: n=%d mx=%d mz=%d k=%d N=%d; W=%d anchor=%d rightonly=%d rule=%d excl=%d prop=%d majority=%d (%d stabs) bound=%d enumerate=%d split=%d P=%d I=%d\n",
           inst, n, mx, mz, K, NN, W, anchor, rightonly, rule, excl, prop, majo, ns, bound, enumerate, splitD, P, I);
    fflush(stdout);
    int ok = 1;
    if (pairmode) {
        if (rightonly) die("-A and -R are separate runs");
        pair_driver();
        ok = 0;
    }
    if (ok && rightonly) {
        for (int q = 0; q < NN; q++) if (st[q] == 0) set_out(q);
        ok = propagate(1);
    }
    if (pairmode) ;
    else if (ok && st[anchor] == 0) ok = propagate(set_in(anchor));
    else if (st[anchor] != 1) ok = 0;
    generated++;
    if (ok) subtree(0);
    double dt = now() - t_start;
    int minw = -1;
    for (int w = 0; w < 1024; w++) if (found_hist[w]) { minw = w; break; }
    printf("RESULT nodes=%lld generated=%lld found=%lld minw=%d time=%.2f stop=%d jobs_done=%lld jobs_total_seen=%lld jobs_skipped=%lld\n",
           nodes, generated, nfound, minw, dt, stopflag, jobs_done, jobid + 1, jobs_skipped);
    for (int w = 0; w < 1024; w++) if (found_hist[w]) printf("  weight %d: %lld\n", w, found_hist[w]);
    if (stopflag == 2) printf("NODE LIMIT REACHED: search incomplete\n");
    if (histo) {
        printf("HISTOGRAM rows |T|, columns u = number of unsatisfied checks (0..20)\n");
        for (int t = 0; t < 64; t++) {
            long long rs = 0;
            for (int u = 0; u < 64; u++) rs += hist2[t][u];
            if (!rs) continue;
            printf("  |T|=%2d total %12lld :", t, rs);
            for (int u = 0; u <= 20; u++) printf(" %lld", hist2[t][u]);
            printf("\n");
        }
    }
    if (outf) fclose(outf);
    if (ckf) fclose(ckf);
    if (sampf) fclose(sampf);
    return 0;
}
