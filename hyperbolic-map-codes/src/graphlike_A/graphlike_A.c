/* graphlike_A: exact d_X and d_Z of a CSS code on every side whose check matrix
 * has all column weights <= 2, by a Horton-set search on Gamma^+ with k-bit
 * logical annotations.  Correctness proof: NOTES.md (Lemmas 1-3, Proposition 5, Lemma 7;
 * counting: Proposition 8).
 *
 *   graphlike_A NAME.hx NAME.hz [--side X|Z|both] [--quiet] [--count | --count-limit NODES]
 *
 * --count also returns the number of minimum-weight logicals of each side (DFS over circuits of
 * length d, at most NODES search nodes, default 1e9; null if the limit is reached).
 *
 * Conventions: d_Z = min weight of ker H_X \ rowsp H_Z (cycles of Gamma^+_X),
 *              d_X = min weight of ker H_Z \ rowsp H_X (cycles of Gamma^+_Z).
 */
#define _GNU_SOURCE
#include <limits.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint64_t u64;

static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "graphlike_A: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(2);
}

static void *xmalloc(size_t s)
{
    void *p = malloc(s ? s : 1);
    if (!p) die("out of memory");
    return p;
}

static void *xcalloc(size_t a, size_t b)
{
    void *p = calloc(a ? a : 1, b ? b : 1);
    if (!p) die("out of memory");
    return p;
}

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}


typedef struct {
    int n, m;
    int *ptr, *idx;
} Sp;

static int cmpint(const void *a, const void *b)
{
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

static int blank(const char *p)
{
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    return *p == 0;
}

static Sp read_sp(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open %s", path);
    char *line = NULL;
    size_t cap = 0;
    if (getline(&line, &cap, f) < 0) die("%s: empty file", path);
    Sp s;
    char extra;
    if (sscanf(line, "%d %d %c", &s.n, &s.m, &extra) != 2 || s.n < 0 || s.m < 0)
        die("%s: bad header line", path);
    s.ptr = xmalloc(((size_t)s.m + 1) * sizeof(int));
    size_t capi = 1024, nnz = 0;
    s.idx = xmalloc(capi * sizeof(int));
    s.ptr[0] = 0;
    for (int i = 0; i < s.m; i++) {
        if (getline(&line, &cap, f) < 0) die("%s: expected %d rows, found %d", path, s.m, i);
        char *p = line, *e;
        for (;;) {
            long v = strtol(p, &e, 10);
            if (e == p) break;
            if (v < 0 || v >= s.n) die("%s: row %d: index %ld out of range", path, i, v);
            if (nnz == capi) {
                capi *= 2;
                s.idx = realloc(s.idx, capi * sizeof(int));
                if (!s.idx) die("out of memory");
            }
            s.idx[nnz++] = (int)v;
            p = e;
        }
        if (!blank(p)) die("%s: row %d: unexpected text '%.20s'", path, i, p);
        s.ptr[i + 1] = (int)nnz;
        qsort(s.idx + s.ptr[i], nnz - (size_t)s.ptr[i], sizeof(int), cmpint);
        for (int t = s.ptr[i] + 1; t < s.ptr[i + 1]; t++)
            if (s.idx[t] == s.idx[t - 1]) die("%s: row %d: repeated index %d", path, i, s.idx[t]);
    }
    while (getline(&line, &cap, f) >= 0)
        if (!blank(line)) die("%s: more than %d rows", path, s.m);
    free(line);
    fclose(f);
    return s;
}


static inline int getbit(const u64 *v, int j) { return (int)((v[j >> 6] >> (j & 63)) & 1); }
static inline void flipbit(u64 *v, int j) { v[j >> 6] ^= 1ULL << (j & 63); }

static int lowbit(const u64 *v, int W)
{
    for (int w = 0; w < W; w++)
        if (v[w]) return w * 64 + __builtin_ctzll(v[w]);
    return -1;
}

static int popc(const u64 *v, int W)
{
    int s = 0;
    for (int w = 0; w < W; w++) s += __builtin_popcountll(v[w]);
    return s;
}

static u64 *dense_rows(const Sp *H, int W)
{
    u64 *a = xcalloc((size_t)H->m * W, sizeof(u64));
    for (int i = 0; i < H->m; i++)
        for (int t = H->ptr[i]; t < H->ptr[i + 1]; t++) flipbit(a + (size_t)i * W, H->idx[t]);
    return a;
}

/* In-place reduced row echelon form; returns the rank, pivot columns in piv. */
static int rref(u64 *a, int rows, int W, int ncols, int *piv)
{
    int r = 0;
    for (int c = 0; c < ncols && r < rows; c++) {
        int w = c >> 6;
        u64 b = 1ULL << (c & 63);
        int p = -1;
        for (int i = r; i < rows; i++)
            if (a[(size_t)i * W + w] & b) { p = i; break; }
        if (p < 0) continue;
        u64 *pr = a + (size_t)r * W;
        if (p != r) {
            u64 *pp = a + (size_t)p * W;
            for (int t = 0; t < W; t++) { u64 x = pr[t]; pr[t] = pp[t]; pp[t] = x; }
        }
        /* the pivot row vanishes in all columns < c, so words < w can be skipped */
        for (int i = 0; i < rows; i++) {
            if (i == r) continue;
            u64 *ri = a + (size_t)i * W;
            if (ri[w] & b)
                for (int t = w; t < W; t++) ri[t] ^= pr[t];
        }
        piv[r++] = c;
    }
    return r;
}

/* Echelon basis kept in insertion order: row i vanishes at the pivots of all rows < i.
 * Reducing in insertion order then clears every pivot bit. */
typedef struct {
    int W, rows, cap;
    u64 *a;
    int *piv;
} Basis;

static void basis_init(Basis *B, int W, int cap)
{
    B->W = W;
    B->rows = 0;
    B->cap = cap > 0 ? cap : 1;
    B->a = xcalloc((size_t)B->cap * W, sizeof(u64));
    B->piv = xmalloc((size_t)B->cap * sizeof(int));
}

static void basis_free(Basis *B)
{
    free(B->a);
    free(B->piv);
}

static int basis_reduce(const Basis *B, u64 *u)
{
    const int W = B->W;
    for (int i = 0; i < B->rows; i++) {
        int p = B->piv[i];
        if ((u[p >> 6] >> (p & 63)) & 1) {
            const u64 *r = B->a + (size_t)i * W;
            for (int w = 0; w < W; w++) u[w] ^= r[w];
        }
    }
    return lowbit(u, W) >= 0;
}

static void basis_push(Basis *B, const u64 *u)
{
    if (B->rows == B->cap) {
        B->cap *= 2;
        B->a = realloc(B->a, (size_t)B->cap * B->W * sizeof(u64));
        B->piv = realloc(B->piv, (size_t)B->cap * sizeof(int));
        if (!B->a || !B->piv) die("out of memory");
    }
    memcpy(B->a + (size_t)B->rows * B->W, u, (size_t)B->W * sizeof(u64));
    B->piv[B->rows++] = lowbit(u, B->W);
}

/* RREF of a check matrix: rows 0..rank-1 of a, pivots piv (an RREF satisfies the Basis invariant). */
typedef struct {
    int rank, W;
    u64 *a;
    int *piv;
} Ech;

static Ech echelon(const Sp *H, int W)
{
    Ech E;
    E.W = W;
    E.a = dense_rows(H, W);
    E.piv = xmalloc(((size_t)H->m + 1) * sizeof(int));
    E.rank = rref(E.a, H->m, W, H->n, E.piv);
    return E;
}

static Basis ech_as_basis(const Ech *E)
{
    Basis B;
    B.W = E->W;
    B.rows = E->rank;
    B.cap = E->rank > 0 ? E->rank : 1;
    B.a = E->a;
    B.piv = E->piv;
    return B;
}

static int sp_parity(const Sp *H, int row, const u64 *v)
{
    int s = 0;
    for (int t = H->ptr[row]; t < H->ptr[row + 1]; t++) s ^= getbit(v, H->idx[t]);
    return s;
}

static int in_kernel(const Sp *H, const u64 *v)
{
    for (int i = 0; i < H->m; i++)
        if (sp_parity(H, i, v)) return 0;
    return 1;
}


enum { ST_OK = 0, ST_NOTGRAPH = 1, ST_K0 = 2, ST_ERR = 3 };

typedef struct {
    int status;
    char reason[256];
    int d;
    int *wit, wlen;
    int verified, lcheck;
    double secs;
    long long nbfs, ncand;
    int cstatus; /* 0 not requested, 1 counted, 2 node limit reached */
    long long count, cnodes;
} SideRes;

static long long count_limit = 0;

static inline void xor3(u64 *restrict out, const u64 *restrict a, const u64 *restrict b, int W)
{
    for (int w = 0; w < W; w++) out[w] = a[w] ^ b[w];
}

static inline int nonzero3(const u64 *restrict a, const u64 *restrict b, const u64 *restrict c, int W)
{
    u64 acc = 0;
    for (int w = 0; w < W; w++) acc |= a[w] ^ b[w] ^ c[w];
    return acc != 0;
}

/* Number of circuits of length d with nonzero annotation, i.e. (Lemma 3 of NOTES.md) the number
 * of logicals of minimum weight d.  Each circuit is enumerated from its least vertex v, inside
 * G[>= v], in both orientations. */
static long long count_min(int V, const int *aptr, const int *adjv, const int *adje, const u64 *Ann, int Wk, int d,
                           int n, const int *cnt, long long limit, long long *nodes, int *aborted)
{
    long long total = 0;
    *nodes = 0, *aborted = 0;
    if (d == 1) {
        for (int j = 0; j < n; j++)
            if (cnt[j] == 0) {
                u64 acc = 0;
                for (int w = 0; w < Wk; w++) acc |= Ann[(size_t)j * Wk + w];
                total += acc != 0;
            }
        return total;
    }
    const int H = d / 2;
    int *stamp = xmalloc((size_t)V * sizeof(int)), *dist = xmalloc((size_t)V * sizeof(int)),
        *queue = xmalloc((size_t)V * sizeof(int)), *onpath = xcalloc((size_t)V, sizeof(int));
    int *pv = xmalloc(((size_t)d + 1) * sizeof(int)), *pe = xmalloc(((size_t)d + 1) * sizeof(int)),
        *pit = xmalloc(((size_t)d + 1) * sizeof(int));
    u64 *PA = xcalloc(((size_t)d + 1) * Wk, sizeof(u64));
    for (int u = 0; u < V; u++) stamp[u] = -1;
    for (int v = 0; v < V && !*aborted; v++) {
        stamp[v] = v, dist[v] = 0, queue[0] = v;
        int qh = 0, qt = 1;
        while (qh < qt) {
            const int x = queue[qh++];
            if (dist[x] >= H) continue;
            for (int t = aptr[x]; t < aptr[x + 1]; t++) {
                const int y = adjv[t];
                if (y < v || stamp[y] == v) continue;
                stamp[y] = v, dist[y] = dist[x] + 1, queue[qt++] = y;
            }
        }
        /* a vertex of a circuit of length d through v lies within distance d/2 of v in G[>= v] */
        int t = 0;
        pv[0] = v, pit[0] = aptr[v], onpath[v] = 1;
        memset(PA, 0, (size_t)Wk * sizeof(u64));
        while (t >= 0) {
            const int u = pv[t];
            if (pit[t] == aptr[u + 1]) {
                onpath[u] = 0;
                t--;
                continue;
            }
            const int idx = pit[t]++, y = adjv[idx], e = adje[idx];
            if (y < v) continue;
            if (y == v) {
                if (t + 1 == d && !(t == 1 && e == pe[1])) {
                    const u64 *a = PA + (size_t)t * Wk, *b = Ann + (size_t)e * Wk;
                    u64 acc = 0;
                    for (int w = 0; w < Wk; w++) acc |= a[w] ^ b[w];
                    total += acc != 0;
                }
                continue;
            }
            if (t + 1 >= d || onpath[y] || stamp[y] != v || t + 1 + dist[y] > d) continue;
            if (++*nodes > limit) {
                *aborted = 1;
                break;
            }
            t++;
            pv[t] = y, pe[t] = e, pit[t] = aptr[y], onpath[y] = 1;
            xor3(PA + (size_t)t * Wk, PA + (size_t)(t - 1) * Wk, Ann + (size_t)e * Wk, Wk);
        }
        for (int i = 0; i <= t; i++) onpath[pv[i]] = 0;
    }
    free(stamp), free(dist), free(queue), free(onpath), free(pv), free(pe), free(pit), free(PA);
    return total / 2;
}

/* Hg: graph-side matrix (its kernel is the cycle space searched), Ho: the other matrix.
 * Finds min weight of ker Hg \ rowsp Ho. */
static void run_side(const Sp *Hg, const Sp *Ho, const Ech *Eg, const Ech *Eo, int n, int k, const char *gname,
                     SideRes *res)
{
    double t0 = now();
    memset(res, 0, sizeof *res);
    res->d = -1;
    const int Wn = (n + 63) / 64;

    int *cnt = xcalloc((size_t)n, sizeof(int)), *c1 = xmalloc((size_t)n * sizeof(int) + 4),
        *c2 = xmalloc((size_t)n * sizeof(int) + 4);
    for (int i = 0; i < Hg->m; i++)
        for (int t = Hg->ptr[i]; t < Hg->ptr[i + 1]; t++) {
            int j = Hg->idx[t];
            if (cnt[j] == 0) c1[j] = i;
            else if (cnt[j] == 1) c2[j] = i;
            cnt[j]++;
        }
    for (int j = 0; j < n; j++)
        if (cnt[j] > 2) {
            res->status = ST_NOTGRAPH;
            snprintf(res->reason, sizeof res->reason, "%s column %d has weight %d > 2", gname, j, cnt[j]);
            free(cnt), free(c1), free(c2);
            res->secs = now() - t0;
            return;
        }
    if (k == 0) {
        res->status = ST_K0;
        snprintf(res->reason, sizeof res->reason, "k = 0");
        free(cnt), free(c1), free(c2);
        res->secs = now() - t0;
        return;
    }

    /* L: k vectors spanning ker Ho modulo rowsp Hg (NOTES.md, Lemma 2) */
    const int ro = Eo->rank;
    int *ispiv = xcalloc((size_t)n, sizeof(int));
    for (int i = 0; i < ro; i++) ispiv[Eo->piv[i]] = 1;
    Basis B;
    basis_init(&B, Wn, Eg->rank + k + 1);
    for (int i = 0; i < Eg->rank; i++) basis_push(&B, Eg->a + (size_t)i * Wn);
    u64 *L = xcalloc((size_t)k * Wn, sizeof(u64));
    int nl = 0;
    u64 *x = xmalloc((size_t)Wn * sizeof(u64));
    for (int f = 0; f < n; f++) {
        if (ispiv[f]) continue;
        /* null vector of the RREF of Ho with free coordinate f */
        memset(x, 0, (size_t)Wn * sizeof(u64));
        flipbit(x, f);
        for (int i = 0; i < ro; i++)
            if (getbit(Eo->a + (size_t)i * Wn, f)) flipbit(x, Eo->piv[i]);
        if (basis_reduce(&B, x)) {
            if (nl == k) die("%s side: more than k = %d independent logicals (inconsistent ranks)", gname, k);
            memcpy(L + (size_t)nl * Wn, x, (size_t)Wn * sizeof(u64));
            nl++;
            basis_push(&B, x);
        }
    }
    basis_free(&B);
    free(ispiv);
    if (nl != k) {
        res->status = ST_ERR;
        snprintf(res->reason, sizeof res->reason, "found %d logicals, expected k = %d", nl, k);
        free(cnt), free(c1), free(c2), free(L), free(x);
        res->secs = now() - t0;
        return;
    }
    /* hypotheses of Lemma 2: rows of L in ker Ho, rank [Hg; L] = rank Hg + k (fresh elimination) */
    {
        int ok = 1;
        for (int t = 0; t < k && ok; t++) ok = in_kernel(Ho, L + (size_t)t * Wn);
        if (ok) {
            u64 *S = xcalloc(((size_t)Hg->m + k) * Wn, sizeof(u64));
            u64 *D = dense_rows(Hg, Wn);
            memcpy(S, D, (size_t)Hg->m * Wn * sizeof(u64));
            memcpy(S + (size_t)Hg->m * Wn, L, (size_t)k * Wn * sizeof(u64));
            int *pv = xmalloc(((size_t)Hg->m + k + 1) * sizeof(int));
            int r = rref(S, Hg->m + k, Wn, n, pv);
            ok = (r == Eg->rank + k);
            free(S), free(D), free(pv);
        }
        res->lcheck = ok;
        if (!ok) {
            res->status = ST_ERR;
            snprintf(res->reason, sizeof res->reason, "logical basis failed its certificate");
            free(cnt), free(c1), free(c2), free(L), free(x);
            res->secs = now() - t0;
            return;
        }
    }

    /* annotation a(e_j) = L e_j, column j of L */
    const int Wk = (k + 63) / 64;
    u64 *Ann = xcalloc((size_t)n * Wk, sizeof(u64));
    for (int t = 0; t < k; t++) {
        const u64 *row = L + (size_t)t * Wn;
        for (int w = 0; w < Wn; w++) {
            u64 b = row[w];
            while (b) {
                int j = w * 64 + __builtin_ctzll(b);
                b &= b - 1;
                Ann[(size_t)j * Wk + (t >> 6)] |= 1ULL << (t & 63);
            }
        }
    }

    /* Gamma^+: vertex 0 = infinity, checks relabelled by decreasing degree (the order of Prop. 5) */
    const int m = Hg->m, V = m + 1;
    int *deg = xcalloc((size_t)m, sizeof(int));
    for (int j = 0; j < n; j++) {
        if (cnt[j] >= 1) deg[c1[j]]++;
        if (cnt[j] == 2) deg[c2[j]]++;
    }
    int *ord = xmalloc((size_t)m * sizeof(int) + 4);
    for (int i = 0; i < m; i++) ord[i] = i;
    {
        int maxd = 0;
        for (int i = 0; i < m; i++)
            if (deg[i] > maxd) maxd = deg[i];
        int *bucket = xcalloc((size_t)maxd + 2, sizeof(int));
        for (int i = 0; i < m; i++) bucket[maxd - deg[i] + 1]++;
        for (int d = 1; d <= maxd + 1; d++) bucket[d] += bucket[d - 1];
        for (int i = 0; i < m; i++) ord[bucket[maxd - deg[i]]++] = i;
        free(bucket);
    }
    int *newid = xmalloc((size_t)m * sizeof(int) + 4);
    for (int p = 0; p < m; p++) newid[ord[p]] = p + 1;
    int *eu = xmalloc((size_t)n * sizeof(int) + 4), *ev = xmalloc((size_t)n * sizeof(int) + 4);
    int *vdeg = xcalloc((size_t)V, sizeof(int));
    for (int j = 0; j < n; j++) {
        if (cnt[j] == 0) { eu[j] = ev[j] = 0; continue; }
        eu[j] = newid[c1[j]];
        ev[j] = cnt[j] == 2 ? newid[c2[j]] : 0;
        vdeg[eu[j]]++, vdeg[ev[j]]++;
    }
    int *aptr = xcalloc((size_t)V + 1, sizeof(int));
    for (int u = 0; u < V; u++) aptr[u + 1] = aptr[u] + vdeg[u];
    int *fill = xmalloc((size_t)V * sizeof(int));
    memcpy(fill, aptr, (size_t)V * sizeof(int));
    int *adjv = xmalloc((size_t)aptr[V] * sizeof(int) + 4), *adje = xmalloc((size_t)aptr[V] * sizeof(int) + 4);
    for (int j = 0; j < n; j++) {
        if (cnt[j] == 0) continue;
        adjv[fill[eu[j]]] = ev[j], adje[fill[eu[j]]++] = j;
        adjv[fill[ev[j]]] = eu[j], adje[fill[ev[j]]++] = j;
    }

    /* loops (weight-0 columns) are the cycles of length 1 */
    int best = INT_MAX, bv = -1;
    int *wbuf = xmalloc((2 * (size_t)n + 2) * sizeof(int)), wn = 0;
    for (int j = 0; j < n && best > 1; j++)
        if (cnt[j] == 0) {
            u64 acc = 0;
            for (int w = 0; w < Wk; w++) acc |= Ann[(size_t)j * Wk + w];
            if (acc) best = 1, bv = 0, wbuf[0] = j, wn = 1;
        }
    int *stamp = xmalloc((size_t)V * sizeof(int)), *dist = xmalloc((size_t)V * sizeof(int)),
        *par = xmalloc((size_t)V * sizeof(int)), *qpos = xmalloc((size_t)V * sizeof(int)),
        *queue = xmalloc((size_t)V * sizeof(int));
    for (int u = 0; u < V; u++) stamp[u] = -1;
    u64 *A = xcalloc((size_t)V * Wk, sizeof(u64));
    long long nbfs = 0, ncand = 0;
    /* Lemma 7: a non-loop candidate has nominal length >= 2, so best <= 2 is final */
    for (int v = 0; v < V && best > 2; v++) {
        int R = best == INT_MAX ? INT_MAX : (best - 1) / 2;
        nbfs++;
        stamp[v] = v, dist[v] = 0, par[v] = -1, qpos[v] = 0, queue[0] = v;
        memset(A + (size_t)v * Wk, 0, (size_t)Wk * sizeof(u64));
        int qh = 0, qt = 1;
        while (qh < qt) {
            const int xv = queue[qh];
            if (dist[xv] > R) break;
            const u64 *Ax = A + (size_t)xv * Wk;
            for (int t = aptr[xv]; t < aptr[xv + 1]; t++) {
                const int y = adjv[t], e = adje[t];
                if (y < v) continue; /* BFS in the subgraph induced on vertices >= v */
                if (stamp[y] != v) {
                    if (dist[xv] < R) {
                        stamp[y] = v, dist[y] = dist[xv] + 1, par[y] = e, qpos[y] = qt, queue[qt++] = y;
                        xor3(A + (size_t)y * Wk, Ax, Ann + (size_t)e * Wk, Wk);
                    }
                    continue;
                }
                /* non-tree edge, handled once at its later-popped endpoint */
                if (qpos[y] >= qh || e == par[xv]) continue;
                const int l = dist[xv] + dist[y] + 1;
                if (l >= best) continue;
                ncand++;
                if (nonzero3(Ax, Ann + (size_t)e * Wk, A + (size_t)y * Wk, Wk)) {
                    best = l, bv = v;
                    R = (best - 1) / 2;
                    /* witness T_v(x) + e + T_v(y); parent pointers never change after discovery */
                    wn = 0;
                    wbuf[wn++] = e;
                    int ends[2] = {xv, y};
                    for (int s = 0; s < 2; s++)
                        for (int u = ends[s]; u != v;) {
                            const int pe = par[u];
                            wbuf[wn++] = pe;
                            u = eu[pe] == u ? ev[pe] : eu[pe];
                        }
                }
            }
            qh++;
        }
    }
    res->nbfs = nbfs, res->ncand = ncand;

    if (bv < 0) {
        res->status = ST_ERR;
        snprintf(res->reason, sizeof res->reason, "no nontrivial cycle found although k = %d", k);
    } else {
        u64 *c = xcalloc((size_t)Wn, sizeof(u64));
        for (int i = 0; i < wn; i++) flipbit(c, wbuf[i]);
        /* recheck: weight, c in ker Hg, L c != 0, and c not in rowsp Ho (independent of L) */
        int w = popc(c, Wn);
        int okker = in_kernel(Hg, c);
        int okL = 0;
        for (int t = 0; t < k && !okL; t++) {
            const u64 *row = L + (size_t)t * Wn;
            int s = 0;
            for (int q = 0; q < Wn; q++) s ^= __builtin_parityll(row[q] & c[q]);
            okL = s;
        }
        u64 *cc = xmalloc((size_t)Wn * sizeof(u64));
        memcpy(cc, c, (size_t)Wn * sizeof(u64));
        Basis Bo = ech_as_basis(Eo);
        int okrow = basis_reduce(&Bo, cc);
        free(cc);
        res->verified = (w == best) && okker && okL && okrow;
        res->d = best;
        res->wlen = w;
        res->wit = xmalloc((size_t)w * sizeof(int) + 4);
        int p = 0;
        for (int j = 0; j < n; j++)
            if (getbit(c, j)) res->wit[p++] = j;
        if (!res->verified) {
            res->status = ST_ERR;
            snprintf(res->reason, sizeof res->reason,
                     "witness recheck failed: weight %d vs %d, kernel %d, L %d, outside rowspace %d", w, best,
                     okker, okL, okrow);
        }
        free(c);
        if (count_limit > 0 && res->verified) {
            int ab;
            res->count = count_min(V, aptr, adjv, adje, Ann, Wk, best, n, cnt, count_limit, &res->cnodes, &ab);
            res->cstatus = ab ? 2 : 1;
        }
    }

    free(cnt), free(c1), free(c2), free(L), free(x), free(Ann), free(deg), free(ord), free(newid);
    free(eu), free(ev), free(vdeg), free(aptr), free(fill), free(adjv), free(adje);
    free(stamp), free(dist), free(par), free(qpos), free(queue), free(A), free(wbuf);
    res->secs = now() - t0;
}

static void print_side(const char *tag, const SideRes *r, int quiet)
{
    if (r->d >= 0) {
        printf("\"d%s\":%d,\"wit%s\":", tag, r->d, tag);
        if (quiet) printf("null");
        else {
            printf("[");
            for (int i = 0; i < r->wlen; i++) printf(i ? ",%d" : "%d", r->wit[i]);
            printf("]");
        }
        printf(",\"ver%s\":%s,\"lcheck%s\":%s,\"bfs%s\":%lld,\"cand%s\":%lld,\"secs%s\":%.3f", tag,
               r->verified ? "true" : "false", tag, r->lcheck ? "true" : "false", tag, r->nbfs, tag, r->ncand, tag,
               r->secs);
        if (r->cstatus == 1) printf(",\"count%s\":%lld,\"cnodes%s\":%lld", tag, r->count, tag, r->cnodes);
        else if (r->cstatus == 2) printf(",\"count%s\":null,\"cnodes%s\":%lld", tag, tag, r->cnodes);
    } else {
        printf("\"d%s\":null,\"wit%s\":null,\"secs%s\":%.3f", tag, tag, tag, r->secs);
    }
    if (r->status != ST_OK) printf(",\"reason%s\":\"%s\"", tag, r->reason);
}

int main(int argc, char **argv)
{
    const char *fx = NULL, *fz = NULL, *sidearg = "both";
    int quiet = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--side") && i + 1 < argc) sidearg = argv[++i];
        else if (!strcmp(argv[i], "--quiet")) quiet = 1;
        else if (!strcmp(argv[i], "--count")) count_limit = 1000000000LL;
        else if (!strcmp(argv[i], "--count-limit") && i + 1 < argc) count_limit = atoll(argv[++i]);
        else if (!fx) fx = argv[i];
        else if (!fz) fz = argv[i];
        else die("unexpected argument %s", argv[i]);
    }
    if (!fx || !fz) die("usage: graphlike_A NAME.hx NAME.hz [--side X|Z|both] [--quiet] [--count | --count-limit NODES]");
    int doX = 0, doZ = 0;
    if (!strcmp(sidearg, "X")) doX = 1;
    else if (!strcmp(sidearg, "Z")) doZ = 1;
    else if (!strcmp(sidearg, "both")) doX = doZ = 1;
    else die("--side must be X, Z or both");

    double t0 = now();
    Sp HX = read_sp(fx), HZ = read_sp(fz);
    if (HX.n != HZ.n) die("n differs: %d vs %d", HX.n, HZ.n);
    const int n = HX.n, Wn = (n + 63) / 64;

    /* commutation H_X H_Z^T = 0 */
    {
        int *mark = xmalloc((size_t)n * sizeof(int) + 4);
        for (int j = 0; j < n; j++) mark[j] = -1;
        for (int a = 0; a < HX.m; a++) {
            for (int t = HX.ptr[a]; t < HX.ptr[a + 1]; t++) mark[HX.idx[t]] = a;
            for (int b = 0; b < HZ.m; b++) {
                int s = 0;
                for (int t = HZ.ptr[b]; t < HZ.ptr[b + 1]; t++) s ^= (mark[HZ.idx[t]] == a);
                if (s) die("X-check %d and Z-check %d anticommute", a, b);
            }
        }
        free(mark);
    }

    Ech EX = echelon(&HX, Wn), EZ = echelon(&HZ, Wn);
    const int k = n - EX.rank - EZ.rank;
    if (k < 0) die("negative k");

    SideRes rX, rZ;
    memset(&rX, 0, sizeof rX);
    memset(&rZ, 0, sizeof rZ);
    /* d_Z: cycles of Gamma^+_X, i.e. ker H_X \ rowsp H_Z */
    rX.d = rZ.d = -1;
    if (doZ) run_side(&HX, &HZ, &EX, &EZ, n, k, "H_X", &rZ);
    else rZ.status = ST_ERR, snprintf(rZ.reason, sizeof rZ.reason, "not requested");
    /* d_X: cycles of Gamma^+_Z, i.e. ker H_Z \ rowsp H_X */
    if (doX) run_side(&HZ, &HX, &EZ, &EX, n, k, "H_Z", &rX);
    else rX.status = ST_ERR, snprintf(rX.reason, sizeof rX.reason, "not requested");

    printf("{\"n\":%d,\"k\":%d,\"mX\":%d,\"mZ\":%d,\"rankX\":%d,\"rankZ\":%d,", n, k, HX.m, HZ.m, EX.rank, EZ.rank);
    print_side("X", &rX, quiet);
    printf(",");
    print_side("Z", &rZ, quiet);
    printf(",\"secs\":%.3f}\n", now() - t0);
    int bad = (doX && rX.status == ST_ERR) || (doZ && rZ.status == ST_ERR);
    return bad ? 1 : 0;
}
