/*
 * graphlike_B: exact distance of a CSS code on each graph-like side.
 *
 *   graphlike_B [--ref] NAME.hx NAME.hz
 *
 * Prints one JSON line {"n","k","dX","dZ","witX","witZ","secs", ...}.
 * Algorithm and correctness proof: NOTES.md (numbers like "Lemma 3" refer there).
 *
 * d_Z: graph G of H_X (vertices = X-checks + infinity, edge = qubit), logical
 * X-operators y_1..y_k (basis of ker H_Z modulo rowsp H_X); then
 *   d_Z = min_j min_v dist((v,0),(v,1)) in the double cover of G with voltage y_j
 * (Corollary 4, Prop. 8). All k covers are searched simultaneously, one bit per cover
 * (bit-sliced BFS, Lemma 10), with the vertex-elimination order of Prop. 9.
 * d_X: same with H_X and H_Z exchanged.
 */
#define _GNU_SOURCE
#include <errno.h>
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
    fprintf(stderr, "graphlike_B: ");
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

static void *xmalloc(size_t sz)
{
    void *p = malloc(sz ? sz : 1);
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
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---------------------------------------------------------------- input */

typedef struct {
    int n, m;
    int *ptr, *idx; /* CSR: row i is idx[ptr[i] .. ptr[i+1]-1] */
} SpMat;

static char *slurp(const char *path, size_t *len)
{
    FILE *f = fopen(path, "rb");
    if (!f) die("cannot open %s: %s", path, strerror(errno));
    size_t cap = 1 << 16, l = 0, r;
    char *buf = xmalloc(cap + 1);
    while ((r = fread(buf + l, 1, cap - l, f)) > 0) {
        l += r;
        if (l == cap) {
            cap *= 2;
            buf = realloc(buf, cap + 1);
            if (!buf) die("out of memory");
        }
    }
    if (ferror(f)) die("read error on %s", path);
    fclose(f);
    buf[l] = 0;
    *len = l;
    return buf;
}

static char *next_line(char **pp)
{
    char *s = *pp;
    if (*s == 0) return NULL;
    char *e = strchr(s, '\n');
    if (e) {
        *e = 0;
        *pp = e + 1;
    } else {
        *pp = s + strlen(s);
    }
    size_t L = strlen(s);
    if (L && s[L - 1] == '\r') s[L - 1] = 0;
    return s;
}

static int parse_ints(char *s, long *out, int maxcnt, const char *path, int lineno)
{
    int cnt = 0;
    for (;;) {
        while (*s == ' ' || *s == '\t') s++;
        if (*s == 0) break;
        char *q;
        errno = 0;
        long v = strtol(s, &q, 10);
        if (q == s || errno) die("%s:%d: bad integer", path, lineno);
        if (*q && *q != ' ' && *q != '\t') die("%s:%d: unexpected character '%c'", path, lineno, *q);
        if (cnt >= maxcnt) die("%s:%d: too many entries", path, lineno);
        out[cnt++] = v;
        s = q;
    }
    return cnt;
}

/* Format: first line "n m", then m lines with the 0-based supports of the checks. */
static SpMat read_checks(const char *path)
{
    size_t len;
    char *buf = slurp(path, &len);
    if (memchr(buf, 0, len)) die("%s: NUL byte in input", path);
    char *p = buf, *line;
    long hdr[3];
    if (!(line = next_line(&p))) die("%s: empty file", path);
    if (parse_ints(line, hdr, 3, path, 1) != 2) die("%s:1: header must be \"n m\"", path);
    if (hdr[0] < 0 || hdr[0] > 50000000 || hdr[1] < 0 || hdr[1] > 50000000)
        die("%s:1: n or m out of range", path);
    SpMat H;
    H.n = (int)hdr[0];
    H.m = (int)hdr[1];
    H.ptr = xmalloc(sizeof(int) * (H.m + 1));
    size_t cap = 1024, nnz = 0;
    H.idx = xmalloc(sizeof(int) * cap);
    long *tmp = xmalloc(sizeof(long) * (H.n + 1));
    int *mark = xmalloc(sizeof(int) * (H.n ? H.n : 1));
    for (int c = 0; c < H.n; c++) mark[c] = -1;
    H.ptr[0] = 0;
    for (int i = 0; i < H.m; i++) {
        if (!(line = next_line(&p))) die("%s: expected %d check lines, found %d", path, H.m, i);
        int cnt = parse_ints(line, tmp, H.n + 1, path, i + 2);
        if ((size_t)cnt + nnz > cap) {
            while ((size_t)cnt + nnz > cap) cap *= 2;
            H.idx = realloc(H.idx, sizeof(int) * cap);
            if (!H.idx) die("out of memory");
        }
        for (int t = 0; t < cnt; t++) {
            long c = tmp[t];
            if (c < 0 || c >= H.n) die("%s:%d: index %ld out of range [0,%d)", path, i + 2, c, H.n);
            if (mark[c] == i) die("%s:%d: repeated index %ld", path, i + 2, c);
            mark[c] = i;
            H.idx[nnz++] = (int)c;
        }
        H.ptr[i + 1] = (int)nnz;
    }
    int ln = H.m + 2;
    while ((line = next_line(&p))) {
        for (char *s = line; *s; s++)
            if (*s != ' ' && *s != '\t') die("%s:%d: trailing data after %d checks", path, ln, H.m);
        ln++;
    }
    free(tmp);
    free(mark);
    free(buf);
    return H;
}

/* ------------------------------------------------------ GF(2) linear algebra */

static inline int getbit(const u64 *v, int c) { return (int)((v[c >> 6] >> (c & 63)) & 1); }
static inline void flipbit(u64 *v, int c) { v[c >> 6] ^= 1ULL << (c & 63); }

static u64 *to_dense(const SpMat *H, int Wn)
{
    u64 *D = xcalloc((size_t)H->m * Wn, sizeof(u64));
    for (int i = 0; i < H->m; i++)
        for (int t = H->ptr[i]; t < H->ptr[i + 1]; t++) flipbit(D + (size_t)i * Wn, H->idx[t]);
    return D;
}

/* Reduced row echelon form of the rows R[0..rows-1] (pointers are permuted).
 * Returns the rank; pivcol[r] is the pivot column of row r. */
static int rref(u64 **R, int rows, int n, int Wn, int *pivcol)
{
    int r = 0;
    for (int c = 0; c < n && r < rows; c++) {
        int w = c >> 6;
        u64 bit = 1ULL << (c & 63);
        int p = -1;
        for (int i = r; i < rows; i++)
            if (R[i][w] & bit) {
                p = i;
                break;
            }
        if (p < 0) continue;
        u64 *t = R[p];
        R[p] = R[r];
        R[r] = t;
        const u64 *pr = R[r]; /* zero in all columns < c, so words < w may be skipped */
        for (int i = 0; i < rows; i++) {
            if (i == r || !(R[i][w] & bit)) continue;
            u64 *ri = R[i];
            for (int q = w; q < Wn; q++) ri[q] ^= pr[q];
        }
        pivcol[r++] = c;
    }
    return r;
}

/* Basis of ker H (vectors in F_2^n), from the RREF: one vector per free column. */
static int kernel_basis(const SpMat *H, int n, int Wn, u64 **out, int *rank_out)
{
    u64 *D = to_dense(H, Wn);
    u64 **R = xmalloc(sizeof(u64 *) * (H->m ? H->m : 1));
    for (int i = 0; i < H->m; i++) R[i] = D + (size_t)i * Wn;
    int *pivcol = xmalloc(sizeof(int) * (H->m + 1));
    int r = rref(R, H->m, n, Wn, pivcol);
    char *ispiv = xcalloc(n, 1);
    for (int i = 0; i < r; i++) ispiv[pivcol[i]] = 1;
    int nk = n - r;
    u64 *K = xcalloc((size_t)(nk ? nk : 1) * Wn, sizeof(u64));
    int t = 0;
    for (int f = 0; f < n; f++) {
        if (ispiv[f]) continue;
        u64 *z = K + (size_t)t * Wn;
        flipbit(z, f);
        for (int i = 0; i < r; i++)
            if (getbit(R[i], f)) flipbit(z, pivcol[i]);
        t++;
    }
    free(ispiv);
    free(pivcol);
    free(R);
    free(D);
    *out = K;
    *rank_out = r;
    return nk;
}

/* Incremental echelon basis; the pivot of a stored vector is its lowest set bit. */
typedef struct {
    int Wn, cnt;
    u64 *vec;
    int *piv; /* piv[c] = index of the stored vector with lowest bit c, or -1 */
} Ech;

static void ech_init(Ech *E, int n, int Wn)
{
    E->Wn = Wn;
    E->cnt = 0;
    E->vec = xcalloc((size_t)(n ? n : 1) * Wn, sizeof(u64));
    E->piv = xmalloc(sizeof(int) * (n ? n : 1));
    for (int c = 0; c < n; c++) E->piv[c] = -1;
}

static void ech_free(Ech *E)
{
    free(E->vec);
    free(E->piv);
}

/* Reduces v (in place) against E; inserts the remainder if it is nonzero.
 * Returns the index of the inserted vector, or -1 if v was in the span. */
static int ech_add(Ech *E, u64 *v)
{
    int Wn = E->Wn, w = 0;
    for (;;) {
        while (w < Wn && v[w] == 0) w++;
        if (w == Wn) return -1;
        int c = (w << 6) + __builtin_ctzll(v[w]);
        int b = E->piv[c];
        if (b < 0) {
            int id = E->cnt++;
            memcpy(E->vec + (size_t)id * Wn, v, sizeof(u64) * Wn);
            E->piv[c] = id;
            return id;
        }
        const u64 *bv = E->vec + (size_t)b * Wn; /* lowest bit of bv is c, so bv[q] = 0 for q < w */
        for (int q = w; q < Wn; q++) v[q] ^= bv[q];
    }
}

/* Basis of ker(B) modulo rowsp(A), assuming rowsp(A) is contained in ker(B) (Lemma 3). */
static int logical_basis(const SpMat *A, const SpMat *B, int n, int Wn, u64 **out, int *rA, int *rB)
{
    Ech E;
    ech_init(&E, n, Wn);
    u64 *tmp = xmalloc(sizeof(u64) * Wn);
    int ra = 0;
    for (int i = 0; i < A->m; i++) {
        memset(tmp, 0, sizeof(u64) * Wn);
        for (int t = A->ptr[i]; t < A->ptr[i + 1]; t++) flipbit(tmp, A->idx[t]);
        if (ech_add(&E, tmp) >= 0) ra++;
    }
    u64 *K;
    int rb, nk = kernel_basis(B, n, Wn, &K, &rb);
    u64 *L = xcalloc((size_t)(nk ? nk : 1) * Wn, sizeof(u64));
    int k = 0;
    for (int t = 0; t < nk; t++) {
        memcpy(tmp, K + (size_t)t * Wn, sizeof(u64) * Wn);
        int id = ech_add(&E, tmp);
        if (id >= 0) memcpy(L + (size_t)(k++) * Wn, E.vec + (size_t)id * Wn, sizeof(u64) * Wn);
    }
    free(K);
    free(tmp);
    ech_free(&E);
    *out = L;
    *rA = ra;
    *rB = rb;
    return k;
}

/* Plain rank by forward elimination on a private copy (used for the independent witness recheck). */
static int rank_rows(const u64 *rows, int m, int n, int Wn)
{
    u64 *D = xmalloc(sizeof(u64) * (size_t)(m ? m : 1) * Wn);
    memcpy(D, rows, sizeof(u64) * (size_t)m * Wn);
    int r = 0;
    for (int c = 0; c < n && r < m; c++) {
        int w = c >> 6;
        u64 bit = 1ULL << (c & 63);
        int p = -1;
        for (int i = r; i < m; i++)
            if (D[(size_t)i * Wn + w] & bit) {
                p = i;
                break;
            }
        if (p < 0) continue;
        if (p != r)
            for (int q = 0; q < Wn; q++) {
                u64 t = D[(size_t)p * Wn + q];
                D[(size_t)p * Wn + q] = D[(size_t)r * Wn + q];
                D[(size_t)r * Wn + q] = t;
            }
        for (int i = r + 1; i < m; i++)
            if (D[(size_t)i * Wn + w] & bit)
                for (int q = 0; q < Wn; q++) D[(size_t)i * Wn + q] ^= D[(size_t)r * Wn + q];
        r++;
    }
    free(D);
    return r;
}

/* ------------------------------------------------------------------- graph */

typedef struct {
    int V, E, inf;
    int *eu, *ev;            /* endpoints of edge (qubit) e; eu == ev for a loop */
    int *aptr, *aw, *ae;     /* CSR adjacency: entries (neighbour, edge) */
} Graph;

/* Graph of H (Lemma 1): returns 0 and a reason if some column has weight > 2. */
static int build_graph(const SpMat *H, Graph *G, char *reason, size_t rlen)
{
    int n = H->n;
    int *cw = xcalloc(n, sizeof(int)), *e0 = xmalloc(sizeof(int) * (n ? n : 1)), *e1 = xmalloc(sizeof(int) * (n ? n : 1));
    for (int i = 0; i < H->m; i++)
        for (int t = H->ptr[i]; t < H->ptr[i + 1]; t++) {
            int c = H->idx[t];
            if (cw[c] == 0) e0[c] = i;
            else if (cw[c] == 1) e1[c] = i;
            cw[c]++;
        }
    for (int c = 0; c < n; c++)
        if (cw[c] > 2) {
            snprintf(reason, rlen, "column %d of the check matrix has weight %d > 2", c, cw[c]);
            free(cw); free(e0); free(e1);
            return 0;
        }
    G->V = H->m + 1;
    G->inf = H->m;
    G->E = n;
    G->eu = xmalloc(sizeof(int) * (n ? n : 1));
    G->ev = xmalloc(sizeof(int) * (n ? n : 1));
    for (int c = 0; c < n; c++) {
        G->eu[c] = cw[c] >= 1 ? e0[c] : G->inf;
        G->ev[c] = cw[c] == 2 ? e1[c] : G->inf;
    }
    int *deg = xcalloc(G->V, sizeof(int));
    for (int c = 0; c < n; c++) {
        deg[G->eu[c]]++;
        if (G->ev[c] != G->eu[c]) deg[G->ev[c]]++;
    }
    G->aptr = xmalloc(sizeof(int) * (G->V + 1));
    G->aptr[0] = 0;
    for (int u = 0; u < G->V; u++) G->aptr[u + 1] = G->aptr[u] + deg[u];
    int tot = G->aptr[G->V];
    G->aw = xmalloc(sizeof(int) * (tot ? tot : 1));
    G->ae = xmalloc(sizeof(int) * (tot ? tot : 1));
    int *fill = xcalloc(G->V, sizeof(int));
    for (int c = 0; c < n; c++) {
        int u = G->eu[c], v = G->ev[c];
        G->aw[G->aptr[u] + fill[u]] = v;
        G->ae[G->aptr[u] + fill[u]++] = c;
        if (v != u) {
            G->aw[G->aptr[v] + fill[v]] = u;
            G->ae[G->aptr[v] + fill[v]++] = c;
        }
    }
    free(fill); free(deg); free(cw); free(e0); free(e1);
    return 1;
}

static void free_graph(Graph *G)
{
    free(G->eu); free(G->ev); free(G->aptr); free(G->aw); free(G->ae);
}

/* ------------------------------------------------------------------ search */

/*
 * Bit-sliced BFS (Prop. 9, Lemma 10). S is the E x Wk voltage matrix: bit j of
 * S[e] is y_j[e]. Roots are taken in the order 'order'; when root order[pos] is
 * processed, the roots order[0..pos-1] are deleted. Bit j of R1[u] after level t
 * is set iff dist_j((v,0),(u,1)) <= t in the double cover of the remaining
 * graph for voltage y_j. Returns the minimum, or INT_MAX if no odd closed walk.
 */
static int search_min(const Graph *G, const u64 *S, int k, int Wk, const int *order,
                      int *arg_pos, int *arg_j, long long *visits)
{
    int V = G->V;
    size_t sz = (size_t)V * Wk;
    u64 *R0 = xmalloc(sizeof(u64) * sz), *R1 = xmalloc(sizeof(u64) * sz);
    u64 *F0 = xcalloc(sz, sizeof(u64)), *F1 = xcalloc(sz, sizeof(u64));
    u64 *N0 = xcalloc(sz, sizeof(u64)), *N1 = xcalloc(sz, sizeof(u64));
    int *front = xmalloc(sizeof(int) * V), *nxt = xmalloc(sizeof(int) * V);
    char *innext = xcalloc(V, 1), *deleted = xcalloc(V, 1);
    u64 lastmask = (k & 63) ? ((1ULL << (k & 63)) - 1) : ~0ULL;
    int best = INT_MAX;
    long long vis = 0;
    *arg_pos = -1;
    *arg_j = -1;
    for (int pos = 0; pos < V; pos++) {
        int v = order[pos];
        memset(R0, 0, sizeof(u64) * sz);
        memset(R1, 0, sizeof(u64) * sz);
        for (int q = 0; q < Wk; q++) R0[(size_t)v * Wk + q] = F0[(size_t)v * Wk + q] = (q == Wk - 1) ? lastmask : ~0ULL;
        int nf = 0;
        front[nf++] = v;
        for (int t = 1; t < best; t++) {
            int nn = 0;
            for (int a = 0; a < nf; a++) {
                int u = front[a];
                const u64 *f0 = F0 + (size_t)u * Wk, *f1 = F1 + (size_t)u * Wk;
                for (int p = G->aptr[u]; p < G->aptr[u + 1]; p++) {
                    int w = G->aw[p];
                    if (deleted[w]) continue;
                    const u64 *s = S + (size_t)G->ae[p] * Wk;
                    u64 *n0 = N0 + (size_t)w * Wk, *n1 = N1 + (size_t)w * Wk;
                    /* lift of edge e from (u,a) ends in (w, a + y_j[e]) */
                    for (int q = 0; q < Wk; q++) {
                        u64 x0 = f0[q], x1 = f1[q], sq = s[q];
                        n0[q] |= (x0 & ~sq) | (x1 & sq);
                        n1[q] |= (x1 & ~sq) | (x0 & sq);
                    }
                    if (!innext[w]) {
                        innext[w] = 1;
                        nxt[nn++] = w;
                    }
                    vis++;
                }
            }
            for (int a = 0; a < nf; a++) {
                memset(F0 + (size_t)front[a] * Wk, 0, sizeof(u64) * Wk);
                memset(F1 + (size_t)front[a] * Wk, 0, sizeof(u64) * Wk);
            }
            nf = 0;
            for (int a = 0; a < nn; a++) {
                int w = nxt[a];
                innext[w] = 0;
                u64 *n0 = N0 + (size_t)w * Wk, *n1 = N1 + (size_t)w * Wk;
                u64 *r0 = R0 + (size_t)w * Wk, *r1 = R1 + (size_t)w * Wk;
                u64 *f0 = F0 + (size_t)w * Wk, *f1 = F1 + (size_t)w * Wk;
                u64 any = 0;
                for (int q = 0; q < Wk; q++) {
                    u64 y0 = n0[q] & ~r0[q], y1 = n1[q] & ~r1[q];
                    r0[q] |= y0;
                    r1[q] |= y1;
                    f0[q] = y0;
                    f1[q] = y1;
                    n0[q] = n1[q] = 0;
                    any |= y0 | y1;
                }
                if (any) front[nf++] = w;
            }
            const u64 *r1v = R1 + (size_t)v * Wk;
            int jj = -1;
            for (int q = 0; q < Wk; q++)
                if (r1v[q]) {
                    jj = (q << 6) + __builtin_ctzll(r1v[q]);
                    break;
                }
            if (jj >= 0) {
                best = t;
                *arg_pos = pos;
                *arg_j = jj;
                break;
            }
            if (nf == 0) break;
        }
        for (int a = 0; a < nf; a++) {
            memset(F0 + (size_t)front[a] * Wk, 0, sizeof(u64) * Wk);
            memset(F1 + (size_t)front[a] * Wk, 0, sizeof(u64) * Wk);
        }
        deleted[v] = 1;
    }
    free(R0); free(R1); free(F0); free(F1); free(N0); free(N1);
    free(front); free(nxt); free(innext); free(deleted);
    *visits = vis;
    return best;
}

/* Scalar BFS in the double cover for voltage y_j, restricted to non-deleted vertices.
 * Returns dist((root,0),(root,1)) and the projected walk (edge list), or -1. */
static int cover_bfs(const Graph *G, const u64 *S, int Wk, int j, int root, const char *deleted, int *walk)
{
    int NV = 2 * G->V;
    int *dist = xmalloc(sizeof(int) * NV), *pe = xmalloc(sizeof(int) * NV), *pv = xmalloc(sizeof(int) * NV);
    int *queue = xmalloc(sizeof(int) * NV);
    for (int x = 0; x < NV; x++) dist[x] = -1;
    int src = 2 * root, dst = 2 * root + 1, head = 0, tail = 0;
    dist[src] = 0;
    queue[tail++] = src;
    while (head < tail && dist[dst] < 0) {
        int x = queue[head++], u = x >> 1, a = x & 1;
        for (int p = G->aptr[u]; p < G->aptr[u + 1]; p++) {
            int w = G->aw[p];
            if (deleted && deleted[w]) continue;
            int e = G->ae[p];
            int y = 2 * w + (a ^ getbit(S + (size_t)e * Wk, j));
            if (dist[y] < 0) {
                dist[y] = dist[x] + 1;
                pe[y] = e;
                pv[y] = x;
                queue[tail++] = y;
            }
        }
    }
    int L = dist[dst];
    if (L >= 0 && walk) {
        int y = dst, i = L;
        while (y != src) {
            walk[--i] = pe[y];
            y = pv[y];
        }
    }
    free(dist); free(pe); free(pv); free(queue);
    return L;
}

/* Reference mode: for every j and every vertex, an unpruned scalar BFS without deletions. */
static int reference_min(const Graph *G, const u64 *S, int k, int Wk)
{
    int best = INT_MAX;
    for (int j = 0; j < k; j++)
        for (int v = 0; v < G->V; v++) {
            int L = cover_bfs(G, S, Wk, j, v, NULL, NULL);
            if (L >= 0 && L < best) best = L;
        }
    return best;
}

/* ------------------------------------------------------------------ driver */

typedef struct {
    int ok;              /* 1: distance computed; 0: see reason */
    int d, ref_d;
    int *wit, wlen;
    long long visits;
    double secs;
    char reason[256];
} SideResult;

static int cmp_deg_desc(const void *a, const void *b, void *ctx)
{
    const Graph *G = ctx;
    int u = *(const int *)a, v = *(const int *)b;
    int du = G->aptr[u + 1] - G->aptr[u], dv = G->aptr[v + 1] - G->aptr[v];
    if (du != dv) return dv - du;
    return u - v;
}

/*
 * Minimum weight of a vector in ker A outside rowsp B, where every column of A
 * has weight <= 2 (A is the graph side). For d_Z: A = H_X, B = H_Z.
 */
static SideResult solve_side(const SpMat *A, const SpMat *B, int n, int Wn, int kexp, int do_ref)
{
    SideResult res;
    memset(&res, 0, sizeof res);
    res.ref_d = -1;
    double t0 = now();
    Graph G;
    if (!build_graph(A, &G, res.reason, sizeof res.reason)) return res;

    /* logical operators of the opposite type: pairing test of Lemma 2 */
    u64 *Y;
    int rA, rB;
    int k = logical_basis(A, B, n, Wn, &Y, &rA, &rB);
    if (k != n - rA - rB) die("internal: logical count %d != n - rank A - rank B = %d", k, n - rA - rB);
    if (k != kexp) die("internal: logical count %d != k = %d from rank_rows", k, kexp);
    for (int j = 0; j < k; j++) /* each y_j must lie in ker B */
        for (int i = 0; i < B->m; i++) {
            int par = 0;
            for (int t = B->ptr[i]; t < B->ptr[i + 1]; t++) par ^= getbit(Y + (size_t)j * Wn, B->idx[t]);
            if (par) die("internal: logical %d not in the kernel", j);
        }
    if (k == 0) {
        snprintf(res.reason, sizeof res.reason, "k = 0: no logical operators");
        free(Y);
        free_graph(&G);
        return res;
    }
    int Wk = (k + 63) >> 6;
    u64 *S = xcalloc((size_t)n * Wk, sizeof(u64));
    for (int j = 0; j < k; j++)
        for (int e = 0; e < n; e++)
            if (getbit(Y + (size_t)j * Wn, e)) flipbit(S + (size_t)e * Wk, j);

    int *order = xmalloc(sizeof(int) * G.V);
    for (int u = 0; u < G.V; u++) order[u] = u;
    qsort_r(order, G.V, sizeof(int), cmp_deg_desc, &G);

    int pos, j;
    int d = search_min(&G, S, k, Wk, order, &pos, &j, &res.visits);
    if (d == INT_MAX) die("internal: k = %d > 0 but no odd closed walk found", k);

    /* witness (NOTES.md section 5): the walk realising d in cover j, with roots before pos deleted */
    char *deleted = xcalloc(G.V, 1);
    for (int p = 0; p < pos; p++) deleted[order[p]] = 1;
    int *walk = xmalloc(sizeof(int) * (size_t)(2 * G.V + 1)); /* a BFS path in the cover has < 2|V| edges */
    int L = cover_bfs(&G, S, Wk, j, order[pos], deleted, walk);
    if (L != d) die("internal: witness walk length %d != %d", L, d);
    u64 *x = xcalloc(Wn, sizeof(u64));
    for (int i = 0; i < L; i++) flipbit(x, walk[i]);

    /* independent recheck: weight, x in ker A, x not in rowsp B (rank test) */
    int wt = 0;
    for (int q = 0; q < Wn; q++) wt += __builtin_popcountll(x[q]);
    if (wt != d) die("witness check failed: weight %d != %d", wt, d);
    for (int i = 0; i < A->m; i++) {
        int par = 0;
        for (int t = A->ptr[i]; t < A->ptr[i + 1]; t++) par ^= getbit(x, A->idx[t]);
        if (par) die("witness check failed: violates check %d of the graph side", i);
    }
    u64 *Bx = to_dense(B, Wn);
    int rb0 = rank_rows(Bx, B->m, n, Wn);
    u64 *Bx1 = xmalloc(sizeof(u64) * (size_t)(B->m + 1) * Wn);
    memcpy(Bx1, Bx, sizeof(u64) * (size_t)B->m * Wn);
    memcpy(Bx1 + (size_t)B->m * Wn, x, sizeof(u64) * Wn);
    int rb1 = rank_rows(Bx1, B->m + 1, n, Wn);
    if (rb1 != rb0 + 1) die("witness check failed: witness lies in the row space of the other side");
    if (rb0 != rB) die("internal: rank mismatch %d != %d", rb0, rB);
    int pair = 0;
    for (int q = 0; q < Wn; q++) pair ^= __builtin_parityll(x[q] & Y[(size_t)j * Wn + q]);
    if (!pair) die("internal: witness has even pairing with logical %d", j);

    res.ok = 1;
    res.d = d;
    res.wit = xmalloc(sizeof(int) * (size_t)(d ? d : 1));
    res.wlen = 0;
    for (int e = 0; e < n; e++)
        if (getbit(x, e)) res.wit[res.wlen++] = e;
    if (do_ref) res.ref_d = reference_min(&G, S, k, Wk);
    res.secs = now() - t0;
    free(Bx); free(Bx1); free(x); free(walk); free(deleted); free(order); free(S); free(Y);
    free_graph(&G);
    return res;
}

static void css_check(const SpMat *HX, const SpMat *HZ, int Wn)
{
    u64 *DZ = to_dense(HZ, Wn);
    for (int i = 0; i < HX->m; i++)
        for (int l = 0; l < HZ->m; l++) {
            int par = 0;
            for (int t = HX->ptr[i]; t < HX->ptr[i + 1]; t++) par ^= getbit(DZ + (size_t)l * Wn, HX->idx[t]);
            if (par) die("not a CSS code: X-check %d and Z-check %d anticommute", i, l);
        }
    free(DZ);
}

static void print_d(const char *name, const SideResult *r)
{
    if (r->ok) printf(",\"d%s\":%d", name, r->d);
    else printf(",\"d%s\":null", name);
}

static void print_wit(const char *name, const SideResult *r)
{
    if (!r->ok) {
        printf(",\"wit%s\":null", name);
        return;
    }
    printf(",\"wit%s\":[", name);
    for (int i = 0; i < r->wlen; i++) printf(i ? ",%d" : "%d", r->wit[i]);
    printf("]");
}

static void print_extra(const char *name, const SideResult *r, int do_ref)
{
    if (!r->ok) printf(",\"reason%s\":\"%s\"", name, r->reason);
    else if (do_ref) printf(",\"ref%s\":%d", name, r->ref_d);
}

int main(int argc, char **argv)
{
    int do_ref = 0;
    const char *files[2];
    int nf = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--ref")) do_ref = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("usage: graphlike_B [--ref] NAME.hx NAME.hz\n");
            return 0;
        } else if (nf < 2) files[nf++] = argv[i];
        else die("too many arguments");
    }
    if (nf != 2) die("usage: graphlike_B [--ref] NAME.hx NAME.hz");
    double t0 = now();
    SpMat HX = read_checks(files[0]), HZ = read_checks(files[1]);
    if (HX.n != HZ.n) die("H_X has %d columns but H_Z has %d", HX.n, HZ.n);
    int n = HX.n, Wn = (n + 63) >> 6;
    css_check(&HX, &HZ, Wn);
    u64 *DX = to_dense(&HX, Wn), *DZ = to_dense(&HZ, Wn);
    int rX = rank_rows(DX, HX.m, n, Wn), rZ = rank_rows(DZ, HZ.m, n, Wn);
    free(DX);
    free(DZ);
    int k = n - rX - rZ;
    /* d_Z lives on the graph of H_X, d_X on the graph of H_Z */
    SideResult rz = solve_side(&HX, &HZ, n, Wn, k, do_ref);
    SideResult rx = solve_side(&HZ, &HX, n, Wn, k, do_ref);
    int bad = 0;
    if (do_ref && ((rx.ok && rx.ref_d != rx.d) || (rz.ok && rz.ref_d != rz.d))) bad = 1;
    printf("{\"n\":%d,\"k\":%d", n, k);
    print_d("X", &rx);
    print_d("Z", &rz);
    print_wit("X", &rx);
    print_wit("Z", &rz);
    printf(",\"secs\":%.3f", now() - t0);
    print_extra("X", &rx, do_ref);
    print_extra("Z", &rz, do_ref);
    printf("}\n");
    if (bad) {
        fprintf(stderr, "graphlike_B: reference mode disagrees with the main search\n");
        return 3;
    }
    return 0;
}
