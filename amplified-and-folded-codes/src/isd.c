/* Information-set search for low Pauli-weight nontrivial logicals of a stabilizer code.
 *
 * Input file (build.py format): "n r L", r generators and L logical representatives as strings of
 * 2n bits (X part, then Z part); the L logicals must span N(S)/S.
 *
 * Vectors are stored interleaved: bit 2j = X part of qubit j, bit 2j+1 = Z part (in the current
 * frame).  The program computes a basis of N(S) = {p : <p,s_i> = 0} and, per basis vector, its
 * label (<p,L_j>)_j; p in N(S) lies in S iff its label is 0 (S = N(S)^perp).
 *
 * One trial: a random single-qubit Clifford frame (a random element of GL_2(F_2) = Sp_2(F_2) on the
 * two bits of each qubit; Pauli weight and symplectic form are invariant), a random pivot order
 * of the 2n bit columns (mode Q: random qubit order, both bits of a qubit consecutive; mode B:
 * uniformly random order of all 2n bits), reduced row echelon form of the basis of N(S) along that
 * order, then Lee-Brickell with p <= 2 (all rows and all sums of two rows).  Candidates with
 * nonzero label are scored by Pauli weight wt_P = #{j : x_j or z_j}.  Words with weight
 * <= best+SLACK are mapped back to the original frame, deduplicated, and their hits counted
 * (at most one hit per word per trial).
 *
 * Usage: isd CODEFILE TRIALS SEED FRAC_B OUTPREFIX
 *   FRAC_B in [0,1]: probability that a trial uses mode B (else mode Q).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXW 16          /* words per vector: n <= 512 */
#define SLACK 2
typedef struct { uint64_t v[MAXW]; } vec_t;

static int n, r, L, W, K;
static vec_t *Sg, *Lg, *Nb, *M;   /* generators, logicals, basis of N(S), work matrix */
static uint64_t *Nlab, *Mlab;

static uint64_t rs[2];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rnd(void) {      /* xoroshiro128+ */
    uint64_t s0 = rs[0], s1 = rs[1], res = s0 + s1;
    s1 ^= s0; rs[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); rs[1] = rotl(s1, 36);
    return res;
}
static uint64_t rndb(uint64_t m) { return rnd() % m; }

static inline int getb(const vec_t *a, int i) { return (a->v[i >> 6] >> (i & 63)) & 1; }
static inline void setb(vec_t *a, int i) { a->v[i >> 6] |= 1ULL << (i & 63); }
static inline void xorv(vec_t *a, const vec_t *b) { for (int w = 0; w < W; w++) a->v[w] ^= b->v[w]; }

static int sympl(const vec_t *a, const vec_t *b) {
    /* <a,b> = sum_j x_a z_b + z_a x_b; interleaved: swap bit pairs of b */
    int s = 0;
    for (int w = 0; w < W; w++) {
        uint64_t bs = ((b->v[w] & 0x5555555555555555ULL) << 1) | ((b->v[w] >> 1) & 0x5555555555555555ULL);
        s ^= __builtin_popcountll(a->v[w] & bs) & 1;
    }
    return s;
}
static inline int pweight(const vec_t *a) {
    int s = 0;
    for (int w = 0; w < W; w++) s += __builtin_popcountll((a->v[w] | (a->v[w] >> 1)) & 0x5555555555555555ULL);
    return s;
}

/* frame: per qubit an element of GL_2(F_2) as 2x2 matrix (a b; c d): x' = a x + b z, z' = c x + d z */
static const int GL2[6][4] = {{1,0,0,1},{0,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1},{1,1,1,0}};
static int INV[6];
static uint64_t mA[MAXW], mB[MAXW], mC[MAXW], mD[MAXW];  /* masks at even positions */
static uint64_t iA[MAXW], iB[MAXW], iC[MAXW], iD[MAXW];

static void set_frame(const int *g) {
    memset(mA, 0, sizeof mA); memset(mB, 0, sizeof mB); memset(mC, 0, sizeof mC); memset(mD, 0, sizeof mD);
    memset(iA, 0, sizeof iA); memset(iB, 0, sizeof iB); memset(iC, 0, sizeof iC); memset(iD, 0, sizeof iD);
    for (int j = 0; j < n; j++) {
        int b = 2 * j, w = b >> 6; uint64_t bit = 1ULL << (b & 63);
        const int *m = GL2[g[j]], *im = GL2[INV[g[j]]];
        if (m[0]) mA[w] |= bit; if (m[1]) mB[w] |= bit; if (m[2]) mC[w] |= bit; if (m[3]) mD[w] |= bit;
        if (im[0]) iA[w] |= bit; if (im[1]) iB[w] |= bit; if (im[2]) iC[w] |= bit; if (im[3]) iD[w] |= bit;
    }
}
static inline void apply(vec_t *a, const uint64_t *A, const uint64_t *B, const uint64_t *C, const uint64_t *D) {
    for (int w = 0; w < W; w++) {
        uint64_t x = a->v[w] & 0x5555555555555555ULL, z = (a->v[w] >> 1) & 0x5555555555555555ULL;
        uint64_t x2 = (A[w] & x) ^ (B[w] & z), z2 = (C[w] & x) ^ (D[w] & z);
        a->v[w] = x2 | (z2 << 1);
    }
}

/* ---- record of light words ---- */
typedef struct { vec_t p; int w; long hits; long first; long lasttrial; } rec_t;
static rec_t *R; static size_t nR, capR; static int best = 1 << 30;
static long trial_no;

static size_t find_rec(const vec_t *p) {
    for (size_t i = 0; i < nR; i++) if (!memcmp(&R[i].p, p, W * sizeof(uint64_t))) return i;
    return (size_t)-1;
}
static void record(const vec_t *pf, int w) {
    vec_t p = *pf;
    apply(&p, iA, iB, iC, iD);   /* back to original frame */
    if (w < best) {
        best = w;
        size_t k = 0;            /* drop records above new threshold */
        for (size_t i = 0; i < nR; i++) if (R[i].w <= best + SLACK) R[k++] = R[i];
        nR = k;
    }
    if (w > best + SLACK) return;
    size_t i = find_rec(&p);
    if (i == (size_t)-1) {
        if (nR == capR) { capR = capR ? 2 * capR : 256; R = realloc(R, capR * sizeof(rec_t)); }
        R[nR].p = p; R[nR].w = w; R[nR].hits = 1; R[nR].first = trial_no; R[nR].lasttrial = trial_no; nR++;
    } else if (R[i].lasttrial != trial_no) { R[i].hits++; R[i].lasttrial = trial_no; }
}

static int cmprec(const void *a, const void *b) {
    const rec_t *x = a, *y = b;
    if (x->w != y->w) return x->w - y->w;
    return (y->hits > x->hits) - (y->hits < x->hits);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s CODEFILE TRIALS SEED FRAC_B OUTPREFIX\n", argv[0]); return 1; }
    long trials = atol(argv[2]); uint64_t seed = strtoull(argv[3], 0, 10); double fracB = atof(argv[4]);
    rs[0] = seed * 0x9E3779B97F4A7C15ULL + 1; rs[1] = (seed + 7) * 0xD1B54A32D192ED03ULL + 3;
    for (int i = 0; i < 20; i++) rnd();
    for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) {
        const int *x = GL2[a], *y = GL2[b];
        int p0 = (x[0]*y[0] + x[1]*y[2]) & 1, p1 = (x[0]*y[1] + x[1]*y[3]) & 1, p2 = (x[2]*y[0] + x[3]*y[2]) & 1, p3 = (x[2]*y[1] + x[3]*y[3]) & 1;
        if (p0 == 1 && p1 == 0 && p2 == 0 && p3 == 1) INV[a] = b;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f || fscanf(f, "%d %d %d", &n, &r, &L) != 3) { fprintf(stderr, "bad input\n"); return 1; }
    W = (2 * n + 63) / 64;
    if (W > MAXW || L > 64) { fprintf(stderr, "too large\n"); return 1; }
    Sg = calloc(r, sizeof(vec_t)); Lg = calloc(L, sizeof(vec_t));
    char *s = malloc(2 * n + 8);
    for (int i = 0; i < r + L; i++) {
        if (fscanf(f, "%s", s) != 1 || (int)strlen(s) != 2 * n) { fprintf(stderr, "bad row\n"); return 1; }
        vec_t *t = i < r ? &Sg[i] : &Lg[i - r];
        for (int j = 0; j < n; j++) { if (s[j] == '1') setb(t, 2 * j); if (s[n + j] == '1') setb(t, 2 * j + 1); }
    }
    fclose(f);
    /* N(S): null space of x -> (<x,s_i>)_i.  <x,s> = x . swap(s); build rows swap(s_i), RREF, null space */
    vec_t *C = calloc(r, sizeof(vec_t));
    for (int i = 0; i < r; i++) for (int w = 0; w < W; w++) {
        uint64_t b = Sg[i].v[w];
        C[i].v[w] = ((b & 0x5555555555555555ULL) << 1) | ((b >> 1) & 0x5555555555555555ULL);
    }
    int *pivc = malloc(2 * n * sizeof(int)), rank = 0;
    char *ispiv = calloc(2 * n, 1);
    for (int col = 0; col < 2 * n && rank < r; col++) {
        int sel = -1;
        for (int i = rank; i < r; i++) if (getb(&C[i], col)) { sel = i; break; }
        if (sel < 0) continue;
        vec_t t = C[sel]; C[sel] = C[rank]; C[rank] = t;
        for (int i = 0; i < r; i++) if (i != rank && getb(&C[i], col)) xorv(&C[i], &C[rank]);
        pivc[rank++] = col; ispiv[col] = 1;
    }
    K = 2 * n - rank;
    Nb = calloc(K, sizeof(vec_t)); Nlab = calloc(K, sizeof(uint64_t));
    int kk = 0;
    for (int fcol = 0; fcol < 2 * n; fcol++) {
        if (ispiv[fcol]) continue;
        setb(&Nb[kk], fcol);
        for (int i = 0; i < rank; i++) if (getb(&C[i], fcol)) setb(&Nb[kk], pivc[i]);
        kk++;
    }
    /* checks: generators in N(S) (isotropy); basis vectors in N(S); labels */
    for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) if (sympl(&Sg[i], &Sg[j])) { fprintf(stderr, "not isotropic\n"); return 1; }
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < r; j++) if (sympl(&Nb[i], &Sg[j])) { fprintf(stderr, "null space error\n"); return 1; }
        for (int j = 0; j < L; j++) if (sympl(&Nb[i], &Lg[j])) Nlab[i] |= 1ULL << j;
    }
    printf("n %d, rank S %d, k %d, dim N(S) %d, trials %ld, seed %llu, fracB %.2f\n", n, rank, n - rank, K, trials, (unsigned long long)seed, fracB);
    M = calloc(K, sizeof(vec_t)); Mlab = calloc(K, sizeof(uint64_t));
    int *g = malloc(n * sizeof(int)), *order = malloc(2 * n * sizeof(int)), *qord = malloc(n * sizeof(int));
    clock_t t0 = clock();
    long first_best = -1; int prev_best = 1 << 30;
    for (trial_no = 1; trial_no <= trials; trial_no++) {
        for (int j = 0; j < n; j++) g[j] = (int)rndb(6);
        set_frame(g);
        for (int i = 0; i < K; i++) { M[i] = Nb[i]; apply(&M[i], mA, mB, mC, mD); Mlab[i] = Nlab[i]; }
        int modeB = (rnd() >> 11) * (1.0 / 9007199254740992.0) < fracB;
        if (modeB) {
            for (int i = 0; i < 2 * n; i++) order[i] = i;
            for (int i = 2 * n - 1; i > 0; i--) { int j = (int)rndb(i + 1), t = order[i]; order[i] = order[j]; order[j] = t; }
        } else {
            for (int i = 0; i < n; i++) qord[i] = i;
            for (int i = n - 1; i > 0; i--) { int j = (int)rndb(i + 1), t = qord[i]; qord[i] = qord[j]; qord[j] = t; }
            for (int i = 0; i < n; i++) { int sw = (int)(rnd() & 1); order[2 * i] = 2 * qord[i] + sw; order[2 * i + 1] = 2 * qord[i] + 1 - sw; }
        }
        int rk = 0;
        for (int oi = 0; oi < 2 * n && rk < K; oi++) {
            int col = order[oi], sel = -1;
            for (int i = rk; i < K; i++) if (getb(&M[i], col)) { sel = i; break; }
            if (sel < 0) continue;
            vec_t t = M[sel]; M[sel] = M[rk]; M[rk] = t;
            uint64_t tl = Mlab[sel]; Mlab[sel] = Mlab[rk]; Mlab[rk] = tl;
            int w = col >> 6; uint64_t bit = 1ULL << (col & 63);
            for (int i = 0; i < K; i++) if (i != rk && (M[i].v[w] & bit)) { xorv(&M[i], &M[rk]); Mlab[i] ^= Mlab[rk]; }
            rk++;
        }
        if (rk != K) { fprintf(stderr, "rank defect\n"); return 1; }
        int thr = best + SLACK;
        int wt[1024];
        for (int i = 0; i < K; i++) {
            wt[i] = pweight(&M[i]);
            if (Mlab[i] && wt[i] <= thr) { record(&M[i], wt[i]); thr = best + SLACK; }
        }
        for (int i = 0; i < K; i++) {
            for (int j = i + 1; j < K; j++) {
                if (!(Mlab[i] ^ Mlab[j])) continue;
                int s2 = 0;
                for (int w = 0; w < W; w++) { uint64_t v = M[i].v[w] ^ M[j].v[w]; s2 += __builtin_popcountll((v | (v >> 1)) & 0x5555555555555555ULL); }
                if (s2 <= thr) {
                    vec_t v; memset(&v, 0, sizeof v); for (int w = 0; w < W; w++) v.v[w] = M[i].v[w] ^ M[j].v[w];
                    record(&v, s2); thr = best + SLACK;
                }
            }
        }
        if (best < prev_best) { prev_best = best; first_best = trial_no;
            printf("trial %ld: new best weight %d (cpu %.1f s)\n", trial_no, best, (double)(clock() - t0) / CLOCKS_PER_SEC); fflush(stdout); }
        if (trial_no % 10000 == 0) { printf("trial %ld: best %d, cpu %.1f s\n", trial_no, best, (double)(clock() - t0) / CLOCKS_PER_SEC); fflush(stdout); }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    qsort(R, nR, sizeof(rec_t), cmprec);
    printf("DONE trials %ld cpu %.1f s (%.3f ms/trial); best weight %d first reached at trial %ld\n", trials, secs, 1000 * secs / trials, best, first_best);
    for (int w = best; w <= best + SLACK; w++) {
        long hits = 0, cnt = 0, firsthit = -1;
        for (size_t i = 0; i < nR; i++) if (R[i].w == w) { cnt++; hits += R[i].hits; if (firsthit < 0 || R[i].first < firsthit) firsthit = R[i].first; }
        printf("weight %d: distinct %ld, hits %ld, mean hits %.2f, first hit trial %ld\n", w, cnt, hits, cnt ? (double)hits / cnt : 0.0, firsthit);
    }
    char fn[1024]; snprintf(fn, sizeof fn, "%s.words", argv[5]);
    FILE *o = fopen(fn, "w");
    for (size_t i = 0; i < nR; i++) {
        fprintf(o, "%d %ld %ld ", R[i].w, R[i].hits, R[i].first);
        for (int j = 0; j < n; j++) { int x = getb(&R[i].p, 2 * j), z = getb(&R[i].p, 2 * j + 1); if (x || z) fprintf(o, "%c%d,", x && z ? 'Y' : (x ? 'X' : 'Z'), j); }
        fprintf(o, "\n");
    }
    fclose(o);
    return 0;
}
