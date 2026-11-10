/* fewnony2: exact minimum Pauli weight of the nontrivial logicals with exactly t non-Y factors, t <= T.
 *
 * Input: build.py .code file ("n r L", r generator rows and L logical rows as 2n-bit strings, X part
 * then Z part; the L logicals together with the generators span N(S)).
 *
 * Every Pauli with non-Y support T (|T| = t) is p = (u + e | u) for exactly one u in F_2^n, e the
 * indicator of T: on j in T, p_j = X if u_j = 0 and Z if u_j = 1; off T, p_j = Y iff u_j = 1.
 * <p, (a_i | b_i)> = (u + e).b_i + u.a_i = u.(a_i + b_i) + e.b_i, so p in N(S) iff M u = B e
 * (M has rows a_i + b_i), and the solutions form u0 + ker M or are empty.  Pauli weight of p is
 * t + |u \ T|; p is in S iff <p, L_l> = u.(a_l + b_l) + e.b_l vanishes for all l.
 * (fewnony.c enumerates the same set 2^t times, splitting T into X- and Z-parts; this version
 * enumerates each p once.)
 *
 * Optional anchor bound A: only subsets T whose smallest element is < A are enumerated.  This is
 * exhaustive when a group of automorphisms of S moves some element of every T into {0..A-1}
 * (for the B510 fold: the translations h -> h + 5m and A = 5).
 *
 * Usage: fewnony2 CODEFILE T [A]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MW 8
typedef struct { uint64_t w[MW]; } bv;
static int n, r, L, NW, RW;

static inline int gb(const bv *v, int i) { return (v->w[i >> 6] >> (i & 63)) & 1; }
static inline void sb(bv *v, int i) { v->w[i >> 6] |= 1ULL << (i & 63); }
static inline void xv(bv *a, const bv *b, int W) { for (int i = 0; i < W; i++) a->w[i] ^= b->w[i]; }
static inline int dotp(const bv *a, const bv *b, int W) { int s = 0; for (int i = 0; i < W; i++) s ^= __builtin_popcountll(a->w[i] & b->w[i]) & 1; return s; }

static bv *A, *B, *LA, *LB;
static bv *Kall; static uint64_t *KlabAll; static long KN;
static int best[17]; static long nbest[17], nsolv[17], nsets[17];
static bv bestp[17];   /* u of a minimum, with T stored separately */
static bv bestT[17];
static uint64_t *keyX, *lX; static bv *uX;
static int T, Aanch;

static void process(int t, const bv *Tm, const bv *u0, uint64_t lab0) {
    bv notT; for (int i = 0; i < NW; i++) notT.w[i] = ~Tm->w[i];
    if (n % 64) notT.w[NW - 1] &= (1ULL << (n % 64)) - 1;
    for (long g = 0; g < KN; g++) {
        uint64_t lb = lab0 ^ KlabAll[g];
        if (!lb) continue;
        int wt = t;
        for (int i = 0; i < NW; i++) wt += __builtin_popcountll((u0->w[i] ^ Kall[g].w[i]) & notT.w[i]);
        if (wt < best[t]) { best[t] = wt; nbest[t] = 0; bestp[t] = *u0; xv(&bestp[t], &Kall[g], NW); bestT[t] = *Tm; }
        if (wt == best[t]) nbest[t]++;
    }
}

static void rec(int t, int depth, int start, bv *Tm, bv *u0, uint64_t key, uint64_t lab) {
    if (depth == t) {
        nsets[t]++;
        if (key == 0) { nsolv[t]++; process(t, Tm, u0, lab); }
        return;
    }
    int jmax = n - (t - depth);
    if (depth == 0 && Aanch > 0 && Aanch - 1 < jmax) jmax = Aanch - 1;
    for (int j = start; j <= jmax; j++) {
        bv T2 = *Tm, u2 = *u0;
        sb(&T2, j); xv(&u2, &uX[j], NW);
        rec(t, depth + 1, j + 1, &T2, &u2, key ^ keyX[j], lab ^ lX[j]);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s CODEFILE T [A]\n", argv[0]); return 1; }
    T = atoi(argv[2]); Aanch = argc > 3 ? atoi(argv[3]) : 0;
    if (T > 16) T = 16;
    FILE *f = fopen(argv[1], "r");
    if (!f || fscanf(f, "%d %d %d", &n, &r, &L) != 3) return 1;
    if (n > 64 * MW || r > 64 * MW || L > 64) { fprintf(stderr, "too large\n"); return 1; }
    NW = (n + 63) / 64; RW = (r + 63) / 64;
    A = calloc(r, sizeof(bv)); B = calloc(r, sizeof(bv)); LA = calloc(L, sizeof(bv)); LB = calloc(L, sizeof(bv));
    char *s = malloc(2 * n + 8);
    for (int i = 0; i < r + L; i++) {
        if (fscanf(f, "%s", s) != 1 || (int)strlen(s) != 2 * n) { fprintf(stderr, "bad row\n"); return 1; }
        for (int j = 0; j < n; j++) {
            if (s[j] == '1') sb(i < r ? &A[i] : &LA[i - r], j);
            if (s[n + j] == '1') sb(i < r ? &B[i] : &LB[i - r], j);
        }
    }
    fclose(f);
    bv *colB = calloc(n, sizeof(bv)), *colM = calloc(n, sizeof(bv));
    for (int i = 0; i < r; i++) for (int j = 0; j < n; j++) {
        if (gb(&B[i], j)) sb(&colB[j], i);
        if (gb(&A[i], j) ^ gb(&B[i], j)) sb(&colM[j], i);
    }
    uint64_t *labY = calloc(n, 8), *labX = calloc(n, 8);
    for (int l = 0; l < L; l++) for (int j = 0; j < n; j++) {
        if (gb(&LA[l], j) ^ gb(&LB[l], j)) labY[j] |= 1ULL << l;
        if (gb(&LB[l], j)) labX[j] |= 1ULL << l;
    }
    /* pivot columns of M, then [C | I_r] reduced: rows 0..rho-1 give alpha (u = sum alpha_l e_pc[l]),
       rows rho..r-1 a basis of the left kernel of M */
    int *pc = malloc(n * sizeof(int)), rho = 0;
    {
        bv *basis = calloc(n, sizeof(bv)); int *lead = malloc(n * sizeof(int));
        for (int j = 0; j < n; j++) {
            bv v = colM[j];
            for (int b = 0; b < rho; b++) if (gb(&v, lead[b])) xv(&v, &basis[b], RW);
            int ld = -1; for (int i = 0; i < r; i++) if (gb(&v, i)) { ld = i; break; }
            if (ld < 0) continue;
            basis[rho] = v; lead[rho] = ld; pc[rho] = j; rho++;
        }
        free(basis); free(lead);
    }
    int CW = (rho + 63) / 64;
    bv *Cr = calloc(r, sizeof(bv)), *Er = calloc(r, sizeof(bv));
    for (int i = 0; i < r; i++) { for (int l = 0; l < rho; l++) if (gb(&colM[pc[l]], i)) sb(&Cr[i], l); sb(&Er[i], i); }
    for (int col = 0, rk = 0; col < rho; col++) {
        int sel = -1;
        for (int i = rk; i < r; i++) if (gb(&Cr[i], col)) { sel = i; break; }
        if (sel < 0) { fprintf(stderr, "pivot error\n"); return 1; }
        bv t = Cr[sel]; Cr[sel] = Cr[rk]; Cr[rk] = t; t = Er[sel]; Er[sel] = Er[rk]; Er[rk] = t;
        for (int i = 0; i < r; i++) if (i != rk && gb(&Cr[i], col)) { xv(&Cr[i], &Cr[rk], CW); xv(&Er[i], &Er[rk], RW); }
        rk++;
    }
    int kap = r - rho;
    if (kap > 64) { fprintf(stderr, "left kernel too large\n"); return 1; }
    char *isp = calloc(n, 1); for (int l = 0; l < rho; l++) isp[pc[l]] = 1;
    int kdim = n - rho;
    if (kdim > 22) { fprintf(stderr, "ker M too large\n"); return 1; }
    bv *K = calloc(kdim, sizeof(bv)); uint64_t *Klab = calloc(kdim, 8);
    for (int j = 0, kk = 0; j < n; j++) if (!isp[j]) {
        sb(&K[kk], j);
        for (int l = 0; l < rho; l++) if (dotp(&Er[l], &colM[j], RW)) sb(&K[kk], pc[l]);
        kk++;
    }
    for (int q = 0; q < kdim; q++) {
        for (int i = 0; i < r; i++) { bv m = A[i]; xv(&m, &B[i], NW); if (dotp(&m, &K[q], NW)) { fprintf(stderr, "kernel error\n"); return 1; } }
        for (int j = 0; j < n; j++) if (gb(&K[q], j)) Klab[q] ^= labY[j];
    }
    KN = 1L << kdim;
    Kall = calloc(KN, sizeof(bv)); KlabAll = calloc(KN, 8);
    for (long g = 1; g < KN; g++) {     /* Kall[g] = sum of K[q] over the bits q of g */
        int q = __builtin_ctzl(g); long h = g & (g - 1);
        Kall[g] = Kall[h]; xv(&Kall[g], &K[q], NW); KlabAll[g] = KlabAll[h] ^ Klab[q];
    }
    keyX = calloc(n, 8); lX = calloc(n, 8); uX = calloc(n, sizeof(bv));
    for (int j = 0; j < n; j++) {
        for (int q = 0; q < kap; q++) if (dotp(&Er[rho + q], &colB[j], RW)) keyX[j] |= 1ULL << q;
        for (int l = 0; l < rho; l++) if (dotp(&Er[l], &colB[j], RW)) sb(&uX[j], pc[l]);
        for (int l = 0; l < rho; l++) if (gb(&uX[j], pc[l])) lX[j] ^= labY[pc[l]];
        lX[j] ^= labX[j];
    }
    printf("n %d r %d L %d: rank(A+B) %d, dim ker(A+B) %d, dim leftker %d, T %d, anchor bound %d\n", n, r, L, rho, kdim, kap, T, Aanch);
    fflush(stdout);
    clock_t t0 = clock();
    for (int t = 0; t <= T; t++) {
        best[t] = 1 << 30; nbest[t] = 0; nsolv[t] = 0; nsets[t] = 0;
        bv Tm, u0; memset(&Tm, 0, sizeof Tm); memset(&u0, 0, sizeof u0);
        rec(t, 0, 0, &Tm, &u0, 0, 0);
        printf("t=%d: subsets %ld, solvable %ld, ", t, nsets[t], nsolv[t]);
        if (best[t] == (1 << 30)) { printf("no nontrivial logical, cpu %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC); fflush(stdout); continue; }
        /* independent re-check of the stored minimum: commutation with every generator and a label */
        bv px = bestp[t], pz = bestp[t]; xv(&px, &bestT[t], NW);
        int comm = 1, wt = 0; uint64_t lb = 0;
        for (int i = 0; i < r; i++) if (dotp(&px, &B[i], NW) ^ dotp(&pz, &A[i], NW)) comm = 0;
        for (int l = 0; l < L; l++) if (dotp(&px, &LB[l], NW) ^ dotp(&pz, &LA[l], NW)) lb |= 1ULL << l;
        for (int i = 0; i < NW; i++) wt += __builtin_popcountll(px.w[i] | pz.w[i]);
        printf("min weight %d (%ld elements of that weight), recheck: commutes %d, label %s, weight %d; cpu %.1f s\n",
               best[t], nbest[t], comm, lb ? "nonzero" : "ZERO", wt, (double)(clock() - t0) / CLOCKS_PER_SEC);
        printf("  witness t=%d: ", t);
        for (int j = 0; j < n; j++) { int x = gb(&px, j), z = gb(&pz, j); if (x || z) printf("%c%d,", x && z ? 'Y' : (x ? 'X' : 'Z'), j); }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
