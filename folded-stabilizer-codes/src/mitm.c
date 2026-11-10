/* mitm.c: exhaustive meet-in-the-middle search for nontrivial logical operators of small Pauli weight
 * in a stabilizer code (verification track verify/fold2; own implementation).
 *
 * Every Pauli P with wt(P) = w <= a+b and zero syndrome is A*B with A = restriction of P to the first
 * min(w,a) qubits of supp P (wt <= a) and B = the rest (wt <= b), and syn(A) = syn(B).  So enumerating all A
 * of weight <= a against a table of all B of weight <= b with equal syndrome meets every element of
 * N(S) of weight <= a+b.  Each product P != I is tested for membership in S (reduction against the RREF of S).
 *
 * Input file: "n m" then m rows "nx x_1..x_nx nz z_1..z_nz" (generators of S, need not be independent).
 * Usage: mitm FILE a b      Output: minimal weight of a nontrivial logical found (<= a+b) or NONE.
 *        mitm FILE a b stab Output: every nonzero element of S of weight <= a+b (with repetitions).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 256
#define MAXW 4 /* 64-bit words per half (X or Z) */
#define MAXSW 4

static int n, mr, rk, a_max, b_max, SW, NW;
static uint64_t S[2 * MAXN][2 * MAXW]; /* rows, X words then Z words */
static int pivot[2 * MAXN];
static uint64_t syn[MAXN][4][MAXSW];

typedef struct {
  uint64_t key;
  uint64_t pk;
} entry;

static entry *T;
static size_t nT, capT;
static uint32_t *H;
static uint64_t Hmask;
static int hbits;
#define EMPTY 0xFFFFFFFFu

static inline uint64_t mkkey(const uint64_t *s) {
  uint64_t k = s[0];
  for (int i = 1; i < SW; i++) k ^= (s[i] + 0x632BE59BD9B4E019ULL * (uint64_t)i) * 0x9E3779B97F4A7C15ULL;
  return k;
}
static inline uint64_t hslot(uint64_t key) { return (key * 0x9E3779B97F4A7C15ULL) >> (64 - hbits); }

static int getbit(const uint64_t *v, int i) { return (v[i >> 6] >> (i & 63)) & 1; }

static void read_input(const char *fn) {
  FILE *f = fopen(fn, "r");
  if (!f) { perror(fn); exit(1); }
  if (fscanf(f, "%d %d", &n, &mr) != 2) exit(2);
  if (n > MAXN || mr > 2 * MAXN) { fprintf(stderr, "too large\n"); exit(2); }
  NW = (n + 63) / 64;
  for (int r = 0; r < mr; r++) {
    memset(S[r], 0, sizeof(S[r]));
    int c, q;
    if (fscanf(f, "%d", &c) != 1) exit(3);
    for (int j = 0; j < c; j++) { if (fscanf(f, "%d", &q) != 1) exit(3); S[r][q >> 6] ^= 1ULL << (q & 63); }
    if (fscanf(f, "%d", &c) != 1) exit(3);
    for (int j = 0; j < c; j++) { if (fscanf(f, "%d", &q) != 1) exit(3); S[r][NW + (q >> 6)] ^= 1ULL << (q & 63); }
  }
  fclose(f);
}

/* RREF over GF(2) on 2*NW words; bit index i<64*NW -> X part, else Z part */
static void rref(void) {
  int nb = 2 * NW * 64;
  rk = 0;
  for (int col = 0; col < nb && rk < mr; col++) {
    int sel = -1;
    for (int r = rk; r < mr; r++) if (getbit(S[r], col)) { sel = r; break; }
    if (sel < 0) continue;
    if (sel != rk) for (int w = 0; w < 2 * NW; w++) { uint64_t t = S[sel][w]; S[sel][w] = S[rk][w]; S[rk][w] = t; }
    for (int r = 0; r < mr; r++)
      if (r != rk && getbit(S[r], col)) for (int w = 0; w < 2 * NW; w++) S[r][w] ^= S[rk][w];
    pivot[rk] = col;
    rk++;
  }
}

static int in_S(const uint64_t *v0) {
  uint64_t v[2 * MAXW];
  memcpy(v, v0, sizeof(uint64_t) * 2 * NW);
  for (int r = 0; r < rk; r++)
    if (getbit(v, pivot[r])) for (int w = 0; w < 2 * NW; w++) v[w] ^= S[r][w];
  for (int w = 0; w < 2 * NW; w++) if (v[w]) return 0;
  return 1;
}

static void build_syn(void) {
  /* syndrome bit r of a single-qubit Pauli on q: symplectic product with basis row r.
     X_q: z_r(q); Z_q: x_r(q); Y_q: both. */
  SW = (rk + 63) / 64;
  if (SW > MAXSW) { fprintf(stderr, "rank too large\n"); exit(2); }
  memset(syn, 0, sizeof(syn));
  for (int q = 0; q < n; q++)
    for (int r = 0; r < rk; r++) {
      int x = getbit(S[r], q), z = getbit(S[r] + NW, q);
      if (z) syn[q][1][r >> 6] ^= 1ULL << (r & 63);
      if (x) syn[q][2][r >> 6] ^= 1ULL << (r & 63);
      if (x ^ z) syn[q][3][r >> 6] ^= 1ULL << (r & 63);
    }
}

/* ---- table of B ---- */
static int cq[16], cl[16];
static uint64_t cs[16][MAXSW];

static void tpush(const uint64_t *s, int depth) {
  if (nT == capT) {
    capT = capT ? 2 * capT : 1 << 20;
    T = realloc(T, capT * sizeof(entry));
    if (!T) { fprintf(stderr, "oom\n"); exit(4); }
  }
  uint64_t pk = 0;
  for (int j = 0; j < depth; j++) pk |= (uint64_t)((cq[j] << 2) | cl[j]) << (10 * j);
  T[nT].key = mkkey(s);
  T[nT].pk = pk;
  nT++;
}

static void tdfs(int depth, int start) {
  tpush(cs[depth], depth);
  if (depth == b_max) return;
  for (int q = start; q < n; q++)
    for (int L = 1; L <= 3; L++) {
      cq[depth] = q; cl[depth] = L;
      for (int w = 0; w < SW; w++) cs[depth + 1][w] = cs[depth][w] ^ syn[q][L][w];
      tdfs(depth + 1, q + 1);
    }
}

static int cmpe(const void *x, const void *y) {
  uint64_t a = ((const entry *)x)->key, b = ((const entry *)y)->key;
  return a < b ? -1 : a > b;
}

static void build_hash(void) {
  qsort(T, nT, sizeof(entry), cmpe);
  hbits = 1;
  while ((1ULL << hbits) < 2 * nT) hbits++;
  Hmask = (1ULL << hbits) - 1;
  H = malloc(sizeof(uint32_t) << hbits);
  if (!H) { fprintf(stderr, "oom\n"); exit(4); }
  memset(H, 0xFF, sizeof(uint32_t) << hbits);
  for (size_t i = 0; i < nT; i++) {
    if (i && T[i].key == T[i - 1].key) continue;
    uint64_t s = hslot(T[i].key);
    while (H[s] != EMPTY) s = (s + 1) & Hmask;
    H[s] = (uint32_t)i;
  }
}

/* ---- queries over A ---- */
static int best, dump_stab;
static uint64_t bestP[2 * MAXW];
static uint64_t nquery, nmatch, nstab, nlog;
static int aq[16], al[16];
static uint64_t as_[16][MAXSW];

static void handle(int depth, size_t i0, uint64_t key) {
  for (size_t i = i0; i < nT && T[i].key == key; i++) {
    uint64_t pk = T[i].pk;
    /* letters of P = A*B */
    int pq[32], pl[32], np = 0;
    int bq[8], bl[8], nb = 0;
    while (pk) { int g = pk & 1023; bq[nb] = g >> 2; bl[nb] = g & 3; nb++; pk >>= 10; }
    if (SW > 1) {
      uint64_t s[MAXSW];
      for (int w = 0; w < SW; w++) s[w] = as_[depth][w];
      for (int j = 0; j < nb; j++) for (int w = 0; w < SW; w++) s[w] ^= syn[bq[j]][bl[j]][w];
      int z = 1;
      for (int w = 0; w < SW; w++) if (s[w]) z = 0;
      if (!z) continue; /* hash-key collision */
    }
    nmatch++;
    int x = 0, y = 0;
    while (x < depth || y < nb) {
      if (y >= nb || (x < depth && aq[x] < bq[y])) { pq[np] = aq[x]; pl[np++] = al[x]; x++; }
      else if (x >= depth || bq[y] < aq[x]) { pq[np] = bq[y]; pl[np++] = bl[y]; y++; }
      else { int L = al[x] ^ bl[y]; if (L) { pq[np] = aq[x]; pl[np++] = L; } x++; y++; }
    }
    if (np == 0) continue;
    if (dump_stab) {
      uint64_t v[2 * MAXW];
      memset(v, 0, sizeof(v));
      for (int j = 0; j < np; j++) {
        if (pl[j] & 1) v[pq[j] >> 6] ^= 1ULL << (pq[j] & 63);
        if (pl[j] & 2) v[NW + (pq[j] >> 6)] ^= 1ULL << (pq[j] & 63);
      }
      if (in_S(v)) {
        printf("S");
        for (int j = 0; j < np; j++) printf(" %d:%c", pq[j], "?XZY"[pl[j]]);
        printf("\n");
      }
      continue;
    }
    if (np >= best) { /* still classify for statistics only when cheap: skip */ continue; }
    uint64_t v[2 * MAXW];
    memset(v, 0, sizeof(v));
    for (int j = 0; j < np; j++) {
      if (pl[j] & 1) v[pq[j] >> 6] ^= 1ULL << (pq[j] & 63);
      if (pl[j] & 2) v[NW + (pq[j] >> 6)] ^= 1ULL << (pq[j] & 63);
    }
    if (in_S(v)) { nstab++; continue; }
    nlog++;
    best = np;
    memcpy(bestP, v, sizeof(v));
  }
}

static void query(int depth) {
  nquery++;
  uint64_t key = mkkey(as_[depth]);
  uint64_t s = hslot(key);
  while (H[s] != EMPTY) {
    if (T[H[s]].key == key) { handle(depth, H[s], key); break; }
    s = (s + 1) & Hmask;
  }
}

static void adfs(int depth, int start) {
  query(depth);
  if (depth == a_max) return;
  for (int q = start; q < n; q++)
    for (int L = 1; L <= 3; L++) {
      aq[depth] = q; al[depth] = L;
      for (int w = 0; w < SW; w++) as_[depth + 1][w] = as_[depth][w] ^ syn[q][L][w];
      adfs(depth + 1, q + 1);
    }
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: mitm FILE a b\n"); return 1; }
  a_max = atoi(argv[2]);
  b_max = atoi(argv[3]);
  dump_stab = (argc > 4 && !strcmp(argv[4], "stab")); /* list all elements of S of weight <= a+b instead */
  if (b_max > 6 || a_max > 15) { fprintf(stderr, "b <= 6\n"); return 1; }
  clock_t t0 = clock();
  read_input(argv[1]);
  rref();
  build_syn();
  memset(cs[0], 0, sizeof(cs[0]));
  tdfs(0, 0);
  build_hash();
  memset(as_[0], 0, sizeof(as_[0]));
  best = a_max + b_max + 1;
  adfs(0, 0);
  double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
  printf("n=%d rank=%d k=%d a=%d b=%d table=%zu queries=%llu matches=%llu stab_hits=%llu logical_hits=%llu cpu=%.2f\n",
         n, rk, n - rk, a_max, b_max, nT, (unsigned long long)nquery, (unsigned long long)nmatch,
         (unsigned long long)nstab, (unsigned long long)nlog, secs);
  if (best <= a_max + b_max) {
    printf("FOUND weight %d:", best);
    printf(" X");
    for (int q = 0; q < n; q++) if (getbit(bestP, q)) printf(" %d", q);
    printf(" | Z");
    for (int q = 0; q < n; q++) if (getbit(bestP + NW, q)) printf(" %d", q);
    printf("\n");
  } else {
    printf("NONE of weight <= %d\n", a_max + b_max);
  }
  return 0;
}
