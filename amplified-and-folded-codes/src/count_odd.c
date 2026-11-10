/* Exact count of the nontrivial logicals of Pauli weight 2h+1 of a stabilizer code with d > 2h.
 *
 * Input: build.py .code file (n r L; generators; logical representatives spanning N(S)/S).
 * Notation as in mitm.c: key(p) = (<p,s_i>)_i, lab(p) = (<p,L_j>)_j.
 *
 * A Pauli p of weight 2h+1 with support q_1 < ... < q_{2h+1} is written uniquely as p = a + b with
 * a = p restricted to {q_1..q_h} (weight h) and b = p restricted to {q_{h+1}..q_{2h+1}} (weight h+1),
 * so max supp a < min supp b.  p is a nontrivial logical iff key(a) = key(b) and lab(a) != lab(b).
 * Conversely every such pair (a, b) gives a nontrivial logical of weight exactly 2h+1.  Hence the
 * number of pairs (a, b), |a| = h, |b| = h+1, max supp a < min supp b, key(a) = key(b),
 * lab(a) != lab(b) is the number of nontrivial logicals of weight 2h+1.
 * The hypothesis d > 2h is checked on the way: weight-h Paulis with equal keys must have equal labels.
 *
 * Usage: count_odd CODEFILE h [P] [MAXLIST]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 128
static int n, r, L, h, P, pass;
static uint64_t kv0[MAXN][4], kv1[MAXN][4];
static uint32_t lv[MAXN][4];
static uint64_t KEYMASK1;

static inline uint64_t mix(uint64_t a, uint64_t b) {
    uint64_t x = a ^ (b * 0x9E3779B97F4A7C15ULL);
    x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 29; x *= 0x94D049BB133111EBULL; x ^= x >> 32;
    return x;
}

typedef struct { uint64_t k0, k1; uint64_t enc; uint32_t lab; int32_t maxs; } ent_t;
static ent_t *tab; static uint64_t tmask; static size_t tcount;
static long label_conflicts;

static void tab_insert(uint64_t k0, uint64_t k1, uint32_t lab, int maxs, uint64_t enc) {
    uint64_t hsh = mix(k0 ^ 0x1234567ULL, k1) & tmask;
    while (tab[hsh].maxs >= 0) {
        if (tab[hsh].k0 == k0 && tab[hsh].k1 == k1 && tab[hsh].lab != lab) label_conflicts++;
        hsh = (hsh + 1) & tmask;
    }
    tab[hsh].k0 = k0; tab[hsh].k1 = k1; tab[hsh].lab = lab; tab[hsh].maxs = maxs; tab[hsh].enc = enc;
    tcount++;
}

static long found; static int maxlist;
static void report(const ent_t *a, uint64_t encb, int wb) {
    if (found >= maxlist) return;
    int xs[MAXN] = {0}, zs[MAXN] = {0};
    for (int i = 0; i < h; i++) { int f = (int)((a->enc >> (9 * i)) & 511), j = f >> 2, t = f & 3; xs[j] ^= t & 1; zs[j] ^= (t >> 1) & 1; }
    for (int i = 0; i < wb; i++) { int f = (int)((encb >> (9 * i)) & 511), j = f >> 2, t = f & 3; xs[j] ^= t & 1; zs[j] ^= (t >> 1) & 1; }
    printf("logical:");
    for (int j = 0; j < n; j++) if (xs[j] || zs[j]) printf(" %c%d", xs[j] && zs[j] ? 'Y' : (xs[j] ? 'X' : 'Z'), j);
    printf("\n");
}

/* build: weight-h Paulis in this pass */
static void build(int depth, int start, uint64_t k0, uint64_t k1, uint32_t lab, uint64_t enc, int last) {
    if (depth == h) {
        if ((int)(mix(k0, k1 & KEYMASK1) % (uint64_t)P) == pass) tab_insert(k0, k1 & KEYMASK1, lab, last, enc);
        return;
    }
    for (int j = start; j < n; j++)
        for (int t = 1; t <= 3; t++)
            build(depth + 1, j + 1, k0 ^ kv0[j][t], k1 ^ kv1[j][t], lab ^ lv[j][t], enc | ((uint64_t)((j << 2) | t) << (9 * depth)), j);
}

/* query: weight-(h+1) Paulis in this pass */
static int minsupp;
static void query(int depth, int start, uint64_t k0, uint64_t k1, uint32_t lab, uint64_t enc) {
    if (depth == h + 1) {
        uint64_t kk1 = k1 & KEYMASK1;
        if ((int)(mix(k0, kk1) % (uint64_t)P) != pass) return;
        uint64_t hsh = mix(k0 ^ 0x1234567ULL, kk1) & tmask;
        while (tab[hsh].maxs >= 0) {
            if (tab[hsh].k0 == k0 && tab[hsh].k1 == kk1 && tab[hsh].maxs < minsupp && tab[hsh].lab != lab) {
                report(&tab[hsh], enc, h + 1);
                found++;
            }
            hsh = (hsh + 1) & tmask;
        }
        return;
    }
    for (int j = start; j < n; j++) {
        if (depth == 0) minsupp = j;
        for (int t = 1; t <= 3; t++)
            query(depth + 1, j + 1, k0 ^ kv0[j][t], k1 ^ kv1[j][t], lab ^ lv[j][t], enc | ((uint64_t)((j << 2) | t) << (9 * depth)));
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s CODEFILE h [P] [MAXLIST]\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r");
    h = atoi(argv[2]); P = argc > 3 ? atoi(argv[3]) : 16; maxlist = argc > 4 ? atoi(argv[4]) : 1000;
    if (!f || fscanf(f, "%d %d %d", &n, &r, &L) != 3) return 1;
    if (n > MAXN || r > 104 || L > 32 || h > 6) { fprintf(stderr, "limits exceeded\n"); return 1; }
    KEYMASK1 = (r > 64) ? ((1ULL << (r - 64)) - 1) : 0;
    char *s = malloc(2 * n + 8);
    for (int i = 0; i < r + L; i++) {
        if (fscanf(f, "%s", s) != 1 || (int)strlen(s) != 2 * n) { fprintf(stderr, "bad row\n"); return 1; }
        for (int j = 0; j < n; j++) {
            int a = s[j] == '1', b = s[n + j] == '1';
            int bit[4] = {0, b, a, a ^ b};
            for (int t = 1; t <= 3; t++) {
                if (!bit[t]) continue;
                if (i < r) { if (i < 64) kv0[j][t] ^= 1ULL << i; else kv1[j][t] ^= 1ULL << (i - 64); }
                else lv[j][t] ^= 1u << (i - r);
            }
        }
    }
    fclose(f);
    /* number of weight-h Paulis */
    double cnt = 1; for (int i = 0; i < h; i++) cnt = cnt * (n - i) / (i + 1) * 3;
    uint64_t tsize = 1; while (tsize < 2.5 * cnt / P) tsize <<= 1;
    tab = malloc(tsize * sizeof(ent_t)); tmask = tsize - 1;
    clock_t t0 = clock();
    for (pass = 0; pass < P; pass++) {
        for (uint64_t i = 0; i < tsize; i++) tab[i].maxs = -1;
        tcount = 0;
        build(0, 0, 0, 0, 0, 0, -1);
        query(0, 0, 0, 0, 0, 0);
        fprintf(stderr, "pass %d/%d: table %zu entries, found so far %ld, label conflicts %ld, cpu %.1f s\n", pass + 1, P, tcount, found, label_conflicts, (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    printf("code %s n %d h %d: label conflicts among weight-%d Paulis %ld (must be 0, i.e. d > %d); nontrivial logicals of weight %d: %ld; cpu %.1f s\n",
           argv[1], n, h, h, label_conflicts, 2 * h, 2 * h + 1, found, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
