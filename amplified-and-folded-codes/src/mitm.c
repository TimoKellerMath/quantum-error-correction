/* Exact meet-in-the-middle test for the minimum distance of a stabilizer code (Pauli weight).
 *
 * Input file (build.py format): "n r L", then r stabilizer generators and L logical
 * representatives as strings of 2n bits (X part, then Z part).  The L logicals must span N(S)/S.
 *
 * For a Pauli p let key(p) = (<p,s_i>)_i (s_i the generators) and lab(p) = (<p,L_j>)_j.
 * For p in N(S): p in S iff lab(p) = 0 (S = N(S)^perp and N(S) = S + span(L_j)).
 * Hence for a, b of weight <= h: a+b is a nontrivial logical iff key(a) = key(b) and
 * lab(a) != lab(b); conversely every nontrivial logical of weight <= 2h is a+b with disjoint
 * supports and |a|,|b| <= h.  So min{wt(a+b) : key(a)=key(b), lab(a)!=lab(b)} equals d when
 * d <= 2h, and the absence of such pairs proves d > 2h.
 *
 * The table of all Paulis of weight <= h is processed in P passes, pass q keeping the records
 * whose key hashes to q (equal keys fall into the same pass), sorted by key and scanned.
 *
 * Usage: mitm CODEFILE h [P]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 128
typedef struct { uint64_t k0, k1, enc; } rec_t;   /* k1: key bits 64.., lab at 40..55, w at 56..59 */

static int n, r, L, h, P;
static uint64_t kv0[MAXN][4], kv1[MAXN][4];
static uint32_t lv[MAXN][4];
static rec_t *buf;
static size_t nbuf, cap;
static int pass;

static inline uint64_t mix(uint64_t a, uint64_t b) {
    uint64_t x = a ^ (b * 0x9E3779B97F4A7C15ULL);
    x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 29; x *= 0x94D049BB133111EBULL; x ^= x >> 32;
    return x;
}

static uint64_t KEYMASK1;

static void enumerate(int depth, int start, uint64_t k0, uint64_t k1, uint32_t lab, uint64_t enc) {
    /* record current Pauli (weight = depth) */
    if ((int)(mix(k0, k1 & KEYMASK1) % (uint64_t)P) == pass) {
        if (nbuf == cap) { cap = cap ? 2 * cap : (1 << 20); buf = realloc(buf, cap * sizeof(rec_t)); if (!buf) { perror("realloc"); exit(1); } }
        buf[nbuf].k0 = k0;
        buf[nbuf].k1 = (k1 & KEYMASK1) | ((uint64_t)lab << 40) | ((uint64_t)depth << 56);
        buf[nbuf].enc = enc;
        nbuf++;
    }
    if (depth == h) return;
    for (int j = start; j < n; j++)
        for (int t = 1; t <= 3; t++)
            enumerate(depth + 1, j + 1, k0 ^ kv0[j][t], k1 ^ kv1[j][t], lab ^ lv[j][t],
                      enc | ((uint64_t)((j << 2) | t) << (9 * depth)));
}

static int cmp(const void *a, const void *b) {
    const rec_t *x = a, *y = b;
    if (x->k0 != y->k0) return x->k0 < y->k0 ? -1 : 1;
    uint64_t u = x->k1 & KEYMASK1, v = y->k1 & KEYMASK1;
    if (u != v) return u < v ? -1 : 1;
    return 0;
}

typedef struct { uint64_t x[2], z[2]; } pauli_t;

static pauli_t decode(const rec_t *a) {
    pauli_t p; memset(&p, 0, sizeof p);
    int w = (int)((a->k1 >> 56) & 15);
    for (int i = 0; i < w; i++) {
        int f = (int)((a->enc >> (9 * i)) & 511), j = f >> 2, t = f & 3;
        if (t & 1) p.x[j >> 6] ^= 1ULL << (j & 63);
        if (t & 2) p.z[j >> 6] ^= 1ULL << (j & 63);
    }
    return p;
}

static int pweight(const pauli_t *p) {
    return __builtin_popcountll(p->x[0] | p->z[0]) + __builtin_popcountll(p->x[1] | p->z[1]);
}

static pauli_t *found; static size_t nfound, capfound; static int dmin = 1 << 30;

static void add_found(pauli_t p, int w) {
    if (w < dmin) { dmin = w; nfound = 0; }
    if (w > dmin) return;
    if (nfound == capfound) { capfound = capfound ? 2 * capfound : 1024; found = realloc(found, capfound * sizeof(pauli_t)); }
    found[nfound++] = p;
}

static int cmpp(const void *a, const void *b) { return memcmp(a, b, sizeof(pauli_t)); }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s CODEFILE h [P]\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 1; }
    h = atoi(argv[2]); P = argc > 3 ? atoi(argv[3]) : 16;
    if (fscanf(f, "%d %d %d", &n, &r, &L) != 3) return 1;
    if (n > MAXN || r > 104 || L > 16 || h > 7) { fprintf(stderr, "limits exceeded\n"); return 1; }
    KEYMASK1 = (r > 64) ? ((1ULL << (r - 64)) - 1) : 0;
    char *s = malloc(2 * n + 8);
    memset(kv0, 0, sizeof kv0); memset(kv1, 0, sizeof kv1); memset(lv, 0, sizeof lv);
    for (int i = 0; i < r + L; i++) {
        if (fscanf(f, "%s", s) != 1 || (int)strlen(s) != 2 * n) { fprintf(stderr, "bad row %d\n", i); return 1; }
        for (int j = 0; j < n; j++) {
            int a = s[j] == '1', b = s[n + j] == '1';
            /* <X_j, (a|b)> = b_j, <Z_j, (a|b)> = a_j, <Y_j, .> = a_j + b_j */
            int bit[4] = {0, b, a, a ^ b};
            for (int t = 1; t <= 3; t++) {
                if (!bit[t]) continue;
                if (i < r) { if (i < 64) kv0[j][t] ^= 1ULL << i; else kv1[j][t] ^= 1ULL << (i - 64); }
                else lv[j][t] ^= 1u << (i - r);
            }
        }
    }
    fclose(f);
    clock_t t0 = clock();
    size_t total = 0, ngroups_bad = 0, maxgroup = 0;
    for (pass = 0; pass < P; pass++) {
        nbuf = 0;
        enumerate(0, 0, 0, 0, 0, 0);
        total += nbuf;
        qsort(buf, nbuf, sizeof(rec_t), cmp);
        for (size_t i = 0; i < nbuf;) {
            size_t e = i + 1;
            while (e < nbuf && cmp(&buf[i], &buf[e]) == 0) e++;
            if (e - i > maxgroup) maxgroup = e - i;
            uint32_t l0 = (uint32_t)((buf[i].k1 >> 40) & 0xFFFF);
            int uniform = 1;
            for (size_t j = i + 1; j < e; j++) if ((uint32_t)((buf[j].k1 >> 40) & 0xFFFF) != l0) { uniform = 0; break; }
            if (!uniform) {
                ngroups_bad++;
                for (size_t a = i; a < e; a++)
                    for (size_t b = a + 1; b < e; b++) {
                        if (((buf[a].k1 ^ buf[b].k1) >> 40 & 0xFFFF) == 0) continue;
                        pauli_t pa = decode(&buf[a]), pb = decode(&buf[b]), p;
                        for (int q = 0; q < 2; q++) { p.x[q] = pa.x[q] ^ pb.x[q]; p.z[q] = pa.z[q] ^ pb.z[q]; }
                        add_found(p, pweight(&p));
                    }
            }
            i = e;
        }
        fprintf(stderr, "pass %d/%d: %zu records, cumulative %.1f s\n", pass + 1, P, nbuf, (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("code %s n %d r %d L %d h %d passes %d: %zu Paulis of weight <= %d, largest key group %zu, groups with distinct labels %zu, cpu %.1f s\n",
           argv[1], n, r, L, h, P, total, h, maxgroup, ngroups_bad, secs);
    if (dmin > 2 * h) {
        printf("RESULT: no nontrivial logical of Pauli weight <= %d, hence d >= %d\n", 2 * h, 2 * h + 1);
        return 0;
    }
    qsort(found, nfound, sizeof(pauli_t), cmpp);
    size_t uniq = 0;
    for (size_t i = 0; i < nfound; i++) if (i == 0 || memcmp(&found[i], &found[i - 1], sizeof(pauli_t))) found[uniq++] = found[i];
    printf("RESULT: d = %d exactly; %zu distinct nontrivial logicals of weight %d (as Pauli operators mod phase)\n", dmin, uniq, dmin);
    for (size_t i = 0; i < uniq && i < 3; i++) {
        printf("witness %zu: ", i);
        for (int j = 0; j < n; j++) {
            int x = found[i].x[j >> 6] >> (j & 63) & 1, z = found[i].z[j >> 6] >> (j & 63) & 1;
            if (x || z) printf("%c%d ", x && z ? 'Y' : (x ? 'X' : 'Z'), j);
        }
        printf("\n");
    }
    return 0;
}
