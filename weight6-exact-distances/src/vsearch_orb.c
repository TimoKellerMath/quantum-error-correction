/*
 * vsearch_orb: exhaustive search for nontrivial Z-logical operators of weight <= W of a CSS code whose
 * X-checks all have weight 6. Implementation A of exact_distances.tex (Algorithm 1 there: anchored
 * search with exclusions, propagation and the majority rule); written independently of cluster2.c.
 *
 * Instance (written by codes.py): "n N k ns", then N lines of 6 role-ordered qubits of the X-checks
 * (role i of check h is the translate by h of role i of check 0; only mode -K uses the roles), ns lines
 * "w q.." (majority family, subsets of rowsp H_Z), k lines "w q.." (basis of ker H_Z modulo rowsp H_X).
 *
 * Modes: -K ORDER  pattern anchoring (31 jobs, ORDER = comma-separated hex role masks or "num"),
 *                  exact by Lemma 2.5 when the translations act simply transitively on the X-checks;
 *        -A        one anchor per block (qubit 0; then qubit N with block 1 excluded);
 *        -L FILE   jobs from FILE, one per line "anchor nex x_1 .. x_nex": the anchor qubit is put into
 *                  T and the listed qubits into X before the search. With anchors = representatives of
 *                  the orbits of a group of code automorphisms, and the qubits of the orbits of all
 *                  earlier anchors excluded, the jobs together are exact (Lemma 2.4).
 * Options: -W w, -e (enumerate all reported words instead of stopping at the first), -o file,
 *          -t tie (0: smallest index, 1: largest index among unsatisfied checks with fewest free qubits),
 *          -j job (run only this job), -M 0 (no majority rule), -B 0 (only ceil(u/3) bound).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, N, K, NS;
static int (*chk)[6];
static int qdeg[4096], qc[4096][8], qr[4096][8];
static int *soff, *sv, *sthr;
static int *qsoff, *qsv;
static uint64_t *lxm;

static unsigned char *st;          /* 0 free, 1 in T, 2 in X */
static unsigned char *inm, *outm;  /* role masks per check */
static int *ulist, *upos, nu;      /* unsatisfied checks */
static int *scnt;
static int Tsz;
static uint64_t cls;
static int *trail, tlen;
static int *cq, ncq, *sqq, nsq;    /* queues (stacks) of checks and saturated stabilizers */
static int conflict;
static unsigned char okT[4096], finT[4096], foutT[4096];
static int *mk, mkgen = 1;

static int W = 0, enumerate = 0, tie = 0, onlyjob = -1, use_maj = 1, use_pack = 1;
static long long nodes = 0, nfound = 0;
static int stopflag = 0;
static FILE *outf = NULL;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void die(const char *m) { fprintf(stderr, "error: %s\n", m); exit(1); }

static void toggle_unsat(int c) {
    if (__builtin_popcount(inm[c]) & 1) { upos[c] = nu; ulist[nu++] = c; }
    else { int p = upos[c], l = ulist[--nu]; ulist[p] = l; upos[l] = p; }
}

static void put_in(int q) {
    st[q] = 1; trail[tlen++] = q; Tsz++; cls ^= lxm[q];
    if (Tsz > W) conflict = 1;
    for (int k = 0; k < qdeg[q]; k++) { int c = qc[q][k]; inm[c] |= (unsigned char)(1u << qr[q][k]); toggle_unsat(c); cq[ncq++] = c; }
    if (use_maj)
        for (int i = qsoff[q]; i < qsoff[q + 1]; i++) {
            int s = qsv[i];
            scnt[s]++;
            if (scnt[s] > sthr[s]) conflict = 1;
            else if (scnt[s] == sthr[s]) sqq[nsq++] = s;
        }
}

static void put_out(int q) {
    st[q] = 2; trail[tlen++] = q;
    for (int k = 0; k < qdeg[q]; k++) { int c = qc[q][k]; outm[c] |= (unsigned char)(1u << qr[q][k]); cq[ncq++] = c; }
}

static void undo(int m) {
    while (tlen > m) {
        int q = trail[--tlen];
        if (st[q] == 1) {
            Tsz--; cls ^= lxm[q];
            for (int k = 0; k < qdeg[q]; k++) { int c = qc[q][k]; inm[c] &= (unsigned char)~(1u << qr[q][k]); toggle_unsat(c); }
            if (use_maj) for (int i = qsoff[q]; i < qsoff[q + 1]; i++) scnt[qsv[i]]--;
        } else {
            for (int k = 0; k < qdeg[q]; k++) { int c = qc[q][k]; outm[c] &= (unsigned char)~(1u << qr[q][k]); }
        }
        st[q] = 0;
    }
}

/* lower bound for |z \ T| over sets of free qubits meeting every unsatisfied check */
static int lbound(void) {
    if (nu == 0) return 0;
    int hmax = 0;
    for (int a = 0; a < nu; a++) {
        int c = ulist[a];
        for (int r = 0; r < 6; r++) {
            int q = chk[c][r];
            if (st[q]) continue;
            int h = 0;
            for (int k = 0; k < qdeg[q]; k++) h += __builtin_popcount(inm[qc[q][k]]) & 1;
            if (h > hmax) hmax = h;
        }
    }
    if (hmax == 0) return 1 << 20;
    int lb = (nu + hmax - 1) / hmax;
    if (use_pack) {
        mkgen++;
        int pk = 0;
        for (int a = 0; a < nu; a++) {
            int c = ulist[a], clash = 0;
            for (int r = 0; r < 6; r++) { int q = chk[c][r]; if (!st[q] && mk[q] == mkgen) { clash = 1; break; } }
            if (clash) continue;
            pk++;
            for (int r = 0; r < 6; r++) { int q = chk[c][r]; if (!st[q]) mk[q] = mkgen; }
        }
        if (pk > lb) lb = pk;
    }
    return lb;
}

static int propagate(void) {
    while (!conflict && (nsq > 0 || ncq > 0)) {
        if (nsq > 0) {
            int s = sqq[--nsq];
            if (scnt[s] > sthr[s]) { conflict = 1; break; }
            if (scnt[s] == sthr[s])
                for (int i = soff[s]; i < soff[s + 1]; i++) if (st[sv[i]] == 0) put_out(sv[i]);
            continue;
        }
        int c = cq[--ncq];
        int key = inm[c] | (outm[c] << 6);
        if (!okT[key]) { conflict = 1; break; }
        unsigned fi = finT[key], fo = foutT[key];
        for (int r = 0; r < 6 && !conflict; r++) {
            int q = chk[c][r];
            if (st[q]) continue;
            if (fi >> r & 1) put_in(q);
            else if (fo >> r & 1) put_out(q);
        }
    }
    int ok = !conflict && Tsz + lbound() <= W;
    conflict = 0; ncq = 0; nsq = 0;
    return ok;
}

static void build_tables(const unsigned char *forb) {
    for (int key = 0; key < 4096; key++) {
        int in = key & 63, out = key >> 6;
        okT[key] = finT[key] = foutT[key] = 0;
        if (in & out) continue;
        int any = 0, inter = 63, uni = 0;
        for (int P = 0; P < 64; P++) {
            if (__builtin_popcount(P) & 1) continue;
            if (P && forb[P]) continue;
            if ((P & in) != in || (P & out)) continue;
            any = 1; inter &= P; uni |= P;
        }
        if (!any) continue;
        okT[key] = 1;
        finT[key] = (unsigned char)(inter & ~in & 63);
        foutT[key] = (unsigned char)(~uni & ~in & ~out & 63);
    }
}

static void report(void) {
    nfound++;
    if (outf) {
        fprintf(outf, "%d :", Tsz);
        for (int q = 0; q < n; q++) if (st[q] == 1) fprintf(outf, " %d", q);
        fprintf(outf, "\n");
        fflush(outf);
    }
    if (!enumerate) stopflag = 1;
}

static void grow(void) {
    if (stopflag) return;
    nodes++;
    if (nu == 0) { if (cls) report(); return; }
    int best = -1, bf = 99;
    for (int a = 0; a < nu; a++) {
        int c = ulist[a], f = 0;
        for (int r = 0; r < 6; r++) if (!st[chk[c][r]]) f++;
        if (f < bf || (f == bf && (tie ? c > best : c < best))) { bf = f; best = c; }
    }
    int cand[6], nc = 0;
    for (int r = 0; r < 6; r++) if (!st[chk[best][r]]) cand[nc++] = chk[best][r];
    int m0 = tlen;
    for (int i = 0; i < nc && !stopflag; i++) {
        int q = cand[i];
        if (st[q] == 2) continue;
        if (st[q] == 1) { grow(); break; }
        int m1 = tlen;
        put_in(q);
        if (propagate()) grow();
        undo(m1);
        if (i < nc - 1) { put_out(q); if (!propagate()) break; }
    }
    undo(m0);
}

static int read_list(FILE *f, int *buf) {
    int w;
    if (fscanf(f, "%d", &w) != 1) die("bad list");
    for (int i = 0; i < w; i++) if (fscanf(f, "%d", &buf[i]) != 1) die("bad list");
    return w;
}

int main(int argc, char **argv) {
    if (argc < 2) die("usage: vsearch INST -W w (-K ORDER | -A) [-e] [-o out] [-t tie] [-j job] [-M 0] [-B 0]");
    const char *order = NULL;
    int anchormode = 0;
    const char *ofn = NULL;
    const char *jobfile = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-W")) W = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-K")) order = argv[++i];
        else if (!strcmp(argv[i], "-A")) anchormode = 1;
        else if (!strcmp(argv[i], "-L")) { jobfile = argv[++i]; anchormode = 2; }
        else if (!strcmp(argv[i], "-e")) enumerate = 1;
        else if (!strcmp(argv[i], "-o")) ofn = argv[++i];
        else if (!strcmp(argv[i], "-t")) tie = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-j")) onlyjob = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-M")) use_maj = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-B")) use_pack = atoi(argv[++i]);
        else die("unknown option");
    }
    if (!order && !anchormode) die("need -K, -A or -L");
    FILE *f = fopen(argv[1], "r");
    if (!f) die("cannot open instance");
    if (fscanf(f, "%d %d %d %d", &n, &N, &K, &NS) != 4) die("bad header");
    if (n != 2 * N || n > 4096 || K > 64) die("size");
    chk = malloc(sizeof(int[6]) * N);
    memset(qdeg, 0, sizeof qdeg);
    for (int c = 0; c < N; c++)
        for (int r = 0; r < 6; r++) {
            if (fscanf(f, "%d", &chk[c][r]) != 1) die("bad check");
            int q = chk[c][r];
            if (q < 0 || q >= n || qdeg[q] >= 8) die("bad qubit");
            qc[q][qdeg[q]] = c; qr[q][qdeg[q]] = r; qdeg[q]++;
        }
    int *buf = malloc(sizeof(int) * (n + 1));
    soff = malloc(sizeof(int) * (NS + 1)); sthr = malloc(sizeof(int) * NS);
    int cap = 16 * NS + 16; sv = malloc(sizeof(int) * cap);
    soff[0] = 0;
    for (int s = 0; s < NS; s++) {
        int w = read_list(f, buf);
        if (soff[s] + w > cap) { cap = 2 * (soff[s] + w); sv = realloc(sv, sizeof(int) * cap); }
        memcpy(sv + soff[s], buf, sizeof(int) * w);
        soff[s + 1] = soff[s] + w;
        sthr[s] = w / 2;
    }
    lxm = calloc(n, sizeof(uint64_t));
    for (int j = 0; j < K; j++) { int w = read_list(f, buf); for (int i = 0; i < w; i++) lxm[buf[i]] |= 1ULL << j; }
    fclose(f);
    /* qubit -> stabilizers */
    qsoff = calloc(n + 1, sizeof(int));
    for (int i = 0; i < soff[NS]; i++) qsoff[sv[i] + 1]++;
    for (int q = 0; q < n; q++) qsoff[q + 1] += qsoff[q];
    qsv = malloc(sizeof(int) * (soff[NS] + 1));
    int *pos = malloc(sizeof(int) * (n + 1)); memcpy(pos, qsoff, sizeof(int) * (n + 1));
    for (int s = 0; s < NS; s++) for (int i = soff[s]; i < soff[s + 1]; i++) qsv[pos[sv[i]]++] = s;
    st = calloc(n, 1); inm = calloc(N, 1); outm = calloc(N, 1);
    ulist = malloc(sizeof(int) * N); upos = malloc(sizeof(int) * N); nu = 0;
    scnt = calloc(NS + 1, sizeof(int));
    trail = malloc(sizeof(int) * (n + 1)); tlen = 0;
    cq = malloc(sizeof(int) * (32 * n + 64)); sqq = malloc(sizeof(int) * (NS + 64));
    mk = calloc(n, sizeof(int));
    if (ofn) outf = fopen(ofn, "w");

    double t0 = now();
    unsigned char forb[64];
    memset(forb, 0, sizeof forb);
    int jobs = 0;
    if (order) {
        int pats[31], np = 0;
        if (!strcmp(order, "num")) {
            pats[np++] = 63;
            for (int w = 4; w >= 2; w -= 2) for (int P = 1; P < 63; P++) if (__builtin_popcount(P) == w) pats[np++] = P;
        } else {
            char b2[512]; strncpy(b2, order, 511); b2[511] = 0;
            unsigned char seen[64] = {0};
            for (char *t = strtok(b2, ","); t; t = strtok(NULL, ",")) {
                int P = (int)strtol(t, NULL, 16);
                if (P <= 0 || P > 63 || (__builtin_popcount(P) & 1) || seen[P] || np >= 31) die("bad order");
                seen[P] = 1; pats[np++] = P;
            }
        }
        if (np != 31) die("order must contain the 31 nonempty even masks");
        for (int j = 0; j < 31 && !stopflag; j++) {
            if (j > 0) forb[pats[j - 1]] = 1;
            if (onlyjob >= 0 && j != onlyjob) continue;
            build_tables(forb);
            long long n0 = nodes, f0 = nfound; double tj = now();
            for (int r = 0; r < 6; r++) { int q = chk[0][r]; if (pats[j] >> r & 1) put_in(q); else put_out(q); }
            if (propagate()) grow();
            undo(0);
            jobs++;
            fprintf(stderr, "job %d pattern %02x: nodes %lld found %lld %.1f s\n", j, pats[j], nodes - n0, nfound - f0, now() - tj);
        }
    } else if (anchormode == 2) {
        build_tables(forb);
        FILE *jf = fopen(jobfile, "r");
        if (!jf) die("cannot open job file");
        int anc, nex, j = 0;
        while (!stopflag && fscanf(jf, "%d %d", &anc, &nex) == 2) {
            int *ex = malloc(sizeof(int) * (nex + 1));
            for (int i = 0; i < nex; i++) if (fscanf(jf, "%d", &ex[i]) != 1) die("bad job file");
            if (onlyjob >= 0 && j != onlyjob) { free(ex); j++; continue; }
            long long n0 = nodes, f0 = nfound; double tj = now();
            for (int i = 0; i < nex; i++) { if (ex[i] == anc) die("anchor excluded"); put_out(ex[i]); }
            put_in(anc);
            if (propagate()) grow();
            undo(0);
            jobs++;
            fprintf(stderr, "list job %d anchor %d excluded %d: nodes %lld found %lld %.1f s\n", j, anc, nex, nodes - n0, nfound - f0, now() - tj);
            free(ex); j++;
        }
        fclose(jf);
    } else {
        build_tables(forb);
        for (int j = 0; j < 2 && !stopflag; j++) {
            if (onlyjob >= 0 && j != onlyjob) continue;
            long long n0 = nodes, f0 = nfound; double tj = now();
            if (j == 0) put_in(0);
            else { for (int q = 0; q < N; q++) put_out(q); put_in(N); }
            if (propagate()) grow();
            undo(0);
            jobs++;
            fprintf(stderr, "anchor job %d: nodes %lld found %lld %.1f s\n", j, nodes - n0, nfound - f0, now() - tj);
        }
    }
    printf("VRESULT inst=%s W=%d mode=%s jobs=%d nodes=%lld found=%lld stop=%d time=%.2f maj=%d pack=%d tie=%d\n",
           argv[1], W, order ? order : (anchormode == 2 ? "list" : "anchor"), jobs, nodes, nfound, stopflag, now() - t0, use_maj, use_pack, tie);
    if (outf) fclose(outf);
    return 0;
}
