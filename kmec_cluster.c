// kmec_cluster.c — two-phase (partition, then reduce) exact counting of k-Markov equivalence classes
// of DAGs, for cluster runs at n = 9. Same counting method as kmec_iso.c and kmec_shard.c.
//
// Phase 1 (partition), one job-array task per part m of M:
//   ./kmec_cluster partition n P M m slices_per_part threads outdir
//   Generates slices m, m+M, m+2M, ... (S of them) out of M*S. A slice combines a geng graph split
//   (geng res/mod) with a directg orientation split (directg -s), which balances the load better
//   than either alone.
//   Threads pull slices from a counter. Each DAG is written as 9 uint16 parent masks (18 bytes) to
//   outdir/part<m>/b<bucket>.bin. The bucket is hash(class invariant) mod P, where the class
//   invariant is the sorted vertex invariants of the DAG's 2-signature.
//   Writes outdir/part<m>/DONE with the number of DAGs generated (A003087(n) summed over all parts).
//
// Phase 2 (reduce), one task per bucket b (optionally a sub-bucket j of m_sub):
//   ./kmec_cluster reduce n kmin kmax P b msub j threads indir
//   Requires all M parts to be complete (DONE files). Reads indir/part*/b<b>.bin, recomputes each
//   DAG's 2-signature hash and checks that it belongs to bucket b. With msub > 1, keeps only DAGs whose
//   kmin-signature hash mod msub == j. The kmin-signature is a class invariant for every k >= kmin and
//   is much finer than the 2-signature, so large buckets split evenly. msub = 1, j = 0 keeps all.
//   Canonicalizes the k-signatures (k = kmin..kmax, kmin >= 2) and prints, for each k,
//   "n k unlabeled labeled" plus a header with the number of DAGs read and kept.
//
// Correctness:
//  * Every DAG in a k-class (k >= 2) has the same 2-signature up to isomorphism, hence the same
//    bucket; for k >= kmin it also has the same kmin-signature, hence the same sub-bucket. So buckets
//    and sub-buckets are disjoint unions of classes, and counts add up. Bucket b is exactly the union
//    of its sub-buckets j = 0..msub-1.
//  * Canonical forms are exact: lexicographic minimum over relabelings within invariant cells;
//    |Aut| = number of minimizing relabelings; labeled count = sum of n!/|Aut|.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <stdatomic.h>
#include <sys/stat.h>
#include <glob.h>
#include <fcntl.h>
#include <unistd.h>

#define MAXN 9
#define MAXS 512
#define NSHARD 4096
#define KB 2

static int n, KMIN, KMAX, NT;
static int npairs, pidx[MAXN][MAXN], pu[40], pv[40], MB;
static int nS; static unsigned Slist[MAXS]; static int Sidx[1 << MAXN]; static int sizeEnd[MAXN];
static uint64_t nfact = 1;
static const uint64_t A003087[10] = {1,1,2,6,31,302,5984,243668,20286025,3424938010ULL};

static inline uint64_t mix64(uint64_t z){ z += 0x9e3779b97f4a7c15ULL; z = (z ^ (z>>30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z>>27)) * 0x94d049bb133111ebULL; return z ^ (z>>31); }
static inline unsigned orset(const unsigned *a, unsigned X){ unsigned r=0; while(X){ int i=__builtin_ctz(X); X&=X-1; r|=a[i]; } return r; }
static void die(const char *m){ fprintf(stderr, "ERROR: %s\n", m); exit(2); }
static void *xcalloc(size_t a, size_t b){ void *p = calloc(a, b); if (!p) die("out of memory"); return p; }

static void setup(void){
    nfact = 1; for (int i=2;i<=n;i++) nfact *= i;
    npairs=0; for (int u=0;u<n;u++) for (int v=u+1;v<n;v++){ pidx[u][v]=pidx[v][u]=npairs; pu[npairs]=u; pv[npairs]=v; npairs++; }
    MB = (npairs + 7) / 8;
    nS=0; for (int s=0;s<=n-2;s++){ for (unsigned S=0;S<(1u<<n);S++) if (__builtin_popcount(S)==s){ Sidx[S]=nS; Slist[nS++]=S; } sizeEnd[s]=nS; }
}

// ---------- d-separation ----------
static inline unsigned bayesball(const unsigned *par, const unsigned *ch, int u, unsigned S, unsigned A){
    unsigned visU = 1u<<u, visD = 0, newU = visU, newD = 0;
    while (newU | newD) {
        unsigned upOpen = newU & ~S, dnChain = newD & ~S, dnColl = newD & A;
        unsigned nu = orset(par, upOpen | dnColl), nd = orset(ch, upOpen | dnChain);
        newU = nu & ~visU; newD = nd & ~visD; visU |= newU; visD |= newD;
    }
    return (visU | visD) & ~S;
}
typedef struct { unsigned par[MAXN], ch[MAXN], anc[MAXN]; } Dag;
static void dag_prepare(Dag *d){
    memset(d->ch, 0, sizeof d->ch);
    for (int x=0;x<n;x++){ unsigned p=d->par[x]; while(p){ int i=__builtin_ctz(p); p&=p-1; d->ch[i]|=1u<<x; } }
    for (int x=0;x<n;x++){ unsigned a=1u<<x, prev=0; while (a!=prev){ prev=a; a |= orset(d->par, a); } d->anc[x]=a; }
}
static void compute_sep(const Dag *d, uint64_t *sep, int from, int to){
    for (int i=from;i<to;i++){
        unsigned S = Slist[i], A = orset(d->anc, S) | S; uint64_t m = 0;
        for (int u=0;u<n;u++){ if (S>>u&1) continue; unsigned R = bayesball(d->par, d->ch, u, S, A);
            for (int v=u+1; v<n; v++) if (!(S>>v&1) && !(R>>v&1)) m |= 1ULL<<pidx[u][v]; }
        sep[i] = m;
    }
}

// ---------- invariants ----------
static int cmp_u64(const void *a, const void *b){ uint64_t x=*(const uint64_t*)a, y=*(const uint64_t*)b; return (x>y)-(x<y); }
static void vertex_invariants(const uint64_t *sep, int k, uint64_t *inv){
    int L = sizeEnd[k];
    uint32_t pairc[MAXN][MAXN][MAXN]; memset(pairc, 0, sizeof pairc);
    uint32_t inS[MAXN][MAXN]; memset(inS, 0, sizeof inS);
    for (int i=0;i<L;i++){
        unsigned S = Slist[i]; int s = __builtin_popcount(S);
        uint64_t m = sep[i]; int cntm = __builtin_popcountll(m);
        while (m){ int p=__builtin_ctzll(m); m&=m-1; pairc[pu[p]][pv[p]][s]++; pairc[pv[p]][pu[p]][s]++; }
        unsigned T=S; while (T){ int b=__builtin_ctz(T); T&=T-1; inS[b][s] += cntm; }
    }
    for (int x=0;x<n;x++){
        uint64_t h = 0xabcdefULL;
        for (int s=0;s<=k;s++){ uint32_t c1=0; for (int y=0;y<n;y++) if (y!=x) c1 += pairc[x][y][s]; h = mix64(h ^ c1 ^ ((uint64_t)inS[x][s]<<32) ^ ((uint64_t)s<<56)); }
        inv[x] = h;
    }
    uint64_t inv2[MAXN];
    for (int round=0; round<2; round++){
        for (int x=0;x<n;x++){
            uint64_t lst[MAXN]; int L2=0;
            for (int y=0;y<n;y++){ if (y==x) continue; uint64_t h = mix64(inv[y]); for (int s=0;s<=k;s++) h = mix64(h ^ pairc[x][y][s] ^ ((uint64_t)s<<48)); lst[L2++]=h; }
            qsort(lst, L2, sizeof lst[0], cmp_u64);
            uint64_t h = mix64(inv[x] ^ 0x5555ULL); for (int i=0;i<L2;i++) h = mix64(h ^ lst[i]);
            inv2[x] = h;
        }
        memcpy(inv, inv2, sizeof(uint64_t)*n);
    }
}
// class hash from the 2-signature (sep must hold at least sizeEnd[KB] entries)
static uint64_t class_hash_k(const uint64_t *sep, int k){
    uint64_t inv[MAXN]; vertex_invariants(sep, k, inv);
    qsort(inv, n, sizeof inv[0], cmp_u64);
    uint64_t h = 0x77ULL ^ ((uint64_t)k << 32); for (int i=0;i<n;i++) h = mix64(h ^ inv[i]);
    return h;
}
static uint64_t class_hash(const uint64_t *sep){ return class_hash_k(sep, KB); }

// =====================================================================================
// Phase 1: partition
// =====================================================================================
typedef struct { FILE *f; pthread_mutex_t mu; } BFile;
static BFile *bfiles; static int P_;
static atomic_int g_slice_next; static int g_part, g_nparts, g_slice_count, g_slice_mod;
static atomic_uint_fast64_t g_generated;

// slice s of mod: geng split G = 64 (when mod is a multiple of 64) times directg split D = mod/G
static void slice_cmd(char *cmd, size_t len, int s, int mod){
    int G = (mod % 64 == 0) ? 64 : 1, D = mod / G;
    snprintf(cmd, len, "bash -o pipefail -c 'geng -q %d %d/%d | directg -q -a -T -s%d/%d'", n, s % G, G, s / G, D);
}
// read one DAG line; returns 1 on a DAG, 0 at a clean end of stream; dies on anything else
static int read_dag(FILE *in, Dag *d){
    int nv, ne, rc = fscanf(in, "%d %d", &nv, &ne);
    if (rc == EOF) { if (ferror(in)) die("read error"); return 0; }
    if (rc != 2 || nv != n || ne < 0 || ne > npairs) die("bad input line");
    memset(d->par, 0, sizeof d->par);
    for (int e=0;e<ne;e++){ int a,b; if (fscanf(in, "%d %d",&a,&b)!=2 || a<0||a>=n||b<0||b>=n) die("parse error"); d->par[b] |= 1u<<a; }
    return 1;
}
#define WBUF 4096
static void *part_worker(void *arg){
    (void)arg;
    uint16_t (*buf)[WBUF][MAXN] = xcalloc((size_t)P_, sizeof *buf);
    int *cnt = xcalloc((size_t)P_, sizeof(int));
    uint64_t sep[MAXS]; uint64_t gen = 0;
    for (;;){
        int s = atomic_fetch_add(&g_slice_next, 1);
        if (s >= g_slice_count) break;
        char cmd[256];
        int res = g_part + g_nparts * s;   /* interleaved: part m takes slices m, m+M, m+2M, ... */
        slice_cmd(cmd, sizeof cmd, res, g_slice_mod);
        FILE *in = popen(cmd, "r"); if (!in) die("popen failed");
        Dag d;
        while (read_dag(in, &d)){
            dag_prepare(&d); compute_sep(&d, sep, 0, sizeEnd[KB]);
            int b = (int)(class_hash(sep) % (uint64_t)P_);
            for (int x=0;x<MAXN;x++) buf[b][cnt[b]][x] = (uint16_t)(x<n ? d.par[x] : 0);
            if (++cnt[b] == WBUF){ pthread_mutex_lock(&bfiles[b].mu); if (fwrite(buf[b], sizeof buf[b][0], WBUF, bfiles[b].f) != WBUF) die("write failed"); pthread_mutex_unlock(&bfiles[b].mu); cnt[b] = 0; }
            gen++;
        }
        int st = pclose(in);
        if (st != 0) { fprintf(stderr, "ERROR: nauty slice %d/%d exited with status %d\n", res, g_slice_mod, st); exit(2); }
    }
    for (int b=0;b<P_;b++) if (cnt[b]){ pthread_mutex_lock(&bfiles[b].mu); if (fwrite(buf[b], sizeof buf[b][0], (size_t)cnt[b], bfiles[b].f) != (size_t)cnt[b]) die("write failed"); pthread_mutex_unlock(&bfiles[b].mu); }
    atomic_fetch_add(&g_generated, gen);
    free(buf); free(cnt); return NULL;
}
static int do_partition(int argc, char **argv){
    if (argc < 9) die("usage: partition n P M m slices_per_part threads outdir");
    n = atoi(argv[2]); P_ = atoi(argv[3]); int M = atoi(argv[4]), m = atoi(argv[5]), S = atoi(argv[6]); NT = atoi(argv[7]); const char *out = argv[8];
    if (n < 4 || n > MAXN || P_ < 1 || P_ > 65536 || M < 1 || m < 0 || m >= M || S < 1 || NT < 1 || NT > 1024) die("bad args (n >= 4)");
    setup();
    char dir[4096]; snprintf(dir, sizeof dir, "%s/part%d", out, m);
    mkdir(out, 0755); if (mkdir(dir, 0755) && errno != EEXIST) die("cannot create part directory");
    { char dn[4200]; snprintf(dn, sizeof dn, "%s/DONE", dir); unlink(dn); }   /* a rerun must not leave a stale DONE */
    bfiles = xcalloc((size_t)P_, sizeof(BFile));
    for (int b=0;b<P_;b++){ char fn[4200]; snprintf(fn, sizeof fn, "%s/b%d.bin", dir, b); bfiles[b].f = fopen(fn, "wb"); if (!bfiles[b].f) die("cannot open bucket file (check ulimit -n)"); fcntl(fileno(bfiles[b].f), F_SETFD, FD_CLOEXEC); pthread_mutex_init(&bfiles[b].mu, NULL); }
    g_part = m; g_nparts = M; g_slice_count = S; g_slice_mod = M*S;
    pthread_t th[1024];
    for (int i=0;i<NT;i++) if (pthread_create(&th[i], NULL, part_worker, NULL)) die("pthread_create failed");
    for (int i=0;i<NT;i++) pthread_join(th[i], NULL);
    for (int b=0;b<P_;b++) if (fclose(bfiles[b].f)) die("close failed");
    char fn[4200]; snprintf(fn, sizeof fn, "%s/DONE", dir);
    FILE *f = fopen(fn, "w"); if (!f) die("cannot write DONE");
    fprintf(f, "n %d P %d M %d m %d slices %d generated %llu\n", n, P_, M, m, S, (unsigned long long)g_generated);
    if (fclose(f)) die("cannot write DONE");
    fprintf(stderr, "partition n=%d part %d/%d: generated %llu DAGs\n", n, m, M, (unsigned long long)g_generated);
    return 0;
}

// =====================================================================================
// Phase 2: reduce
// =====================================================================================
typedef struct { uint8_t *rec; size_t cap, cnt; pthread_mutex_t mu; } Shard;
typedef struct { Shard sh[NSHARD]; size_t KB_, RB; } HSet;
static HSet *sets[MAXN];
static atomic_int g_autmismatch, g_wrongbucket;

static uint64_t hbytes(const uint8_t *k, size_t len){ uint64_t h = 0x1234567ULL; size_t i=0; for (; i+8<=len; i+=8){ uint64_t w; memcpy(&w, k+i, 8); h = mix64(h ^ w); } uint64_t w=0; memcpy(&w, k+i, len-i); return mix64(h ^ w ^ len); }
static void shard_put(Shard *s, HSet *hs, const uint8_t *key, uint32_t aut, uint64_t h);
static void shard_grow(Shard *s, HSet *hs){
    size_t oc = s->cap; uint8_t *old = s->rec;
    s->cap = oc ? oc*2 : 64; s->cnt = 0; s->rec = xcalloc(s->cap, hs->RB);
    for (size_t i=0;i<oc;i++){ uint8_t *r = old + i*hs->RB; uint32_t a; memcpy(&a, r+hs->KB_, 4); if (a) shard_put(s, hs, r, a, hbytes(r, hs->KB_)); }
    free(old);
}
static void shard_put(Shard *s, HSet *hs, const uint8_t *key, uint32_t aut, uint64_t h){
    if ((s->cnt+1)*10 >= s->cap*7) shard_grow(s, hs);
    size_t m = s->cap-1, i = (h>>12) & m;
    for (;;) {
        uint8_t *r = s->rec + i*hs->RB; uint32_t a; memcpy(&a, r+hs->KB_, 4);
        if (!a) { memcpy(r, key, hs->KB_); memcpy(r+hs->KB_, &aut, 4); s->cnt++; return; }
        if (!memcmp(r, key, hs->KB_)) { if (a != aut) atomic_fetch_add(&g_autmismatch, 1); return; }
        i = (i+1) & m;
    }
}
static void hset_put(HSet *hs, const uint8_t *key, uint32_t aut){
    uint64_t h = hbytes(key, hs->KB_);
    Shard *s = &hs->sh[h & (NSHARD-1)];
    pthread_mutex_lock(&s->mu); shard_put(s, hs, key, aut, h); pthread_mutex_unlock(&s->mu);
}

typedef struct {
    const uint64_t *sep; int L;
    int order[MAXN], cellEnd[MAXN], ncell, pos[MAXN], freecell[MAXN];
    uint64_t best[MAXS], cur[MAXS]; uint64_t nbest; int have;
} Canon;
// is the transposition (x y) an automorphism of the signature sep[0..L)?
static int swap_is_aut(const uint64_t *sep, int L, int x, int y){
    int sw[MAXN]; for (int i=0;i<n;i++) sw[i]=i; sw[x]=y; sw[y]=x;
    for (int i=0;i<L;i++){
        unsigned S = Slist[i], S2 = 0; unsigned T=S; while (T){ int b=__builtin_ctz(T); T&=T-1; S2 |= 1u<<sw[b]; }
        uint64_t m = sep[i], m2 = 0; while (m){ int p=__builtin_ctzll(m); m&=m-1; m2 |= 1ULL<<pidx[sw[pu[p]]][sw[pv[p]]]; }
        if (sep[Sidx[S2]] != m2) return 0;
    }
    return 1;
}
static void eval_perm(Canon *c){
    for (int i=0;i<c->L;i++){
        unsigned S = Slist[i], S2 = 0; while (S){ int b=__builtin_ctz(S); S&=S-1; S2 |= 1u<<c->pos[b]; }
        uint64_t m = c->sep[i], m2 = 0; while (m){ int p=__builtin_ctzll(m); m&=m-1; m2 |= 1ULL<<pidx[c->pos[pu[p]]][c->pos[pv[p]]]; }
        c->cur[Sidx[S2]] = m2;
    }
    int cmp = 0;
    if (!c->have) cmp = -1;
    else for (int i=0;i<c->L;i++){ if (c->cur[i] != c->best[i]) { cmp = c->cur[i] < c->best[i] ? -1 : 1; break; } }
    if (cmp < 0) { memcpy(c->best, c->cur, 8*(size_t)c->L); c->nbest = 1; c->have = 1; }
    else if (cmp == 0) c->nbest++;
}
static void perm_cells(Canon *c, int ci, int start);
static void perm_within(Canon *c, int ci, int i, int hi){
    if (i == hi) { perm_cells(c, ci+1, hi); return; }
    for (int j=i;j<hi;j++){
        int t=c->order[i]; c->order[i]=c->order[j]; c->order[j]=t;
        c->pos[c->order[i]] = i;
        perm_within(c, ci, i+1, hi);
        t=c->order[i]; c->order[i]=c->order[j]; c->order[j]=t;
    }
}
static void perm_cells(Canon *c, int ci, int start){
    if (ci == c->ncell) { eval_perm(c); return; }
    if (c->freecell[ci]) { perm_cells(c, ci+1, c->cellEnd[ci]); return; }   /* one ordering suffices */
    perm_within(c, ci, start, c->cellEnd[ci]);
}
static void canonicalize(Canon *c, const uint64_t *sep, int k){
    int L = sizeEnd[k]; c->sep = sep; c->L = L;
    int allzero = 1; for (int i=0;i<L;i++) if (sep[i]) { allzero = 0; break; }
    if (allzero) { memset(c->best, 0, 8*(size_t)L); c->nbest = nfact; c->have = 1; return; }
    uint64_t inv[MAXN]; vertex_invariants(sep, k, inv);
    for (int i=0;i<n;i++) c->order[i]=i;
    for (int i=1;i<n;i++){ int v=c->order[i], j=i; while (j>0 && inv[c->order[j-1]] > inv[v]) { c->order[j]=c->order[j-1]; j--; } c->order[j]=v; }
    c->ncell=0; for (int i=0;i<n;i++){ if (i==n-1 || inv[c->order[i]] != inv[c->order[i+1]]) c->cellEnd[c->ncell++] = i+1; }
    // Free cells: if every adjacent transposition of a cell's vertices is an automorphism of the
    // signature, the whole symmetric group on that cell is. Then all orderings of the cell give the
    // same permuted signature, so one ordering is evaluated and the count is multiplied by |cell|!.
    uint64_t mult = 1; int start = 0;
    for (int ci=0; ci<c->ncell; ci++){
        int end = c->cellEnd[ci], sz = end - start;
        c->freecell[ci] = 0;
        if (sz >= 2) {
            int ok = 1;
            for (int i=start; i+1<end && ok; i++) ok = swap_is_aut(sep, L, c->order[i], c->order[i+1]);
            if (ok) { c->freecell[ci] = 1; for (int f=2; f<=sz; f++) mult *= (uint64_t)f; }
        }
        start = end;
    }
    for (int i=0;i<n;i++) c->pos[c->order[i]] = i;
    c->have = 0; c->nbest = 0;
    perm_cells(c, 0, 0);
    c->nbest *= mult;
}
static void pack(const uint64_t *sig, int L, uint8_t *out){ for (int i=0;i<L;i++) for (int b=0;b<MB;b++) out[i*MB+b] = (uint8_t)(sig[i] >> (8*b)); }

static uint16_t (*recs)[MAXN]; static size_t nrecs;
static int R_P, R_b, R_msub, R_j;
static atomic_size_t g_next; static atomic_uint_fast64_t g_kept;
static void *reduce_worker(void *arg){
    (void)arg;
    Canon *c = xcalloc(1, sizeof(Canon)); uint64_t sep[MAXS]; uint8_t key[MAXS*8]; uint64_t kept = 0;
    for (;;){
        size_t start = atomic_fetch_add(&g_next, 1024); if (start >= nrecs) break;
        size_t end = start + 1024 < nrecs ? start + 1024 : nrecs;
        for (size_t r=start; r<end; r++){
            Dag d; for (int x=0;x<n;x++) d.par[x] = recs[r][x];
            dag_prepare(&d); compute_sep(&d, sep, 0, sizeEnd[KB]);
            uint64_t h = class_hash(sep);
            if (h % (uint64_t)R_P != (uint64_t)R_b) { atomic_fetch_add(&g_wrongbucket, 1); continue; }
            compute_sep(&d, sep, sizeEnd[KB], sizeEnd[KMAX]);
            if (R_msub > 1 && class_hash_k(sep, KMIN) % (uint64_t)R_msub != (uint64_t)R_j) continue;
            kept++;
            for (int k=KMIN;k<=KMAX;k++){ canonicalize(c, sep, k); pack(c->best, sizeEnd[k], key); hset_put(sets[k], key, (uint32_t)c->nbest); }
        }
    }
    atomic_fetch_add(&g_kept, kept); free(c); return NULL;
}
static int do_reduce(int argc, char **argv){
    if (argc < 11) die("usage: reduce n kmin kmax P b msub j threads indir");
    n = atoi(argv[2]); KMIN = atoi(argv[3]); KMAX = atoi(argv[4]); R_P = atoi(argv[5]); R_b = atoi(argv[6]); R_msub = atoi(argv[7]); R_j = atoi(argv[8]); NT = atoi(argv[9]); const char *in = argv[10];
    if (n < 4 || n > MAXN || KMIN < KB || KMAX > n-2 || KMIN > KMAX || R_P < 1 || R_b < 0 || R_b >= R_P || R_msub < 1 || R_j < 0 || R_j >= R_msub || NT < 1 || NT > 1024) die("bad args (n >= 4, kmin >= 2)");
    setup();
    {   // the partition must be complete: M DONE files, all agreeing on n and P, generated total = A003087(n)
        char dp[4200]; snprintf(dp, sizeof dp, "%s/part*/DONE", in);
        glob_t gd; if (glob(dp, 0, NULL, &gd) != 0) die("no partition DONE files");
        unsigned long long tot = 0; int Mfirst = -1;
        for (size_t i=0;i<gd.gl_pathc;i++){
            FILE *f = fopen(gd.gl_pathv[i], "r"); if (!f) die("cannot read DONE");
            int dn, dP, dM, dm, dS; unsigned long long dg;
            if (fscanf(f, "n %d P %d M %d m %d slices %d generated %llu", &dn, &dP, &dM, &dm, &dS, &dg) != 6) die("bad DONE file");
            fclose(f);
            if (dn != n || dP != R_P) die("DONE file does not match n or P");
            if (Mfirst < 0) Mfirst = dM; else if (dM != Mfirst) die("DONE files disagree on M");
            tot += dg;
        }
        if ((int)gd.gl_pathc != Mfirst) { fprintf(stderr, "ERROR: %zu of %d partition parts done\n", gd.gl_pathc, Mfirst); exit(2); }
        if (tot != A003087[n]) { fprintf(stderr, "ERROR: partition total %llu != A003087(%d)\n", tot, n); exit(2); }
        globfree(&gd);
    }
    char pat[4200]; snprintf(pat, sizeof pat, "%s/part*/b%d.bin", in, R_b);
    glob_t g; if (glob(pat, 0, NULL, &g) != 0) die("no bucket files found");
    size_t total = 0;
    for (size_t i=0;i<g.gl_pathc;i++){ struct stat st; if (stat(g.gl_pathv[i], &st)) die("stat failed"); if (st.st_size % (MAXN*2)) die("bucket file size not a multiple of the record size"); total += (size_t)st.st_size / (MAXN*2); }
    recs = malloc(total ? total * sizeof *recs : 1); if (!recs) die("out of memory reading records");
    nrecs = 0;
    for (size_t i=0;i<g.gl_pathc;i++){
        struct stat st; if (stat(g.gl_pathv[i], &st)) die("stat failed"); size_t want_r = (size_t)st.st_size / (MAXN*2);
        FILE *f = fopen(g.gl_pathv[i], "rb"); if (!f) die("open failed");
        size_t got = fread(recs[nrecs], sizeof *recs, want_r, f);
        if (got != want_r || ferror(f)) die("short read of a bucket file");
        nrecs += got; fclose(f);
    }
    if (nrecs != total) die("short read");
    fprintf(stderr, "reduce n=%d k=%d..%d bucket %d/%d sub %d/%d: %zu part files, %zu DAGs\n", n, KMIN, KMAX, R_b, R_P, R_j, R_msub, g.gl_pathc, nrecs);
    for (int k=KMIN;k<=KMAX;k++){ sets[k] = xcalloc(1, sizeof(HSet)); sets[k]->KB_ = (size_t)sizeEnd[k]*MB; sets[k]->RB = sets[k]->KB_ + 4; for (int i=0;i<NSHARD;i++) pthread_mutex_init(&sets[k]->sh[i].mu, NULL); }
    pthread_t th[1024];
    for (int i=0;i<NT;i++) if (pthread_create(&th[i], NULL, reduce_worker, NULL)) die("pthread_create failed");
    for (int i=0;i<NT;i++) pthread_join(th[i], NULL);
    int bad = 0;
    if (g_wrongbucket) { fprintf(stderr, "ERROR: %d DAGs in the wrong bucket file\n", (int)g_wrongbucket); bad = 1; }
    if (g_autmismatch) { fprintf(stderr, "ERROR: %d automorphism-count mismatches\n", (int)g_autmismatch); bad = 1; }
    printf("# n %d kmin %d kmax %d P %d b %d msub %d j %d read %zu kept %llu\n", n, KMIN, KMAX, R_P, R_b, R_msub, R_j, nrecs, (unsigned long long)g_kept);
    for (int k=KMIN;k<=KMAX;k++){
        size_t U = 0; unsigned __int128 T = 0; HSet *hs = sets[k];
        for (int i=0;i<NSHARD;i++){ Shard *s=&hs->sh[i]; U += s->cnt;
            for (size_t q=0;q<s->cap;q++){ uint32_t a; memcpy(&a, s->rec + q*hs->RB + hs->KB_, 4); if (a) { if (nfact % a) { fprintf(stderr,"ERROR: |Aut| %u does not divide n!\n", a); bad = 1; } T += nfact / a; } } }
        printf("%d %d %zu %llu\n", n, k, U, (unsigned long long)T);
    }
    fflush(stdout);
    return bad ? 3 : 0;
}

// =====================================================================================
// Direct mode: generate all DAGs and count k = kmin..kmax (any k >= 0) without buckets.
//   ./kmec_cluster direct n kmin kmax threads slices
// For small k, where the number of classes is small. Threads pull directg -s slices from a counter.
// =====================================================================================
static atomic_int g_dslice; static int g_dslices; static atomic_uint_fast64_t g_dgen;
static void *direct_worker(void *arg){
    (void)arg;
    Canon *c = xcalloc(1, sizeof(Canon)); uint64_t sep[MAXS]; uint8_t key[MAXS*8]; uint64_t gen = 0;
    for (;;){
        int s = atomic_fetch_add(&g_dslice, 1); if (s >= g_dslices) break;
        char cmd[256]; slice_cmd(cmd, sizeof cmd, s, g_dslices);
        FILE *in = popen(cmd, "r"); if (!in) die("popen failed");
        Dag d;
        while (read_dag(in, &d)){
            dag_prepare(&d); compute_sep(&d, sep, 0, sizeEnd[KMAX]);
            for (int k=KMIN;k<=KMAX;k++){ canonicalize(c, sep, k); pack(c->best, sizeEnd[k], key); hset_put(sets[k], key, (uint32_t)c->nbest); }
            gen++;
        }
        int st = pclose(in); if (st != 0) { fprintf(stderr, "ERROR: nauty slice %d/%d exited with status %d\n", s, g_dslices, st); exit(2); }
    }
    atomic_fetch_add(&g_dgen, gen); free(c); return NULL;
}
static int print_sets(const char *hdr){
    int bad = 0;
    if (g_autmismatch) { fprintf(stderr, "ERROR: %d automorphism-count mismatches\n", (int)g_autmismatch); bad = 1; }
    printf("%s\n", hdr);
    for (int k=KMIN;k<=KMAX;k++){
        size_t U = 0; unsigned __int128 T = 0; HSet *hs = sets[k];
        for (int i=0;i<NSHARD;i++){ Shard *sh=&hs->sh[i]; U += sh->cnt;
            for (size_t q=0;q<sh->cap;q++){ uint32_t a; memcpy(&a, sh->rec + q*hs->RB + hs->KB_, 4); if (a) { if (nfact % a) { fprintf(stderr,"ERROR: |Aut| %u does not divide n!\n", a); bad = 1; } T += nfact / a; } } }
        printf("%d %d %zu %llu\n", n, k, U, (unsigned long long)T);
    }
    fflush(stdout);
    return bad;
}
static int do_direct(int argc, char **argv){
    if (argc < 7) die("usage: direct n kmin kmax threads slices");
    n = atoi(argv[2]); KMIN = atoi(argv[3]); KMAX = atoi(argv[4]); NT = atoi(argv[5]); g_dslices = atoi(argv[6]);
    if (n < 2 || n > MAXN || KMIN < 0 || KMAX > n-2 || KMIN > KMAX || NT < 1 || NT > 1024 || g_dslices < 1) die("bad args");
    setup();
    for (int k=KMIN;k<=KMAX;k++){ sets[k] = xcalloc(1, sizeof(HSet)); sets[k]->KB_ = (size_t)sizeEnd[k]*MB; sets[k]->RB = sets[k]->KB_ + 4; for (int i=0;i<NSHARD;i++) pthread_mutex_init(&sets[k]->sh[i].mu, NULL); }
    pthread_t th[1024];
    for (int i=0;i<NT;i++) if (pthread_create(&th[i], NULL, direct_worker, NULL)) die("pthread_create failed");
    for (int i=0;i<NT;i++) pthread_join(th[i], NULL);
    char hdr[256]; snprintf(hdr, sizeof hdr, "# n %d kmin %d kmax %d direct slices %d generated %llu", n, KMIN, KMAX, g_dslices, (unsigned long long)g_dgen);
    int bad = 0;
    if (g_dgen != A003087[n]) { fprintf(stderr, "ERROR: generated %llu DAGs, expected A003087(%d) = %llu\n", (unsigned long long)g_dgen, n, (unsigned long long)A003087[n]); bad = 1; }
    bad |= print_sets(hdr);
    return bad ? 3 : 0;
}

int main(int argc, char **argv){
    if (argc >= 2 && !strcmp(argv[1], "partition")) return do_partition(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "reduce")) return do_reduce(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "direct")) return do_direct(argc, argv);
    fprintf(stderr, "usage:\n  %s partition n P M m slices_per_part threads outdir\n  %s reduce n kmin kmax P b msub j threads indir\n  %s direct n kmin kmax threads slices\n", argv[0], argv[0], argv[0]);
    return 1;
}
