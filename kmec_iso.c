// kmec_iso.c — exact counts of k-Markov equivalence classes of DAGs via isomorphism classes.
//
// Input: one unlabeled DAG per line in nauty's directg -T format ("nv ne v1 w1 v2 w2 ..."),
// for example:  geng 7 | directg -a -T | ./kmec_iso 7 0 5 14
// Every isomorphism class of DAGs must appear exactly once (directg guarantees this).
//
// For each DAG and each k, the k-signature (the d-separation pattern over all S with |S| <= k)
// is put in canonical form: the lexicographic minimum over all node relabelings that respect
// isomorphism-invariant vertex cells, together with the number of relabelings attaining the
// minimum. That number is the size of the signature's automorphism group, Aut.
// Canonical signatures are stored exactly (no hashing).
//   unlabeled count U(n,k) = number of distinct canonical signatures
//   labeled count   T(n,k) = sum over distinct canonical signatures of n! / |Aut|
// The labeled count is exact, so it independently checks the hashed counts from kmec2.c.
//
// Usage: ./kmec_iso n kmin kmax threads < dags.txt
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>

#define MAXN 8
#define MAXS 256
#define NSHARD 4096

static int n, KMIN, KMAX, NT;
static int npairs, pidx[MAXN][MAXN], pu[32], pv[32];
static int nS; static unsigned Slist[MAXS]; static int Sidx[1 << MAXN]; static int sizeEnd[MAXN];
static uint64_t nfact = 1;

static inline uint64_t mix64(uint64_t z){ z += 0x9e3779b97f4a7c15ULL; z = (z ^ (z>>30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z>>27)) * 0x94d049bb133111ebULL; return z ^ (z>>31); }
static inline unsigned orset(const unsigned *a, unsigned X){ unsigned r=0; while(X){ int i=__builtin_ctz(X); X&=X-1; r|=a[i]; } return r; }

// ---------- exact sharded hash set: records of (L key words + 1 aut word) ----------
typedef struct { uint32_t *rec; size_t cap, cnt; pthread_mutex_t mu; } Shard;
typedef struct { Shard sh[NSHARD]; int L; } HSet;
static HSet *sets[MAXN];
static atomic_int g_autmismatch;

static uint64_t hkey(const uint32_t *k, int L){ uint64_t h = 0x1234567ULL; for (int i=0;i<L;i++) h = mix64(h ^ k[i] ^ ((uint64_t)i<<40)); return h; }
static void shard_put(Shard *s, int L, const uint32_t *key, uint32_t aut, uint64_t h);
static void shard_grow(Shard *s, int L){
    size_t oc = s->cap; uint32_t *old = s->rec;
    s->cap = oc ? oc*2 : 64; s->cnt = 0;
    s->rec = calloc(s->cap * (size_t)(L+1), 4);
    if (!s->rec) { fprintf(stderr, "out of memory\n"); exit(2); }
    for (size_t i=0;i<oc;i++){ uint32_t *r = old + i*(size_t)(L+1); if (r[L]) shard_put(s, L, r, r[L], hkey(r, L)); }
    free(old);
}
static void shard_put(Shard *s, int L, const uint32_t *key, uint32_t aut, uint64_t h){
    if ((s->cnt+1)*10 >= s->cap*7) shard_grow(s, L);
    size_t m = s->cap-1, i = (h>>12) & m;
    for (;;) {
        uint32_t *r = s->rec + i*(size_t)(L+1);
        if (!r[L]) { memcpy(r, key, 4*(size_t)L); r[L] = aut; s->cnt++; return; }
        if (!memcmp(r, key, 4*(size_t)L)) { if (r[L] != aut) atomic_fetch_add(&g_autmismatch, 1); return; }
        i = (i+1) & m;
    }
}
static void hset_put(HSet *hs, const uint32_t *key, uint32_t aut){
    uint64_t h = hkey(key, hs->L);
    Shard *s = &hs->sh[h & (NSHARD-1)];
    pthread_mutex_lock(&s->mu); shard_put(s, hs->L, key, aut, h); pthread_mutex_unlock(&s->mu);
}

// ---------- d-separation (Bayes-ball, same rules as kmec2.c) ----------
static inline unsigned bayesball(const unsigned *par, const unsigned *ch, int u, unsigned S, unsigned A){
    unsigned visU = 1u<<u, visD = 0, newU = visU, newD = 0;
    while (newU | newD) {
        unsigned upOpen = newU & ~S, dnChain = newD & ~S, dnColl = newD & A;
        unsigned nu = orset(par, upOpen | dnColl), nd = orset(ch, upOpen | dnChain);
        newU = nu & ~visU; newD = nd & ~visD; visU |= newU; visD |= newD;
    }
    return (visU | visD) & ~S;
}

// ---------- per-thread canonicalization state ----------
typedef struct {
    const uint32_t *sep; int L;
    int order[MAXN], cellEnd[MAXN], ncell;   // vertices sorted by invariant; cell boundaries
    int pos[MAXN];                           // pos[x] = new label of old vertex x
    uint32_t best[MAXS], cur[MAXS]; uint64_t nbest; int have;
} Canon;

static void eval_perm(Canon *c){
    for (int i=0;i<c->L;i++){
        unsigned S = Slist[i], S2 = 0; while (S){ int b=__builtin_ctz(S); S&=S-1; S2 |= 1u<<c->pos[b]; }
        uint32_t m = c->sep[i], m2 = 0; while (m){ int p=__builtin_ctz(m); m&=m-1; m2 |= 1u<<pidx[c->pos[pu[p]]][c->pos[pv[p]]]; }
        c->cur[Sidx[S2]] = m2;
    }
    int cmp = c->have ? memcmp(c->cur, c->best, 4*(size_t)c->L) : -1;
    if (cmp < 0) { memcpy(c->best, c->cur, 4*(size_t)c->L); c->nbest = 1; c->have = 1; }
    else if (cmp == 0) c->nbest++;
}
// permute the vertices of cell ci among that cell's positions (recursive over cells, then within a cell)
static void perm_cells(Canon *c, int ci, int start);
static void perm_within(Canon *c, int ci, int lo, int i, int hi){
    if (i == hi) { perm_cells(c, ci+1, hi); return; }
    for (int j=i;j<hi;j++){
        int t=c->order[i]; c->order[i]=c->order[j]; c->order[j]=t;
        c->pos[c->order[i]] = i;
        perm_within(c, ci, lo, i+1, hi);
        t=c->order[i]; c->order[i]=c->order[j]; c->order[j]=t;
    }
}
static void perm_cells(Canon *c, int ci, int start){
    if (ci == c->ncell) { eval_perm(c); return; }
    perm_within(c, ci, start, start, c->cellEnd[ci]);
}

static int cmp_u64pair(const void *a, const void *b){ const uint64_t *x=a,*y=b; return (x[0]>y[0])-(x[0]<y[0]); }

static void canonicalize(Canon *c, const uint32_t *sep, int k){
    int L = sizeEnd[k]; c->sep = sep; c->L = L;
    // invariant round 1: for each vertex x and each size s, separated statements with x as an endpoint, and with x in S
    uint64_t inv[MAXN], inv2[MAXN];
    uint32_t pairc[MAXN][MAXN][MAXN]; memset(pairc, 0, sizeof pairc); // pairc[x][y][s]: S of size s with x _|_ y | S
    uint32_t inS[MAXN][MAXN]; memset(inS, 0, sizeof inS);
    for (int i=0;i<L;i++){
        unsigned S = Slist[i]; int s = __builtin_popcount(S);
        uint32_t m = sep[i]; int cntm = __builtin_popcount(m);
        while (m){ int p=__builtin_ctz(m); m&=m-1; pairc[pu[p]][pv[p]][s]++; pairc[pv[p]][pu[p]][s]++; }
        unsigned T=S; while (T){ int b=__builtin_ctz(T); T&=T-1; inS[b][s] += cntm; }
    }
    for (int x=0;x<n;x++){
        uint64_t h = 0xabcdefULL;
        for (int s=0;s<=k;s++){ uint32_t c1=0; for (int y=0;y<n;y++) if (y!=x) c1 += pairc[x][y][s]; h = mix64(h ^ c1 ^ ((uint64_t)inS[x][s]<<32) ^ ((uint64_t)s<<56)); }
        inv[x] = h;
    }
    // refinement rounds: combine each vertex's invariant with the sorted multiset of (neighbor invariant, pair pattern)
    for (int round=0; round<2; round++){
        for (int x=0;x<n;x++){
            uint64_t lst[MAXN][1]; int L2=0;
            for (int y=0;y<n;y++){ if (y==x) continue; uint64_t h = mix64(inv[y]); for (int s=0;s<=k;s++) h = mix64(h ^ pairc[x][y][s] ^ ((uint64_t)s<<48)); lst[L2++][0]=h; }
            qsort(lst, L2, sizeof lst[0], cmp_u64pair);
            uint64_t h = mix64(inv[x] ^ 0x5555ULL); for (int i=0;i<L2;i++) h = mix64(h ^ lst[i][0]);
            inv2[x] = h;
        }
        memcpy(inv, inv2, sizeof inv);
    }
    // order vertices by invariant (ties = one cell)
    for (int i=0;i<n;i++) c->order[i]=i;
    for (int i=1;i<n;i++){ int v=c->order[i], j=i; while (j>0 && inv[c->order[j-1]] > inv[v]) { c->order[j]=c->order[j-1]; j--; } c->order[j]=v; }
    c->ncell=0; for (int i=0;i<n;i++){ if (i==n-1 || inv[c->order[i]] != inv[c->order[i+1]]) c->cellEnd[c->ncell++] = i+1; }
    for (int i=0;i<n;i++) c->pos[c->order[i]] = i;
    c->have = 0; c->nbest = 0;
    perm_cells(c, 0, 0);
}

// ---------- input and work distribution ----------
static unsigned (*dags)[MAXN]; static size_t ndags;
static atomic_size_t g_next;
static atomic_uint_fast64_t g_done;

static void *worker(void *arg){
    (void)arg;
    Canon *c = malloc(sizeof(Canon));
    uint32_t sep[MAXS];
    for (;;) {
        size_t start = atomic_fetch_add(&g_next, 4096);
        if (start >= ndags) break;
        size_t end = start + 4096 < ndags ? start + 4096 : ndags;
        for (size_t d=start; d<end; d++){
            unsigned par[MAXN], ch[MAXN]={0}, anc[MAXN];
            memcpy(par, dags[d], sizeof par);
            for (int x=0;x<n;x++){ unsigned p=par[x]; while(p){ int i=__builtin_ctz(p); p&=p-1; ch[i]|=1u<<x; } }
            for (int x=0;x<n;x++){ unsigned a=1u<<x, prev=0; while (a!=prev){ prev=a; a |= orset(par, a); } anc[x]=a; }
            for (int i=0;i<sizeEnd[KMAX];i++){
                unsigned S = Slist[i], A = orset(anc, S) | S; uint32_t m = 0;
                for (int u=0;u<n;u++){ if (S>>u&1) continue; unsigned R = bayesball(par, ch, u, S, A);
                    for (int v=u+1; v<n; v++) if (!(S>>v&1) && !(R>>v&1)) m |= 1u<<pidx[u][v]; }
                sep[i] = m;
            }
            for (int k=KMIN;k<=KMAX;k++){ canonicalize(c, sep, k); hset_put(sets[k], c->best, (uint32_t)c->nbest); }
            uint64_t done = atomic_fetch_add(&g_done, 1) + 1;
            if ((done & ((1u<<21)-1)) == 0) fprintf(stderr, "progress: %llu / %zu DAGs\n", (unsigned long long)done, ndags);
        }
    }
    free(c); return NULL;
}

int main(int argc, char **argv){
    if (argc < 5) { fprintf(stderr, "usage: %s n kmin kmax threads < directg_T_output\n", argv[0]); return 1; }
    n = atoi(argv[1]); KMIN = atoi(argv[2]); KMAX = atoi(argv[3]); NT = atoi(argv[4]);
    if (n < 2 || n > MAXN || KMIN < 0 || KMAX > n-2 || KMIN > KMAX || NT < 1 || NT > 256) { fprintf(stderr, "bad args\n"); return 1; }
    for (int i=2;i<=n;i++) nfact *= i;
    npairs=0; for (int u=0;u<n;u++) for (int v=u+1;v<n;v++){ pidx[u][v]=pidx[v][u]=npairs; pu[npairs]=u; pv[npairs]=v; npairs++; }
    nS=0; for (int s=0;s<=n-2;s++){ for (unsigned S=0;S<(1u<<n);S++) if (__builtin_popcount(S)==s){ Sidx[S]=nS; Slist[nS++]=S; } sizeEnd[s]=nS; }
    // read DAGs: "nv ne v1 w1 ..." (arc v -> w)
    size_t cap = 1<<20; dags = malloc(cap * sizeof *dags); ndags = 0;
    int nv, ne;
    while (scanf("%d %d", &nv, &ne) == 2){
        if (nv != n) { fprintf(stderr, "input graph has %d nodes, expected %d\n", nv, n); return 1; }
        if (ndags == cap) { cap *= 2; dags = realloc(dags, cap * sizeof *dags); }
        memset(dags[ndags], 0, sizeof dags[ndags]);
        for (int e=0;e<ne;e++){ int a,b; if (scanf("%d %d",&a,&b)!=2) { fprintf(stderr,"parse error\n"); return 1; } dags[ndags][b] |= 1u<<a; }
        ndags++;
    }
    fprintf(stderr, "read %zu unlabeled DAGs on %d nodes\n", ndags, n);
    for (int k=KMIN;k<=KMAX;k++){ sets[k] = calloc(1, sizeof(HSet)); sets[k]->L = sizeEnd[k]; for (int i=0;i<NSHARD;i++) pthread_mutex_init(&sets[k]->sh[i].mu, NULL); }
    pthread_t th[256];
    for (int i=0;i<NT;i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i=0;i<NT;i++) pthread_join(th[i], NULL);
    if (g_autmismatch) fprintf(stderr, "WARNING: %d automorphism-count mismatches\n", (int)g_autmismatch);
    printf("# n k unlabeled labeled   (unlabeled DAGs read: %zu)\n", ndags);
    for (int k=KMIN;k<=KMAX;k++){
        size_t U = 0; unsigned __int128 T = 0; int L = sets[k]->L;
        for (int i=0;i<NSHARD;i++){ Shard *s=&sets[k]->sh[i]; U += s->cnt;
            for (size_t j=0;j<s->cap;j++){ uint32_t a = s->rec[j*(size_t)(L+1)+L]; if (a) { if (nfact % a) fprintf(stderr,"WARNING: |Aut| %u does not divide n!\n", a); T += nfact / a; } } }
        printf("%d %d %zu %llu\n", n, k, U, (unsigned long long)T); fflush(stdout);
    }
    return g_autmismatch ? 3 : 0;
}
