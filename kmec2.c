// kmec2.c — fast, multithreaded count of k-Markov equivalence classes of
// labeled DAGs on n nodes (n <= 8 in principle; n = 7 is the target).
//
// Definition: DAGs G, H are k-Markov equivalent iff for all distinct u,v and
// all S ⊆ V\{u,v} with |S| <= k:  u ⊥_G v | S  <=>  u ⊥_H v | S.
//
// Method (brute force by definition, no characterization theorem used):
//  * Enumerate every labeled DAG exactly once, node by node: node j chooses a
//    parent set P and child set C among nodes < j; valid iff no child reaches a
//    parent (descendant masks are maintained incrementally).
//  * For each DAG, for every S (in a fixed order, by increasing |S|) and every
//    u ∉ S, compute the set of nodes d-connected to u given S with a bitmask
//    Bayes-ball (reachability over (node, direction) states).
//  * The d-separation pattern for each S (a bitmask over pairs) is streamed into
//    a 128-bit fingerprint; the running fingerprint after all S of size <= k is
//    the k-signature. Distinct fingerprints are counted per k in a sharded,
//    mutex-striped hash set.
//  * Fingerprints are hashes: two DAGs with different signatures collide with
//    probability ~2^-128 per pair; for N <= 1e9 classes the total collision
//    probability is < 1e-20. (kmec.c is the exact-signature version, n <= 6.)
//
// Usage: ./kmec2 n [kmax] [threads] [seed]
// Different seeds give independent fingerprint functions; identical counts across
// seeds are the empirical check against hash collisions.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>

#define MAXN 8
#define NSHARD 4096

static int n, K, NT;
static uint64_t SEED = 0;
static int npairs, pidx[MAXN][MAXN];
static int nS; static unsigned Slist[256];
static int sizeEnd[MAXN]; // Slist index one past the last S of size s

// ---------------- fingerprint ----------------
static inline uint64_t mix64(uint64_t z){ z += 0x9e3779b97f4a7c15ULL; z = (z ^ (z>>30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z>>27)) * 0x94d049bb133111ebULL; return z ^ (z>>31); }

// ---------------- sharded hash set of 128-bit keys ----------------
typedef struct { uint64_t *k; size_t cap, cnt; pthread_mutex_t mu; } Shard;
typedef struct { Shard sh[NSHARD]; } HSet;
static HSet *sets[MAXN];

static void shard_put_nolock(Shard *s, uint64_t a, uint64_t b);
static void shard_grow(Shard *s){
    size_t oc = s->cap; uint64_t *ok = s->k;
    s->cap = oc ? oc*2 : 1024; s->k = calloc(s->cap*2, 8); s->cnt = 0;
    if (!s->k) { fprintf(stderr,"OOM\n"); exit(2); }
    for (size_t i=0;i<oc;i++) if (ok[2*i]|ok[2*i+1]) shard_put_nolock(s, ok[2*i], ok[2*i+1]);
    free(ok);
}
static void shard_put_nolock(Shard *s, uint64_t a, uint64_t b){
    if ((s->cnt+1)*10 >= s->cap*7) shard_grow(s);
    size_t m = s->cap-1, i = (b>>12) & m;
    for (;;) {
        uint64_t *e = s->k + 2*i;
        if (!(e[0]|e[1])) { e[0]=a; e[1]=b; s->cnt++; return; }
        if (e[0]==a && e[1]==b) return;
        i = (i+1) & m;
    }
}
static void hset_put(HSet *h, uint64_t a, uint64_t b){
    if (!(a|b)) a = 1; // reserve all-zero as empty
    Shard *s = &h->sh[b & (NSHARD-1)];
    pthread_mutex_lock(&s->mu); shard_put_nolock(s, a, b); pthread_mutex_unlock(&s->mu);
}

// ---------------- per-DAG work ----------------
static inline unsigned orpar(const unsigned *par, unsigned X){ unsigned r=0; while(X){ int i=__builtin_ctz(X); X&=X-1; r|=par[i]; } return r; }

// nodes d-connected to u given S (excluding S; includes u itself)
static inline unsigned bayesball(const unsigned *par, const unsigned *ch, int u, unsigned S, unsigned A){
    unsigned visU = 1u<<u, visD = 0, newU = visU, newD = 0;
    while (newU | newD) {
        unsigned upOpen = newU & ~S;       // arrived from a child, not conditioned: pass to parents and children
        unsigned dnChain = newD & ~S;      // arrived from a parent, not conditioned: pass to children
        unsigned dnColl = newD & A;        // arrived from a parent, collider active: pass to parents
        unsigned nu = orpar(par, upOpen | dnColl);
        unsigned nd = orpar(ch, upOpen | dnChain);
        newU = nu & ~visU; newD = nd & ~visD;
        visU |= newU; visD |= newD;
    }
    return (visU | visD) & ~S;
}

static atomic_uint_fast64_t g_ndag;

static void process_dag(const unsigned *par, const unsigned *reach){
    unsigned ch[MAXN] = {0}, anc[MAXN] = {0};
    for (int x=0;x<n;x++){ unsigned p=par[x]; while(p){ int i=__builtin_ctz(p); p&=p-1; ch[i]|=1u<<x; } }
    for (int x=0;x<n;x++){ unsigned r=reach[x]; while(r){ int i=__builtin_ctz(r); r&=r-1; anc[i]|=1u<<x; } }
    // SEED = 0 (default) reproduces the original build used for results_n7_seedA.txt.
    uint64_t h1 = 0x1234567887654321ULL ^ (SEED ? mix64(SEED) : 0), h2 = 0x0fedcba987654321ULL ^ (SEED ? mix64(~SEED) : 0);
    int si = 0;
    for (int s=0; s<=K; s++) {
        for (; si < sizeEnd[s]; si++) {
            unsigned S = Slist[si];
            unsigned A = orpar(anc, S) | S; // ancestors of S (incl. S)
            uint64_t sepmask = 0;
            for (int u=0; u<n; u++) {
                if (S>>u & 1) continue;
                unsigned R = bayesball(par, ch, u, S, A);
                for (int v=u+1; v<n; v++) if (!(S>>v & 1) && !(R>>v & 1)) sepmask |= 1ULL << pidx[u][v];
            }
            h1 = mix64(h1 ^ sepmask);
            h2 = mix64(h2 + (sepmask * (0xd6e8feb86659fd93ULL ^ (SEED<<1))) + (uint64_t)si);
        }
        hset_put(sets[s], h1, h2);
    }
    uint64_t c = atomic_fetch_add_explicit(&g_ndag, 1, memory_order_relaxed) + 1;
    if ((c & ((1u<<24)-1)) == 0) { fprintf(stderr, "progress: %llu DAGs\n", (unsigned long long)c); }
}

// ---------------- enumeration (node by node) with dynamic task split ----------------
#define SPLIT_DEPTH 5   // prefixes = valid DAGs on the first SPLIT_DEPTH nodes
static atomic_uint_fast64_t g_next_task;

typedef struct { uint64_t counter, mine; } TaskCtx;

static void rec(int j, unsigned *par, unsigned *reach, TaskCtx *t){
    if (j == n) { process_dag(par, reach); return; }
    if (j == SPLIT_DEPTH || (SPLIT_DEPTH >= n && j == n-1)) {
        // task boundary: only descend into prefixes we have claimed
        uint64_t id = t->counter++;
        if (id != t->mine) return;
        t->mine = atomic_fetch_add(&g_next_task, 1);
    }
    unsigned full = (1u<<j) - 1;
    // P = parents of j, C = children of j, disjoint subsets of [0,j)
    for (unsigned P = 0; P <= full; P++) {
        unsigned rest = full & ~P;
        unsigned C = rest;
        for (;;) {
            unsigned D = orpar(reach, C); // descendants of children
            if (!(D & P)) {
                unsigned np[MAXN], nr[MAXN];
                memcpy(np, par, sizeof(unsigned)*n); memcpy(nr, reach, sizeof(unsigned)*n);
                np[j] = P;
                unsigned c = C; while (c) { int i=__builtin_ctz(c); c&=c-1; np[i] |= 1u<<j; }
                nr[j] = (1u<<j) | D;
                for (int x=0;x<j;x++) if (nr[x] & P) nr[x] |= nr[j];
                rec(j+1, np, nr, t);
            }
            if (C == 0) break;
            C = (C-1) & rest;
        }
    }
}

static void *worker(void *arg){
    (void)arg;
    TaskCtx t = { 0, atomic_fetch_add(&g_next_task, 1) };
    unsigned par[MAXN]={0}, reach[MAXN]={0};
    rec(0, par, reach, &t);
    return NULL;
}

int main(int argc, char **argv){
    if (argc < 2) { fprintf(stderr,"usage: %s n [kmax] [threads]\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    K = argc > 2 ? atoi(argv[2]) : n-2;
    NT = argc > 3 ? atoi(argv[3]) : 1;
    SEED = argc > 4 ? strtoull(argv[4], NULL, 0) : 0;
    if (n < 2 || n > MAXN || K < 0 || K > n-2 || NT < 1 || NT > 256) { fprintf(stderr,"bad args\n"); return 1; }
    npairs = 0; for (int u=0;u<n;u++) for (int v=u+1;v<n;v++) pidx[u][v] = npairs++;
    nS = 0;
    for (int s=0;s<=n-2;s++){ for (unsigned S=0;S<(1u<<n);S++) if (__builtin_popcount(S)==s && __builtin_popcount(S) <= n-2) { Slist[nS]=S; nS++; } sizeEnd[s]=nS; }
    for (int s=0;s<=K;s++){ sets[s] = calloc(1,sizeof(HSet)); for (int i=0;i<NSHARD;i++) pthread_mutex_init(&sets[s]->sh[i].mu, NULL); }
    pthread_t th[256];
    for (int i=0;i<NT;i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i=0;i<NT;i++) pthread_join(th[i], NULL);
    fprintf(stderr, "n=%d labeled DAGs=%llu\n", n, (unsigned long long)g_ndag);
    for (int s=0;s<=K;s++){ size_t c=0; for (int i=0;i<NSHARD;i++) c += sets[s]->sh[i].cnt; printf("%d %d %zu\n", n, s, c); fflush(stdout); }
    return 0;
}
