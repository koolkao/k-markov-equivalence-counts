// kmec.c — count k-Markov equivalence classes of labeled DAGs on n nodes.
// Two DAGs are k-Markov equivalent iff they entail the same d-separation
// statements X _||_ Y | S with |S| <= k  (Kocaoglu 2023 definition).
// Brute force: enumerate every labeled DAG, compute its full separation
// signature, mask to |S|<=k, count distinct masked signatures (exact hash set).
// Usage: ./kmec n      (n <= 7)
// Output: one line per k: n k count
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 8
#define MAXW 16  // 64-bit words per signature (enough for n<=7: 21*32=672 bits)

static int n, npairs, W;
static int pu[64], pv[64];
// triple index list: for pair p and subset S of V\{u,v}
static int ntrip;
static int tp[4096], ts[4096], tsz[4096];

static int popc(unsigned x){ return __builtin_popcount(x); }

// d-separation u _||_ v | S in DAG with parent masks par[]
static int dsep(const unsigned *par, int u, int v, unsigned S) {
    // ancestral closure of {u,v} U S
    unsigned A = (1u<<u)|(1u<<v)|S, prev = 0;
    while (A != prev) { prev = A; for (int i=0;i<n;i++) if (A>>i&1) A |= par[i]; }
    // moral graph restricted to A
    unsigned adj[MAXN] = {0};
    for (int i=0;i<n;i++) if (A>>i&1) {
        unsigned p = par[i] & A;
        for (int j=0;j<n;j++) if (p>>j&1) { adj[i] |= 1u<<j; adj[j] |= 1u<<i; adj[j] |= p & ~(1u<<j); }
    }
    // BFS from u avoiding S
    unsigned allowed = A & ~S, seen = 1u<<u, frontier = 1u<<u;
    while (frontier) {
        unsigned nf = 0;
        for (int i=0;i<n;i++) if (frontier>>i&1) nf |= adj[i];
        nf &= allowed & ~seen;
        if (nf>>v&1) return 0;
        seen |= nf; frontier = nf;
    }
    return 1;
}

// ---- exact hash set of signatures (W words each) ----
typedef struct { uint64_t *keys; uint8_t *used; size_t cap, cnt; int w; } HSet;
static uint64_t hashw(const uint64_t *k, int w){ uint64_t h=1469598103934665603ULL; for(int i=0;i<w;i++){ h^=k[i]; h*=1099511628211ULL; h^=h>>29;} return h; }
static void hs_init(HSet *s, int w){ s->w=w; s->cap=1<<16; s->cnt=0; s->keys=calloc(s->cap*w,8); s->used=calloc(s->cap,1); }
static void hs_insert(HSet *s, const uint64_t *k);
static void hs_grow(HSet *s){
    HSet t; t.w=s->w; t.cap=s->cap*2; t.cnt=0; t.keys=calloc(t.cap*t.w,8); t.used=calloc(t.cap,1);
    for (size_t i=0;i<s->cap;i++) if (s->used[i]) hs_insert(&t, s->keys+i*s->w);
    free(s->keys); free(s->used); *s=t;
}
static void hs_insert(HSet *s, const uint64_t *k){
    if (s->cnt*2 >= s->cap) hs_grow(s);
    size_t i = hashw(k,s->w) & (s->cap-1);
    while (s->used[i]) { if (!memcmp(s->keys+i*s->w,k,8*s->w)) return; i=(i+1)&(s->cap-1); }
    s->used[i]=1; memcpy(s->keys+i*s->w,k,8*s->w); s->cnt++;
}

int main(int argc, char **argv){
    if (argc < 2) { fprintf(stderr,"usage: %s n   (2 <= n <= 6 recommended)\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    if (n < 2 || n > 7) { fprintf(stderr,"2<=n<=7\n"); return 1; }
    npairs=0; for (int i=0;i<n;i++) for (int j=i+1;j<n;j++){ pu[npairs]=i; pv[npairs]=j; npairs++; }
    // triples sorted by |S| so that masking to |S|<=k is a prefix
    ntrip=0;
    for (int sz=0; sz<=n-2; sz++)
        for (int p=0;p<npairs;p++){
            unsigned rest = ((1u<<n)-1) & ~(1u<<pu[p]) & ~(1u<<pv[p]);
            for (unsigned S=0; S<(1u<<n); S++) if ((S & ~rest)==0 && popc(S)==sz){ tp[ntrip]=p; ts[ntrip]=S; tsz[ntrip]=sz; ntrip++; }
        }
    W = (ntrip+63)/64;
    int K = n-2;
    int prefix[MAXN]; // number of triples with |S|<=k
    for (int k=0;k<=K;k++){ prefix[k]=0; for(int t=0;t<ntrip;t++) if(tsz[t]<=k) prefix[k]=t+1; }
    HSet sets[MAXN]; for (int k=0;k<=K;k++) hs_init(&sets[k], (prefix[k]+63)/64);

    // enumerate DAGs: each pair is none / i->j / j->i ; reject cycles
    uint64_t total = 1; for (int p=0;p<npairs;p++) total*=3;
    uint64_t ndag = 0;
    int digit[64]={0};
    for (uint64_t code=0; code<total; code++){
        if (code) { int p=0; while (++digit[p]==3){ digit[p]=0; p++; } }
        unsigned par[MAXN]={0};
        for (int p=0;p<npairs;p++){ if(digit[p]==1) par[pv[p]]|=1u<<pu[p]; else if(digit[p]==2) par[pu[p]]|=1u<<pv[p]; }
        // acyclicity via repeated source removal
        unsigned left=(1u<<n)-1; int ok=1;
        while (left){ unsigned src=0; for(int i=0;i<n;i++) if((left>>i&1) && !(par[i]&left)) src|=1u<<i; if(!src){ok=0;break;} left&=~src; }
        if (!ok) continue;
        ndag++;
        uint64_t sig[MAXW]={0};
        for (int t=0;t<ntrip;t++) if (dsep(par,pu[tp[t]],pv[tp[t]],ts[t])) sig[t>>6] |= 1ULL<<(t&63);
        for (int k=0;k<=K;k++){
            uint64_t m[MAXW]; int w=sets[k].w; memcpy(m,sig,8*w);
            int r = prefix[k] & 63; if (r) m[w-1] &= (1ULL<<r)-1;
            hs_insert(&sets[k], m);
        }
    }
    fprintf(stderr,"n=%d labeled DAGs=%llu\n", n, (unsigned long long)ndag);
    for (int k=0;k<=K;k++) printf("%d %d %zu\n", n, k, sets[k].cnt);
    return 0;
}
