// kmec_unl.c — count k-Markov equivalence classes of DAGs on n nodes up to
// isomorphism (unlabeled), n <= 6, exact (no hashing of signatures).
// Step 1: enumerate labeled DAGs, compute exact signature sig[S] = bitmask over
//         pairs {u,v} (u,v ∉ S) of d-separated pairs, for every S with |S|<=n-2.
// Step 2: per k, collect distinct labeled k-signatures (restriction to |S|<=k).
// Step 3: canonicalize each distinct one as the lexicographic minimum over all
//         n! vertex relabelings; count distinct canonical forms.
// Also counts unlabeled DAGs (canonical form of the adjacency) as a check vs A003087.
// Usage: ./kmec_unl n
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 6
#define MAXS 64
static int n, npairs, pidx[MAXN][MAXN], pu[16], pv[16];
static int nS; static unsigned Slist[MAXS]; static int Spos[1<<MAXN]; static int sizeEnd[MAXN];
static int nperm; static int perm[720][MAXN];
static int permS[720][MAXS];          // position of sigma(S)
static uint16_t permP[720][1<<8], permP2[720][1<<8]; // pair-mask permutation (low/high byte)

typedef struct { uint16_t s[MAXS]; } Sig;

// exact hash set of Sig prefixes of length L
typedef struct { Sig *k; uint8_t *u; size_t cap, cnt; int L; } HS;
static uint64_t hsig(const Sig *a, int L){ uint64_t h=1469598103934665603ULL; for(int i=0;i<L;i++){ h^=a->s[i]; h*=1099511628211ULL; h^=h>>31; } return h; }
static void hs_init(HS *h,int L){ h->L=L; h->cap=1<<12; h->cnt=0; h->k=calloc(h->cap,sizeof(Sig)); h->u=calloc(h->cap,1); }
static int hs_put(HS *h,const Sig *a);
static void hs_grow(HS *h){ HS t; hs_init(&t,h->L); free(t.k); free(t.u); t.cap=h->cap*2; t.k=calloc(t.cap,sizeof(Sig)); t.u=calloc(t.cap,1);
    for(size_t i=0;i<h->cap;i++) if(h->u[i]) hs_put(&t,&h->k[i]); free(h->k); free(h->u); *h=t; }
static int hs_put(HS *h,const Sig *a){
    if(h->cnt*2>=h->cap) hs_grow(h);
    size_t i=hsig(a,h->L)&(h->cap-1);
    while(h->u[i]){ if(!memcmp(h->k[i].s,a->s,2*h->L)) return 0; i=(i+1)&(h->cap-1); }
    h->u[i]=1; memset(&h->k[i],0,sizeof(Sig)); memcpy(h->k[i].s,a->s,2*h->L); h->cnt++; return 1;
}

static int dsep(const unsigned *par,int u,int v,unsigned S){
    unsigned A=(1u<<u)|(1u<<v)|S,prev=0;
    while(A!=prev){prev=A;for(int i=0;i<n;i++) if(A>>i&1) A|=par[i];}
    unsigned adj[MAXN]={0};
    for(int i=0;i<n;i++) if(A>>i&1){ unsigned p=par[i]&A; for(int j=0;j<n;j++) if(p>>j&1){adj[i]|=1u<<j;adj[j]|=1u<<i;adj[j]|=p&~(1u<<j);} }
    unsigned allowed=A&~S,seen=1u<<u,fr=1u<<u;
    while(fr){unsigned nf=0;for(int i=0;i<n;i++) if(fr>>i&1) nf|=adj[i]; nf&=allowed&~seen; if(nf>>v&1) return 0; seen|=nf; fr=nf;}
    return 1;
}

static void gen_perms(int d,int *p,int used){
    if(d==n){ memcpy(perm[nperm++],p,sizeof(int)*n); return; }
    for(int i=0;i<n;i++) if(!(used>>i&1)){ p[d]=i; gen_perms(d+1,p,used|1<<i); }
}
static unsigned mapset(int q,unsigned S){ unsigned r=0; for(int i=0;i<n;i++) if(S>>i&1) r|=1u<<perm[q][i]; return r; }

static void canon(const Sig *a,int L,Sig *best){
    int first=1;
    for(int q=0;q<nperm;q++){
        Sig b; memset(&b,0,sizeof b);
        for(int i=0;i<L;i++){ uint16_t m=a->s[i]; b.s[permS[q][i]] = permP[q][m&255] | permP2[q][m>>8]; }
        if(first || memcmp(b.s,best->s,2*L)<0){ *best=b; first=0; }
    }
}
// memcmp on uint16 is byte-order dependent but consistent, which is all a canonical form needs.

int main(int argc,char**argv){
    if(argc<2){fprintf(stderr,"usage: %s n   (2 <= n <= 6)\n",argv[0]);return 1;}
    n=atoi(argv[1]); if(n<2||n>MAXN){fprintf(stderr,"2<=n<=6\n");return 1;}
    npairs=0; for(int u=0;u<n;u++) for(int v=u+1;v<n;v++){ pidx[u][v]=pidx[v][u]=npairs; pu[npairs]=u; pv[npairs]=v; npairs++; }
    nS=0; for(int s=0;s<=n-2;s++){ for(unsigned S=0;S<(1u<<n);S++) if(__builtin_popcount(S)==s){ Spos[S]=nS; Slist[nS++]=S; } sizeEnd[s]=nS; }
    int p[MAXN]; gen_perms(0,p,0);
    for(int q=0;q<nperm;q++){
        for(int i=0;i<nS;i++) permS[q][i]=Spos[mapset(q,Slist[i])];
        for(int half=0;half<2;half++) for(int m=0;m<256;m++){ uint16_t r=0;
            for(int b=0;b<8;b++) if(m>>b&1){ int pp=half*8+b; if(pp<npairs) r|=1u<<pidx[perm[q][pu[pp]]][perm[q][pv[pp]]]; }
            if(half==0) permP[q][m]=r; else permP2[q][m]=r; }
    }
    int K=n-2;
    HS lab[MAXN], unl[MAXN], dagU; for(int k=0;k<=K;k++){ hs_init(&lab[k],sizeEnd[k]); hs_init(&unl[k],sizeEnd[k]); }
    hs_init(&dagU, n); // canonical adjacency: s[i] = parent mask of node i
    int digit[16]={0}; uint64_t total=1; for(int i=0;i<npairs;i++) total*=3; uint64_t ndag=0;
    for(uint64_t code=0;code<total;code++){
        if(code){int i=0; while(++digit[i]==3){digit[i]=0;i++;}}
        unsigned par[MAXN]={0};
        for(int i=0;i<npairs;i++){ if(digit[i]==1) par[pv[i]]|=1u<<pu[i]; else if(digit[i]==2) par[pu[i]]|=1u<<pv[i]; }
        unsigned left=(1u<<n)-1; int ok=1;
        while(left){unsigned src=0;for(int i=0;i<n;i++) if((left>>i&1)&&!(par[i]&left)) src|=1u<<i; if(!src){ok=0;break;} left&=~src;}
        if(!ok) continue;
        ndag++;
        // canonical DAG (for A003087 check)
        Sig bd; int firstd=1;
        for(int q=0;q<nperm;q++){ Sig b; memset(&b,0,sizeof b);
            for(int x=0;x<n;x++) b.s[perm[q][x]] = (uint16_t)mapset(q,par[x]);
            if(firstd||memcmp(b.s,bd.s,2*n)<0){bd=b;firstd=0;} }
        hs_put(&dagU,&bd);
        Sig s; memset(&s,0,sizeof s);
        for(int i=0;i<nS;i++){ unsigned S=Slist[i]; uint16_t m=0;
            for(int pp=0;pp<npairs;pp++) if(!(S>>pu[pp]&1)&&!(S>>pv[pp]&1)&&dsep(par,pu[pp],pv[pp],S)) m|=1u<<pp;
            s.s[i]=m; }
        for(int k=0;k<=K;k++) hs_put(&lab[k],&s);
    }
    fprintf(stderr,"n=%d labeled DAGs=%llu unlabeled DAGs=%zu\n",n,(unsigned long long)ndag,dagU.cnt);
    for(int k=0;k<=K;k++){
        for(size_t i=0;i<lab[k].cap;i++) if(lab[k].u[i]){ Sig c; canon(&lab[k].k[i],lab[k].L,&c); hs_put(&unl[k],&c); }
        printf("%d %d %zu %zu\n",n,k,lab[k].cnt,unl[k].cnt); // n k labeled unlabeled fflush(stdout);
    }
    return 0;
}
