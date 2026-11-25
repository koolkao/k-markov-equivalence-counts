// Cross-check: kmec2.c bitmask Bayes-ball vs moralized-ancestral-graph d-separation on random DAGs, n = 3..8.
// Build: cc -O2 -o tests/test_bayesball tests/test_bayesball.c && ./tests/test_bayesball   (run from the repo root)
#define main k2main
#include "../kmec2.c"
#undef main
static int dsep(const unsigned *par, int u, int v, unsigned S) {
    unsigned A = (1u<<u)|(1u<<v)|S, prev = 0;
    while (A != prev) { prev = A; for (int i=0;i<n;i++) if (A>>i&1) A |= par[i]; }
    unsigned adj[MAXN] = {0};
    for (int i=0;i<n;i++) if (A>>i&1) {
        unsigned p = par[i] & A;
        for (int j=0;j<n;j++) if (p>>j&1) { adj[i] |= 1u<<j; adj[j] |= 1u<<i; adj[j] |= p & ~(1u<<j); }
    }
    unsigned allowed = A & ~S, seen = 1u<<u, frontier = 1u<<u;
    while (frontier) { unsigned nf = 0; for (int i=0;i<n;i++) if (frontier>>i&1) nf |= adj[i];
        nf &= allowed & ~seen; if (nf>>v&1) return 0; seen |= nf; frontier = nf; }
    return 1;
}
int main(void){
  srand(12345); long bad=0, tot=0;
  for (n=3;n<=8;n++) for (int it=0; it<4000; it++){
    int perm[8]; for(int i=0;i<n;i++)perm[i]=i; for(int i=n-1;i>0;i--){int j=rand()%(i+1);int t=perm[i];perm[i]=perm[j];perm[j]=t;}
    int dens = 1+rand()%9; unsigned par[MAXN]={0}, ch[MAXN]={0}, anc[MAXN]={0};
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++) if(rand()%10<dens){ par[perm[j]]|=1u<<perm[i]; ch[perm[i]]|=1u<<perm[j]; }
    for(int i=0;i<n;i++){ unsigned a=1u<<i,p=0; while(a!=p){p=a; for(int x=0;x<n;x++) if(a>>x&1) a|=par[x];} anc[i]=a; }
    for (unsigned S=0;S<(1u<<n);S++){ unsigned A=orpar(anc,S)|S;
      for(int u=0;u<n;u++){ if(S>>u&1)continue; unsigned R=bayesball(par,ch,u,S,A);
        for(int v=0;v<n;v++){ if(v==u||(S>>v&1))continue; tot++; if(((R>>v)&1)==dsep(par,u,v,S)) bad++; } } }
  }
  printf("tests=%ld mismatches=%ld\n",tot,bad); return bad != 0; }
