"""Independent Python recount (own d-separation and canonicalization) of labeled and unlabeled k-MEC counts, n <= 5. Written by an independent reviewer."""
import itertools, sys
def run(n):
    V=range(n); pairs=list(itertools.combinations(V,2))
    def dconn(par,ch,u,S):
        # path-state reachability (node, arrived-via 'up'/'down'), independent of the C code
        anc=set(S); st=list(S)
        while st:
            x=st.pop()
            for p in par[x]:
                if p not in anc: anc.add(p); st.append(p)
        seen=set(); st=[(u,'up')]; R=set()
        while st:
            x,d=st.pop()
            if (x,d) in seen: continue
            seen.add((x,d))
            if x not in S: R.add(x)
            if d=='up' and x not in S:
                for p in par[x]: st.append((p,'up'))
                for c in ch[x]: st.append((c,'down'))
            elif d=='down':
                if x not in S:
                    for c in ch[x]: st.append((c,'down'))
                if x in anc:
                    for p in par[x]: st.append((p,'up'))
        return R
    trip=[(u,v,S) for sz in range(n-1) for u,v in pairs for S in itertools.combinations([w for w in V if w not in(u,v)],sz)]
    perms=list(itertools.permutations(V))
    lab={k:set() for k in range(n-1)}
    for choice in itertools.product((0,1,2),repeat=len(pairs)):
        par=[set() for _ in V]; ch=[set() for _ in V]
        for (u,v),c in zip(pairs,choice):
            if c==1: par[v].add(u); ch[u].add(v)
            elif c==2: par[u].add(v); ch[v].add(u)
        left=set(V); ok=True
        while left:
            src=[x for x in left if not (par[x]&left)]
            if not src: ok=False;break
            left-=set(src)
        if not ok: continue
        cache={}
        sig=[]
        for u,v,S in trip:
            key=(u,S)
            if key not in cache: cache[key]=dconn(par,ch,u,set(S))
            if v not in cache[key]: sig.append((u,v,S))
        for k in lab: lab[k].add(frozenset(t for t in sig if len(t[2])<=k))
    for k in lab:
        can=set()
        for s in lab[k]:
            best=None
            for p in perms:
                t=tuple(sorted((min(p[u],p[v]),max(p[u],p[v]),tuple(sorted(p[x] for x in S))) for u,v,S in s))
                if best is None or t<best: best=t
            can.add(best)
        print(n,k,len(lab[k]),len(can),flush=True)
for n in range(2,int(sys.argv[1])+1): run(n)
