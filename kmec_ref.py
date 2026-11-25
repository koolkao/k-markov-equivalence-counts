"""Independent reference implementation (networkx d-separation), small n only.
Enumerates DAGs differently (via random-free itertools over edge orientations
then networkx acyclicity check) and uses networkx.is_d_separator."""
import itertools, sys
import networkx as nx

def count(n):
    V = range(n)
    pairs = list(itertools.combinations(V, 2))
    triples = []
    for u, v in pairs:
        rest = [w for w in V if w not in (u, v)]
        for sz in range(len(rest) + 1):
            for S in itertools.combinations(rest, sz):
                triples.append((u, v, frozenset(S)))
    classes = {k: set() for k in range(n - 1)}
    ndag = 0
    for choice in itertools.product((0, 1, 2), repeat=len(pairs)):
        G = nx.DiGraph(); G.add_nodes_from(V)
        for (u, v), c in zip(pairs, choice):
            if c == 1: G.add_edge(u, v)
            elif c == 2: G.add_edge(v, u)
        if not nx.is_directed_acyclic_graph(G):
            continue
        ndag += 1
        stmts = frozenset((u, v, S) for (u, v, S) in triples if nx.is_d_separator(G, {u}, {v}, set(S)))
        for k in classes:
            classes[k].add(frozenset(t for t in stmts if len(t[2]) <= k))
    print(f"n={n} dags={ndag}", file=sys.stderr)
    for k in classes:
        print(n, k, len(classes[k]))

if __name__ == "__main__":
    count(int(sys.argv[1]))
