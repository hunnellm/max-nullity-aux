"""Incomplete exact searches: exhaustion is never a nonexistence certificate."""
import random
from fractions import Fraction
from collections import Counter
from .graphs import normalize
from .algebra import Algebra,primitive,pairing
from .certificates import make_certificate,verify_certificate

def _order(G,r,rng,policy):
    left=list(range(G.n));removed=[]
    while left:
        degree={v:sum(not(G.adjacency[v]>>w&1) for w in left if w!=v) for v in left}
        if policy=='min_degree': candidates=[v for v in left if degree[v]==min(degree.values())]
        elif policy=='eligible':
            candidates=[v for v in left if degree[v]<r]
            if not candidates:return None
        else:raise ValueError('order_policy must be min_degree or eligible')
        v=rng.choice(candidates);removed.append(v);left.remove(v)
    return removed[::-1]

def search(graph, rank, *, signature=None, trials=100, seed=0, backend='auto',
           order_policy='min_degree', coefficients=(-2,-1,1,2), vertex_retries=8,
           deterministic_fallback=True, vertex_order=None):
    """Find A=V^T D V of rank <= rank; default D=I.

    signature=(p,q), p+q=rank. Returns {'status':'found','certificate':...}
    or {'status':'exhausted','certificate':None}; no assertion of impossibility.
    The graph's indexed vertex order and seed are part of reproducibility.
    """
    G=normalize(graph,vertex_order);r=int(rank)
    if r!=rank or r<0:raise ValueError('rank must be a nonnegative integer')
    if trials<0 or vertex_retries<1:raise ValueError('invalid search budget')
    if signature is None:p,q=r,0
    else:
        p,q=signature
        if int(p)!=p or int(q)!=q or p<0 or q<0 or p+q!=r:raise ValueError('signature must be nonnegative and sum to rank')
        p,q=int(p),int(q)
    D=[1]*p+[-1]*q;coefficients=tuple(Fraction(x) for x in coefficients)
    if not coefficients or not any(coefficients):raise ValueError('coefficient set must contain a nonzero value')
    if r==0:
        if any(G.adjacency):return dict(status='exhausted',certificate=None,statistics={'rank_zero_obstruction':1})
        c=make_certificate(G,[[] for _ in range(G.n)],[]);return dict(status='found',certificate=c,statistics={})
    alg=Algebra(backend);rng=random.Random(seed);stats=Counter()
    for trial in range(trials):
        stats['trials']+=1;order=_order(G,r,rng,order_policy)
        if order is None:stats['order_obstruction']+=1;continue
        vectors={};trace=[]
        for v in order:
            if not G.adjacency[v]:vectors[v]=[Fraction(0)]*r;trace.append({'vertex':v,'isolated':True});continue
            prior=list(vectors);nonedges=[u for u in prior if not(G.adjacency[v]>>u&1)];edges=[u for u in prior if G.adjacency[v]>>u&1]
            basis=alg.kernel([[D[k]*vectors[u][k] for k in range(r)] for u in nonedges],r)
            if not basis:stats['zero_kernel']+=1;break
            f=[[pairing(vectors[u],b,D) for b in basis] for u in edges]
            if any(not any(row) for row in f):stats['forced_zero_edge']+=1;break
            def candidate(cs):
                return primitive([sum((cs[j]*basis[j][i] for j in range(len(basis))),Fraction(0)) for i in range(r)])
            def allowed(x):return any(x) and all(pairing(vectors[u],x,D) for u in edges)
            x=None
            for _ in range(vertex_retries):
                y=candidate([rng.choice(coefficients) for b in basis])
                if allowed(y):x=y;break
                stats['rejected_samples']+=1
            if x is None and deterministic_fallback:
                # At most m(d-1) forbidden integer parameters; first coefficient is 1.
                for t in range(len(edges)*(len(basis)-1)+1):
                    y=candidate([Fraction(t**j) for j in range(len(basis))])
                    if allowed(y):x=y;stats['deterministic_fallbacks']+=1;break
                if x is None:raise RuntimeError('local feasibility invariant failed')
            if x is None:stats['sampling_exhausted']+=1;break
            vectors[v]=x;trace.append({'vertex':v,'prior_nonneighbors':nonedges,'kernel_dimension':len(basis)})
        else:
            metadata=dict(seed=seed,trial=trial,backend=alg.name,signature=[p,q],order=order,
                order_policy=order_policy,coefficients=list(map(str,coefficients)),vertex_retries=vertex_retries,
                deterministic_fallback=deterministic_fallback,trace=trace,algorithm_version='0.1.0')
            c=make_certificate(G,[vectors[i] for i in range(G.n)],D,metadata);verify_certificate(c)
            return dict(status='found',certificate=c,statistics=dict(stats))
    return dict(status='exhausted',certificate=None,statistics=dict(stats))
