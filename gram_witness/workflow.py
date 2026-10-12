"""Bound-guided search; exact forcing is exponential and opt-in for large graphs."""
from .graphs import normalize
from .forcing import forcing_number
from .search import search

def analyze(graph, target_rank, *, trials=100, seed=0, backend='auto',
            signatures=None, forcing_limit=18):
    G=normalize(graph)
    if G.n>forcing_limit:raise ValueError('exact forcing budget exceeded; call search directly or raise forcing_limit')
    z,zset=forcing_number(G.adjacency);zp,zpset=forcing_number(G.adjacency,True)
    report=dict(Z=z,Zplus=zp,zero_forcing_set=zset,psd_forcing_set=zpset,
                general_rank_lower=G.n-z,psd_rank_lower=G.n-zp,attempts=[])
    if target_rank<G.n-z:
        report['status']='ruled_out_by_Z';return report
    if target_rank>=G.n-zp:
        a=search(G,target_rank,trials=trials,seed=seed,backend=backend)
        report['attempts'].append(a)
        if a['certificate'] is not None:report['status']='found_psd';return report
    else:report['psd_ruled_out']=True
    if signatures is None:
        # Balanced first; q and r-q are equivalent for graph pattern after negation.
        signatures=[(target_rank-q,q) for q in range(target_rank//2,0,-1)]
    for k,(p,q) in enumerate(signatures):
        a=search(G,target_rank,signature=(p,q),trials=trials,seed=seed+k+1,backend=backend,
                 order_policy='eligible',coefficients=(-1,0,1))
        report['attempts'].append(a)
        if a['certificate'] is not None:report['status']='found_symmetric';return report
    report['status']='exhausted';return report
