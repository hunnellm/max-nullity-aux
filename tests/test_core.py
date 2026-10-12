import copy,json
from pathlib import Path
import networkx as nx
import pytest
import sympy as sp
from gram_witness import search,verify_certificate,normalize,forcing_number
from gram_witness.graphs import from_graph6
from gram_witness.algebra import inertia

ROOT=Path(__file__).resolve().parents[1]
@pytest.mark.parametrize('n',[0,1,5,63])
def test_graph6(n):
    G=nx.path_graph(n);s=nx.to_graph6_bytes(G,header=False).decode().strip()
    assert normalize(G).graph6()==s
    assert from_graph6(s).adjacency==normalize(G).adjacency

def test_inertia_zero_diagonal_and_degenerate():
    assert inertia([[0,2],[2,0]])==(1,1,0)
    assert inertia([[1,1],[1,1]])==(1,0,1)
    assert inertia([[0,0],[0,0]])==(0,0,2)
    # A nonzero isotropic column in an indefinite space can have zero Gram rank.
    c=dict(graph6='@',D=[1,-1],V=[['1'],['1']],A=[['0']],rank=0)
    assert verify_certificate(c)['nullity']==1

def test_certificate_search_and_mutation():
    a=search(nx.path_graph(4),3,seed=2,backend='sympy');b=search(nx.path_graph(4),3,seed=2,backend='sympy')
    assert a==b and a['status']=='found'
    assert verify_certificate(a['certificate'])['psd']
    c=copy.deepcopy(a['certificate']);c['A'][0][1]='0'
    with pytest.raises(ValueError):verify_certificate(c)
    assert search(nx.empty_graph(3),0)['certificate']['rank']==0
    assert search(nx.complete_graph(3),0)['status']=='exhausted'

def test_validation_and_labels():
    G=nx.Graph();G.add_edge('left','right')
    c=search(G,1,seed=0,vertex_order=['right','left'])['certificate']
    assert c['vertex_labels_repr']==["'right'","'left'"]
    with pytest.raises(ValueError):search(nx.DiGraph([(0,1)]),2)
    with pytest.raises(ValueError):search(nx.Graph([(0,0)]),2)
    with pytest.raises(ValueError):search(G,2,signature=(3,-1))
    assert search(G,1,trials=0)['status']=='exhausted'

def reference_closure(G,S,psd):
    B=set(S)
    while True:
        old=B.copy();W=set(G)-B
        components=list(nx.connected_components(G.subgraph(W))) if psd else [W]
        additions=set()
        for C in components:
            for u in B:
                ws=set(G[u])&C
                if len(ws)==1:additions|=ws
        B|=additions
        if B==old:return B

def test_forcing_against_independent_reference():
    from itertools import combinations
    # Exhaustive starting sets on every graph of order at most five.
    for G in nx.graph_atlas_g():
        if len(G)>5:continue
        a=normalize(G).adjacency
        for psd in (False,True):
            expected=next(k for k in range(len(G)+1) if any(len(reference_closure(G,S,psd))==len(G) for S in combinations(G,k)))
            assert forcing_number(a,psd)[0]==expected
    G=nx.star_graph(4);a=normalize(G).adjacency
    assert forcing_number(a)[0]==3 and forcing_number(a,True)[0]==1

def test_all_saved_certificates():
    for fn in ('witnesses.json','psd8.json'):
        for c in json.loads((ROOT/'data'/fn).read_text()):
            report=verify_certificate(c)
            # Independent exact rank calculation with a different algorithm/library.
            A=sp.Matrix([[sp.Rational(x) for x in row] for row in c['A']])
            assert A.rank()==report['rank']

def test_necessity_bounds():
    rows=json.loads((ROOT/'data'/'forcing.json').read_text())
    for r in rows:
        G=from_graph6(r['graph6'])
        assert forcing_number(G.adjacency)[0]==r['Z']
        assert forcing_number(G.adjacency,True)[0]==r['Zplus']
    assert [r['Zplus'] for r in rows[-2:]]==[4,4]

def test_displayed_psd_traces():
    for record in json.loads((ROOT/'data'/'psd_forcing_traces.json').read_text()):
        G=nx.from_graph6_bytes(record['graph6'].encode());B=set(record['initial'])
        for step in record['forces']:
            C=set(step['white_component']);u,v=step['source'],step['target']
            assert C in list(nx.connected_components(G.subgraph(set(G)-B)))
            assert u in B and (set(G[u])&C)=={v}
            B.add(v)
        assert len(B)==12

def test_fallback_and_incomplete_order():
    G=nx.path_graph(4)
    seen=False
    for seed in range(12):
        r=search(G,3,seed=seed,coefficients=(0,1),vertex_retries=1)
        if r['statistics'].get('deterministic_fallbacks',0):
            assert verify_certificate(r['certificate'])['valid'];seen=True;break
    assert seen
    # Eligible ordering can fail for structural reasons; that status is not nonexistence.
    r=search(nx.empty_graph(4),1,order_policy='eligible',trials=1)
    assert r['status']=='exhausted' and r['statistics']['order_obstruction']==1
