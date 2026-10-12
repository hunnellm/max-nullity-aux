"""Portable rational certificates; verifier has no optional dependencies."""
from fractions import Fraction
from .graphs import normalize
from .algebra import inertia,pairing

def make_certificate(G,vectors,D,metadata=None):
    n=G.n;r=len(D)
    A=[[pairing(vectors[i],vectors[j],D) for j in range(n)] for i in range(n)]
    ip,im,iz=inertia(A)
    return dict(schema_version=1,graph6=G.graph6(),vertex_labels_repr=[repr(x) for x in G.labels],
        D=list(D),V=[[str(vectors[j][i]) for j in range(n)] for i in range(r)],
        A=[[str(x) for x in row] for row in A],rank=ip+im,inertia=[ip,im,iz],metadata=metadata or {})

def verify_certificate(certificate, graph=None):
    c=certificate;G=normalize(c['graph6']);n=G.n
    if graph is not None and normalize(graph).adjacency!=G.adjacency:raise ValueError('graph indexing/pattern mismatch')
    D=[Fraction(x) for x in c['D']];V=[[Fraction(x) for x in row] for row in c['V']];A=[[Fraction(x) for x in row] for row in c['A']]
    if any(d not in (1,-1) for d in D):raise ValueError('D entries must be +1 or -1')
    if len(V)!=len(D) or any(len(row)!=n for row in V):raise ValueError('incorrect V dimensions')
    if len(A)!=n or any(len(row)!=n for row in A):raise ValueError('incorrect A dimensions')
    vs=[[V[k][i] for k in range(len(D))] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j]!=A[j][i] or A[i][j]!=pairing(vs[i],vs[j],D):raise ValueError('factorization or symmetry failure')
            if i!=j and bool(A[i][j])!=bool(G.adjacency[i]>>j&1):raise ValueError('off-diagonal pattern failure')
    ip,im,iz=inertia(A)
    if c['rank']!=ip+im:raise ValueError('rank mismatch')
    if 'inertia' in c and c['inertia']!=[ip,im,iz]:raise ValueError('inertia mismatch')
    return dict(valid=True,rank=ip+im,nullity=iz,inertia=[ip,im,iz],psd=im==0)
