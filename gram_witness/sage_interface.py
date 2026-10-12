"""Thin native-Sage wrapper; no SymPy needed for this interface."""
from .search import search
from .certificates import verify_certificate

def sage_search(G, rank, **kwargs):
    from sage.all import QQ,matrix,diagonal_matrix
    kwargs.setdefault('backend','sage');result=search(G,rank,**kwargs)
    if result['certificate'] is not None:
        c=result['certificate'];n=len(c['A']);r=len(c['D'])
        result['A']=matrix(QQ,n,n,[QQ(x) for row in c['A'] for x in row])
        result['V']=matrix(QQ,r,n,[QQ(x) for row in c['V'] for x in row])
        result['D']=diagonal_matrix(QQ,c['D'])
    return result

def sage_verify(c):
    """Independent native Sage matrix/rank check, plus exact inertia verification."""
    from sage.all import QQ,matrix,diagonal_matrix
    report=verify_certificate(c);n=len(c['A']);r=len(c['D'])
    A=matrix(QQ,n,n,[QQ(x) for row in c['A'] for x in row]);V=matrix(QQ,r,n,[QQ(x) for row in c['V'] for x in row]);D=diagonal_matrix(QQ,c['D'])
    if A!=V.transpose()*D*V or A.rank()!=c['rank']:raise ValueError('native Sage verification failed')
    return report
