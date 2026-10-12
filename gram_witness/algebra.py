"""Exact kernels: native Sage when available; otherwise SymPy (optionally FLINT)."""
from fractions import Fraction
from math import gcd,lcm
from functools import reduce

class Algebra:
    def __init__(self,backend='auto'):
        if backend not in ('auto','sage','sympy'):raise ValueError('backend must be auto, sage, or sympy')
        if backend in ('auto','sage'):
            try:
                from sage.all import QQ,matrix
                self.QQ,self.matrix,self.name=QQ,matrix,'sage';return
            except ImportError:
                if backend=='sage':raise ImportError('run with sage -python for the Sage backend')
        import sympy
        self.sp,self.name=sympy,'sympy'
    def kernel(self,rows,ncols):
        if not rows:return [[Fraction(int(i==j)) for i in range(ncols)] for j in range(ncols)]
        if self.name=='sage':
            A=self.matrix(self.QQ,len(rows),ncols,[self.QQ(str(x)) for row in rows for x in row])
            return [[Fraction(str(x)) for x in v] for v in A.right_kernel_matrix().rows()]
        A=self.sp.Matrix(rows)
        return [[Fraction(str(x)) for x in v] for v in A.nullspace()]

def primitive(v):
    den=lcm(*(x.denominator for x in v));v=[int(x*den) for x in v]
    g=reduce(gcd,v,0)
    return [Fraction(x//g) for x in v] if g else [Fraction(0) for _ in v]

def pairing(x,y,D):return sum((d*a*b for d,a,b in zip(D,x,y)),Fraction(0))

def inertia(A):
    """Exact symmetric congruence elimination, including zero-diagonal 2x2 pivots."""
    A=[[Fraction(x) for x in row] for row in A];p=q=0
    while A:
        n=len(A);pivot=next((i for i in range(n) if A[i][i]),None)
        if pivot is not None:
            idx=[pivot]+[i for i in range(n) if i!=pivot];A=[[A[i][j] for j in idx] for i in idx]
            d=A[0][0];p+=d>0;q+=d<0
            A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,n)] for i in range(1,n)]
        else:
            edge=next(((i,j) for j in range(n) for i in range(j) if A[i][j]),None)
            if edge is None:return int(p),int(q),n
            u,v=edge;idx=[u,v]+[i for i in range(n) if i not in edge];A=[[A[i][j] for j in idx] for i in idx];b=A[0][1]
            p+=1;q+=1
            A=[[A[i][j]-(A[i][0]*A[1][j]+A[i][1]*A[0][j])/b for j in range(2,n)] for i in range(2,n)]
    return int(p),int(q),0
