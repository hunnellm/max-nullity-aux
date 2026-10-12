"""Exact exponential ordinary and positive semidefinite zero forcing."""
from itertools import combinations

def closure(adjacency, initial, psd=False):
    n=len(adjacency); full=(1<<n)-1; blue=initial
    while True:
        old=blue; white=full^blue
        if psd:
            components=[]; left=white
            while left:
                comp=left&-left
                while True:
                    prev=comp
                    for v in range(n):
                        if comp>>v&1: comp |= adjacency[v]&white
                    if comp==prev: break
                components.append(comp); left &= ~comp
        else: components=[white]
        for comp in components:
            for v in range(n):
                if blue>>v&1:
                    w=adjacency[v]&comp&~blue
                    if w and not w&(w-1): blue |= w
        if blue==old: return blue

def forcing_number(adjacency, psd=False):
    n=len(adjacency); full=(1<<n)-1
    for k in range(n+1):
        for S in combinations(range(n),k):
            mask=sum(1<<i for i in S)
            if closure(adjacency,mask,psd)==full:
                return k,list(S)
    raise RuntimeError('unreachable')
