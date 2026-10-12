"""Graph normalization without a required graph-library dependency."""
from dataclasses import dataclass

@dataclass(frozen=True)
class GraphData:
    adjacency: tuple
    labels: tuple
    @property
    def n(self): return len(self.adjacency)
    @property
    def edges(self):
        return [(i,j) for j in range(self.n) for i in range(j) if self.adjacency[i]>>j&1]
    def graph6(self):
        n=self.n
        if n<63: head=[n]
        elif n<258048: head=[63]+[(n>>k)&63 for k in (12,6,0)]
        else: raise ValueError('graph6 encoding limited to n<258048')
        bits=[(self.adjacency[i]>>j)&1 for j in range(n) for i in range(j)]
        bits += [0]*((-len(bits))%6)
        vals=[sum(bits[t+k]<<(5-k) for k in range(6)) for t in range(0,len(bits),6)]
        return ''.join(chr(x+63) for x in head+vals)

def from_graph6(text):
    text=text.strip().removeprefix('>>graph6<<'); vals=[ord(c)-63 for c in text]
    if not vals or any(x<0 or x>63 for x in vals): raise ValueError('invalid graph6')
    pos=1; n=vals[0]
    if n==63:
        if len(vals)<4: raise ValueError('truncated graph6')
        if vals[1]==63: raise ValueError('36-bit graph6 orders not supported')
        n=(vals[1]<<12)|(vals[2]<<6)|vals[3];pos=4
    needed=(n*(n-1)//2+5)//6
    if len(vals)-pos!=needed: raise ValueError('incorrect graph6 length')
    a=[0]*n;k=0
    for j in range(n):
        for i in range(j):
            bit=(vals[pos+k//6]>>(5-k%6))&1;k+=1
            if bit:a[i]|=1<<j;a[j]|=1<<i
    if needed and k%6 and vals[-1]&((1<<(6-k%6))-1):raise ValueError('nonzero graph6 padding')
    return GraphData(tuple(a),tuple(range(n)))

def normalize(graph, vertex_order=None):
    if isinstance(graph,GraphData):
        if vertex_order is not None:raise ValueError('GraphData is already indexed')
        return graph
    if isinstance(graph,str):
        if vertex_order is not None:raise ValueError('graph6 fixes vertex order')
        return from_graph6(graph)
    if graph.is_directed(): raise ValueError('only undirected graphs supported')
    if hasattr(graph,'is_multigraph') and graph.is_multigraph():raise ValueError('use a simple graph')
    if hasattr(graph,'has_multiple_edges') and graph.has_multiple_edges():raise ValueError('use a simple graph')
    vertices=list(graph.nodes()) if hasattr(graph,'nodes') else list(graph.vertices(sort=False))
    labels=vertices if vertex_order is None else list(vertex_order)
    if len(labels)!=len(vertices) or set(labels)!=set(vertices):raise ValueError('vertex_order must list every vertex exactly once')
    index={v:i for i,v in enumerate(labels)};a=[0]*len(labels)
    for v in labels:
        for w in graph.neighbors(v):
            if v==w:raise ValueError('loops are not supported')
            a[index[v]] |= 1<<index[w]
    return GraphData(tuple(a),tuple(labels))
