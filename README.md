# graph-gram-witness

Exact rational positive semidefinite and indefinite Gram witnesses for the minimum rank of a simple graph. This is a standalone package. It does not claim to compute minimum rank in general.

For a graph G, search constructs A=V^T D V with the exact off-diagonal zero/nonzero pattern of G. D=I gives PSD witnesses. D=diag(I_p,-I_q) permits indefinite witnesses. An exact certificate proves an upper bound on minimum rank. **A failed search proves nothing about existence.**

## Installation

Regular Python 3.10 or newer:

```sh
python -m pip install -e '.[python]'
# Optional exact FLINT arithmetic used by compatible SymPy installations:
python -m pip install -e '.[fast]'
```

Inside an existing Sage installation:

```sh
sage -pip install -e .
```

The native Sage backend does not require SymPy. Certificate verification and the graph6/forcing routines use only Python's standard library. SymPy is imported only when its kernel backend is requested or auto-selected outside Sage.

## Python

```python
from gram_witness import search, verify_certificate
result = search('KrcmR?DQGI_k', 7, signature=(3,4),
                trials=5000, seed=10122026,
                order_policy='eligible', coefficients=(-1,0,1))
if result['status'] == 'found':
    print(verify_certificate(result['certificate']))
else:
    print('Search exhausted; existence is unresolved.')
```

Graph6 strings, NetworkX simple graphs, and Sage simple graphs are accepted. A supplied graph's current vertex iteration order is used unless `vertex_order` is given. Certificate matrices are indexed 0,...,n-1; `vertex_labels_repr` records the association for display. Arbitrary Python labels are not deserialized from repr. Use the same explicit vertex ordering to compare indexed matrices with a labelled graph.

## Sage

```python
from gram_witness.sage_interface import sage_search, sage_verify
G = graphs.PathGraph(5)
r = sage_search(G, 4, trials=100, seed=12)
if r['certificate'] is not None:
    A, V, D = r['A'], r['V'], r['D']
    assert A == V.transpose()*D*V
    print(A.rank(), sage_verify(r['certificate']))
```

`search(..., backend='sage')` also works. Sage's exact rational `right_kernel_matrix()` supplies the nullspaces. `sage_verify` additionally checks the factorization and rank with native Sage matrices. The default `auto` backend selects Sage if importable, otherwise SymPy. Seeds are reproducible for a fixed backend, version, graph order, coefficient set, and search parameters; different kernel bases can give different searches across backends.

## Recommended bound-guided workflow

```python
from gram_witness.workflow import analyze
report = analyze('KrcmR?DQGI_k', target_rank=7,
                 trials=5000, seed=12)
```

1. Compute ordinary Z and PSD Z_+ exactly, for small graphs.
2. A rank target r<n-Z is impossible for every symmetric witness.
3. A PSD target r<n-Z_+ is impossible; skip that PSD search.
4. Otherwise try PSD first, because its rank and feasibility constraints are simpler.
5. If needed, try balanced indefinite signatures first, followed by more unbalanced signatures. This order is a heuristic, not a theorem.
6. Verify every returned matrix exactly. Never turn `exhausted` into a lower rank bound.

The exact forcing enumeration is exponential and has an 18-vertex default guard in `analyze`. Direct `search` has no forcing dependency. Signature (p,q) means p positive and q negative **ambient** directions, not necessarily the inertia of the returned matrix: the verifier computes actual inertia. The graph6 reader supports the one-byte and 18-bit order encodings; sparse6 and digraph6 are unsupported.

## Mathematical improvements over the exploratory scripts

- Exact local-feasibility test: every required-edge functional must be nonzero on the current kernel.
- Several samples per vertex, followed by a deterministic moment-curve fallback when requested. This guarantees a feasible *local* choice if one exists for the current partial assignment, not a full representation.
- Isolated vertices may receive zero vectors.
- Handles all-zero symmetric blocks and zero-diagonal 2x2 pivots during exact inertia verification.
- Explicit failure statistics and retained seed/order/signature metadata.
- No Sage top-level exploratory code is executed on import.

The search only samples rational representations in a fixed diagonal form. It is not complete over the real numbers, and a greedy partial representation can become trapped even when a global witness exists. Exhausting a finite coefficient set, order policy, signature list, or trial budget is not an impossibility proof.

## Data and proved outcome

`data/witnesses.json`: ten exact certificates from the earlier survivor search.
`data/psd8.json`: PSD rank-eight certificates for the two exceptional PSD cases.
`data/forcing.json` and `data/psd_forcing_traces.json`: bounds and readable forcing traces.
`docs/heuristics_and_necessity.tex`: independently compilable research note, including explicit factors, proofs, historical search order, and a worked general family.

For the graphs `KrcmR@G@YOhC` and `KrcmR?DQGI_k`,

    Z=6, Z_+=4, M_+=4, and 5<=M<=6.

Thus PSD rank <=7 is impossible and their rank-seven symmetric witnesses must be indefinite. Six other survivors have M=M_+=6. The last two have 5<=M_+<=M<=6; no strict separation is asserted for them. All ten have Z-M<=1.

## Verification and tests

```sh
python -m gram_witness verify data/witnesses.json
python -m gram_witness verify data/psd8.json
python -m pip install -e '.[test]'
python -m pytest -q
# In Sage:
sage -pip install pytest
sage -python -m pytest -q tests/test_sage.py
```

The delivered Python run checks independent SymPy ranks, certificate corruption, degeneracy, graph indexing, and agreement of forcing routines with a separate NetworkX implementation on all graphs of order at most five. The Sage-specific test is supplied but **was not run** in the authoring environment, where Sage was unavailable. The GitHub Actions workflow runs Python tests; it makes no claim to validate Sage.

## CLI

```sh
python -m gram_witness search 'KrcmR?DQGI_k' --rank 7 --negative 4 --trials 5000 --seed 12 --output search.json
python -m gram_witness analyze 'KrcmR?DQGI_k' --rank 7 --trials 5000 --seed 12
```

`verify` accepts a raw certificate, a certificate list, or a successful `search` result JSON. The generated `search.json` is therefore directly verifiable.

## Provenance

The investigation was motivated by https://github.com/hunnellm/MZ_matroid (commit cd9d1eb927f7e784233112f33423610eb8793dba) and its witness/lifting and forcing tools. The package is a fresh implementation of the exact bilinear construction, not a port of the entire Sage repository. It does not redistribute upstream source.
