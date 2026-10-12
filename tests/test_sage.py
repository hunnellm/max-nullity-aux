"""Run with sage -python -m pytest tests/test_sage.py."""
import json
from pathlib import Path
import pytest
sage=pytest.importorskip('sage.all')
from gram_witness.sage_interface import sage_search,sage_verify
from gram_witness import search

def test_native_sage():
    G=sage.graphs.PathGraph(4)
    a=sage_search(G,3,seed=2,trials=100)
    assert a['certificate'] is not None
    assert a['A']==a['V'].transpose()*a['D']*a['V']
    assert a['A'].rank()<=3
    assert sage_verify(a['certificate'])['psd']
    for fn in ('witnesses.json','psd8.json'):
        for c in json.loads((Path(__file__).resolve().parents[1]/'data'/fn).read_text()):
            assert sage_verify(c)['valid']
