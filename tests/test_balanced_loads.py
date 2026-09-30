"""Semantic boundary controls for the selected balanced-load certificate."""
from pathlib import Path
import copy, gzip, importlib.util, json
import pytest

CAPSULE = Path(__file__).resolve().parents[1] / 'evidence/balanced_loads'
spec = importlib.util.spec_from_file_location('balanced_case', CAPSULE / 'verify_case.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)

@pytest.fixture(scope='module')
def certificate():
    with gzip.open(CAPSULE / 'results/joint_core106.json.gz', 'rt', encoding='utf-8') as stream:
        return json.load(stream)

def test_rank_cannot_understate_the_exact_core(certificate):
    damaged = copy.deepcopy(certificate)
    damaged['rank'] = [0] * len(damaged['rank'])
    with pytest.raises(ValueError, match='rank inequality'):
        checker.verify(CAPSULE, damaged)

def test_signed_lift_must_satisfy_rows_not_only_mass(certificate):
    damaged = copy.deepcopy(certificate)
    # Removing all loads violates a stable-PCR lower row independently of
    # the summary mass comparison.
    damaged['gamma'] = [[0, 1] for _ in damaged['gamma']]
    with pytest.raises(ValueError, match='lift row'):
        checker.verify(CAPSULE, damaged)

def test_return_family_cannot_omit_a_stable_arc(certificate):
    damaged = copy.deepcopy(certificate)
    damaged['families'][0]['walks'].pop()
    with pytest.raises(ValueError, match='complete return family'):
        checker.verify(CAPSULE, damaged)

def test_return_ceiling_cannot_be_weakened(certificate):
    damaged = copy.deepcopy(certificate)
    damaged['families'][0]['walks'][0]['rhs'] -= 1
    with pytest.raises(ValueError, match='ceiling'):
        checker.verify(CAPSULE, damaged)

def test_disconnected_return_is_not_a_circulation_certificate(certificate):
    damaged = copy.deepcopy(certificate)
    damaged['families'][0]['walks'][0]['return'].pop()
    with pytest.raises(ValueError, match='closed|walk'):
        checker.verify(CAPSULE, damaged)
