"""Use the sibling pinned-input pipeline without machine-specific paths."""
from pathlib import Path
from fractions import Fraction
import gzip
import json
import os
import pickle
import sys

HERE = Path(__file__).resolve().parent
POSITIVE = Path(os.environ.get('FORTUITY_POSITIVE_SECTOR', str(HERE.parent / 'positive_sector'))).resolve()
if not (POSITIVE / 'probe_stable.py').is_file():
    raise FileNotFoundError('Expected the companion positive_sector/probe_stable.py; set FORTUITY_POSITIVE_SECTOR if it is elsewhere.')
sys.path.insert(0, str(POSITIVE))
from probe_stable import m, np, nmod_mat


def load_stable_rows():
    path = POSITIVE / 'stable_rows.pkl'
    if not path.is_file():
        raise FileNotFoundError('Run the companion positive_sector stable-trace regeneration before the native rank replays.')
    return pickle.loads(path.read_bytes())


def load_omega():
    """Read the portable rational cochain, keeping every original coefficient."""
    path = POSITIVE / 'omega_rational.json.gz'
    if not path.is_file():
        raise FileNotFoundError('Run the companion positive_sector rational_omega.py export first.')
    data = json.loads(gzip.decompress(path.read_bytes()))
    assert data['weight'] == [3, 3, 5, 0, 0] and data['length'] == 9
    omega = {m.decode(word): Fraction(coefficient) for coefficient, word in data['terms']}
    assert len(omega) == 4264
    return omega
