"""Replay the radial and one-instanton endpoint receipts without changing the originals."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

root = Path(__file__).resolve().parent
checks = [('r007_free_radial.py','r007-free-radial.json'),
          ('r007_full_radial.py','r007-full-radial.json'),
          ('r007_classical_adjoint.py','r007-classical-adjoint.json'),
          ('r007_adjoint_matrix_check.py','r007-adjoint-independent.json'),
          ('r007_boundary_projection.py','r007-boundary-projection.json'),
          ('r007_profile_projection.py','r007-profile-projection.json'),
          ('r007_polarized_propagators.py','r007-polarized-propagators.json')]
with TemporaryDirectory() as tmp:
    for script, name in checks:
        out = Path(tmp)/name
        subprocess.run([sys.executable, str(root/script), '--out', str(out)], check=True,
                       stdout=subprocess.DEVNULL)
        if json.loads(out.read_bytes()) != json.loads((root/name).read_bytes()):
            raise RuntimeError('receipt disagreement: '+name)
        print('verified', name, flush=True)
