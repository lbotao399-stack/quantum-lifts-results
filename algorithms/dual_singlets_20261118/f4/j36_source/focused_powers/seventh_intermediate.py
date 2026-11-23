"""Test the classical cokernel at the first correction to the seventh glueball."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "eight_albert_review"))
from first_page_tower import Equations, det_mod, independent_multiplicities
from source_sizes import multiplicity, trace


def calculate():
    # Horizontal Pieri strip of length three on the two equal dotted rows.
    diagrams = [(9 - k, 6, k) if k else (9, 6) for k in range(4)]
    terminal = [multiplicity(shape, 15, trace) for shape in diagrams]
    model = Equations(7)
    blocks = independent_multiplicities(model)
    packet = json.loads((HERE.parent / "eight_albert_review/tower_probes.json").read_bytes())["7"]
    matrix = [model.row(probe)[0] for probe in packet["probes"]]
    assert blocks == [43, 16, 9] and len(matrix) == len(matrix[0]) == sum(blocks)
    determinants = {p: det_mod(matrix, p) for p in (1009, 1013, 1000003)}
    assert all(determinants.values())
    return dict(terminal_diagrams=diagrams, terminal_multiplicities=terminal,
                intermediate_classical_rank=sum(terminal), terminal_cokernel=0,
                homogeneous_block_dimensions=blocks, minor_orders=len(matrix),
                ordered_homogeneous_determinants=determinants)


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2))
