"""Check both inclusions in the scalar support calculation on the chosen SU(3) Cartan."""

from hashlib import sha256
from io import BytesIO
import json
from math import comb
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from discriminant_pencil import verify as pencil
from verify import blob

PIN = "8ba864d4eb9413131922b75aa653ab119d10a04a"


def checked_archive_member(archive, suffix, digest):
    for path in archive.namelist():
        if path.endswith("/" + suffix):
            data = archive.read(path)
            if sha256(data).hexdigest() == digest:
                return data
    raise ValueError("missing pinned member " + suffix)


def verify():
    assert __debug__
    directions = pencil()
    assert directions["directions"] == (0, 1, 2, 3, 4, 6)
    clearing = {}
    for power in (6, 70):
        coeff = [(-1) ** k * comb(power + k - 1, k) for k in range(4)]
        product = [
            sum(coeff[j] * comb(power, k - j) for j in range(k + 1)) for k in range(4)
        ]
        assert product == [1, 0, 0, 0]
        clearing[str(power)] = {
            "degree": 3,
            "inverse_coefficients": coeff,
            "uniform_power": power + 3,
        }
    with TemporaryDirectory() as temp:
        dest = Path(temp)
        with ZipFile(
            BytesIO(blob("dot-r007-exact-input/phase34-public-q0-q1-replay.zip"))
        ) as source:
            bare = checked_archive_member(
                source,
                "dmk_o330_hw.json",
                "8bee15e4f741312ca2c05e711b6aec8afcd69fe47dddf78e7671a05f5bab9f34",
            )
            (dest / "source.json").write_bytes(bare)
        with ZipFile(
            BytesIO(blob("dot-j24-formal-completion/j24-phase7-public-full-chain.zip"))
        ) as lower:
            name = "o0_classical_certificate.json"
            frozen = checked_archive_member(
                lower,
                name,
                "78c4eab761234f9c136a42c1481a3c92cc03887e8345545c0e6810995e019eda",
            )
            (dest / name).write_bytes(frozen)
            verifier = "verify_o0_dual.py"
            (dest / verifier).write_bytes(
                lower.read("j24-phase7-public-full-chain/" + verifier)
            )
        (dest / "block.py").write_bytes(blob("dot-r007-exact-input/block_witness.py"))
        dual = json.loads(
            subprocess.check_output([sys.executable, "-B", verifier], cwd=dest)
        )
        block = json.loads(
            subprocess.check_output(
                [sys.executable, "-B", "block.py", "source.json", name], cwd=dest
            )
        )
    assert (
        dual["full_source_columns_replayed"] == 2808 and dual["dual_pairing"] == "-30"
    )
    assert block["real_imag_pairing"] == ["0", "31104"]
    assert block["literal_multiple_of_saved_O0"] == ["0", "-5184/5"]
    # The target block is balanced. The shift only removes elementary scalars.
    # Select all possible background flavour counts that can reach that block.
    source_weight = (7, 4, 4)
    target_weight = (4, 4, 4)
    choices = [(i, j, 3 - i - j) for i in range(4) for j in range(4 - i)]
    contributing = [
        v
        for v in choices
        if tuple(a - b for a, b in zip(source_weight, v)) == target_weight
    ]
    assert contributing == [(3, 0, 0)]
    return dict(
        pinned_block_input=PIN,
        source_terms=directions["source_terms"],
        lower_page_primitive_checks=dual["full_source_columns_replayed"],
        dual_on_frozen_SU2_class=dual["dual_pairing"],
        literal_block_multiplier_imaginary=block["literal_multiple_of_saved_O0"][1],
        collision_dual_imaginary=31104,
        collision_dual_exponents=contributing[0],
        clearing=clearing,
        directions=directions["directions"],
        secant_determinants=[x["determinant"] for x in directions["matching_tests"]],
        five_probe_counterexample=directions["five_test_control"],
        physical_instanton_pairing="open",
    )


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
