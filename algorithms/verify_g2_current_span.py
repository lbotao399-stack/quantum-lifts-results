"""Ranks of physical current products in two complete G2 singlet source sectors."""
import argparse
import contextlib
import io
import json
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

from verify_g2_current_product import run as four_letter_source
from verify_g2_s3_currents import probe
from g2_full_classical_screen_v1 import map_ranks, sector


def run():
    classical = []
    with TemporaryDirectory() as folder:
        for half_charge in (2, 3):
            for t in range(half_charge + 1):
                mu = (t, t, *(3 * (half_charge - t,)))
                length = 3 * half_charge - 1 - t
                charge = (1, 1, *(3 * (2 * half_charge,)))
                assert sector(mu, length) == charge
                rows = []
                for prime in (1000003, 1000033):
                    dest = Path(folder) / f'{half_charge}-{t}-{prime}'
                    dest.mkdir()
                    with contextlib.redirect_stdout(io.StringIO()):
                        data = map_ranks(mu, length, 'G2', prime, dest)
                    rows.append([data['dimensions']['0'], data['incoming_rank_lower_bound'],
                                 data['outgoing_rank_lower_bound']])
                assert rows[0] == rows[1]
                dim, incoming, outgoing = rows[0]
                classical.append(dict(charge=charge, weight=mu, length=length, dim=dim,
                                      incoming=incoming, outgoing=outgoing, upper=dim-incoming-outgoing))
    upper = [sum(r['upper'] for r in classical if r['charge'][2] == n) for n in (4, 6)]
    assert upper == [1, 2]
    four = four_letter_source()
    assert four['physical_singlet_cartan_terms'] == [4, 0, 0, 0, 0, 0]
    minor = []
    for matter, mask in ((True, 31), (False, 13)):
        first, second, checks = probe(matter)
        row = [Fraction(first.get(mask, 0)), Fraction(second.get(mask, 0))]
        assert all(c.denominator == 1 for c in row) and checks == 497
        minor.append([int(c) for c in row])
    assert minor == [[-64, 0], [16, 576]]
    determinant = minor[0][0]*minor[1][1] - minor[0][1]*minor[1][0]
    assert determinant == -36864
    return dict(primes=[1000003, 1000033], classical_blocks=classical, source_dimension_upper=upper,
                physical_current_product_values=[[[4]], minor], determinants=[4, determinant],
                matter_free_four_letter_jet=True, matter_term_only_in_first_row_of_six_letter_minor=True,
                dimension_after_quotient_by_current_products=[0, 0],
                physical_input='Exact stress-tensor-multiplet currents and a supersymmetric primary operator product with the usual classical limit.',
                instanton_integral_evaluated=False)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a=p.parse_args();assert __debug__
    z=run()
    a.output.write_text(json.dumps(z, indent=2)+'\n')
    print(z['source_dimension_upper'],z['determinants'],flush=True)
