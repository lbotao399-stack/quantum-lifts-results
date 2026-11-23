"""Compare the two ordered native contractions with the rational classical section."""

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement, product
import json
from pathlib import Path

from classical_score import HERE, TARGET_SCALE, load, rational_classical_row
from exchange_scores import graph_key, selected
from other_crossing import NORMALIZATION, POLES


def native_scores(banks, data, prefix):
    graphs = list(selected(banks))
    poles = {}
    for shard in range(2):
        for i, values in load(data / f"{prefix}-poles-{shard}.json"):
            assert i % 2 == shard and i not in poles
            poles[i] = sum(w * Fraction(v) for w, v in zip(POLES, values, strict=True))
    assert set(poles) == set(range(len(graphs)))
    return {g: poles[i] for i, g in enumerate(graphs)}


def pair(banks, poles, classical):
    result = {"classical_section": str(classical)}
    for side in ("left", "right"):
        assert not banks[side]["ordinary"]
        value = sum(
            (g["w"] * poles[graph_key(g)] for g in banks[side]["banks"]), Fraction(0)
        )
        result[side] = str(value / NORMALIZATION)
    assert result["left"] == result["right"]
    result["corrected"] = str(Fraction(result["left"]) - classical)
    return result


def run(data):
    point_ids, weights, den, _ = rational_classical_row(
        HERE, HERE / "classical-score-row.json.gz", False
    )

    def section(columns):
        assert len(columns) == 720
        assert [i for i, _ in columns] == list(range(720))
        ordered = [columns[i][1] for i in point_ids]
        return [
            Fraction(
                sum(w * x for w, x in zip(weights, col, strict=True)),
                den * TARGET_SCALE,
            )
            for col in zip(*ordered, strict=True)
        ]

    ordinary = load(data / "unreduced.json.gz")
    assert [tuple(item["indices"]) for item in ordinary] == list(
        product(range(6), repeat=2)
    )
    columns = load(data / "sweep-columns.json.gz")
    assert all(
        len(v) == 36 and i == j and v[3] == a and v[18] == b
        for (i, v), (j, a, b) in zip(
            columns, load(HERE / "reverse-columns.json"), strict=True
        )
    )
    values = section(columns)
    poles = native_scores(ordinary, data, "sweep")
    closed = []
    skipped = []
    for item, classical in zip(ordinary, values, strict=True):
        if any(side["ordinary"] for side in item["banks"].values()):
            skipped.append(item["indices"])
        else:
            closed.append(
                {"indices": item["indices"], **pair(item["banks"], poles, classical)}
            )

    middle = load(data / "middle02.json.gz")
    vectors = middle["kernel"]
    assert middle["indices"] == [0, 2] and middle["ordinary_rank"] + len(vectors) == 20
    assert list(map(tuple, middle["middle_slots"])) == list(
        combinations_with_replacement(range(4), 3)
    )
    control = middle["middle_slots"].index([0, 1, 3])
    varied_columns = load(data / "middle-columns.json.gz")
    assert all(
        i == j and len(raw) == 20 and raw[control] == sweep[2]
        for (i, raw), (j, sweep) in zip(varied_columns, columns, strict=True)
    )
    varied_values = section(varied_columns)
    varied_poles = native_scores(vectors, data, "middle")
    variations = []
    for i, v in enumerate(vectors):
        score = sum(
            x * y for x, y in zip(v["coefficients"], varied_values, strict=True)
        )
        variations.append({"kernel_column": i, **pair(v["banks"], varied_poles, score)})
    return {
        "exact_classical_points": len(columns),
        "fixed_middle_distinct_pole_pairs": len(poles),
        "fixed_middle": closed,
        "fixed_middle_with_unreduced_scalar_terms": skipped,
        "varied_middle_distinct_pole_pairs": len(varied_poles),
        "ordinary_scalar_coordinates": middle["ordinary_coordinates"],
        "ordinary_scalar_rank": middle["ordinary_rank"],
        "varied_middle": variations,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "exchange_data")
    print(json.dumps(run(parser.parse_args().data), indent=2))
