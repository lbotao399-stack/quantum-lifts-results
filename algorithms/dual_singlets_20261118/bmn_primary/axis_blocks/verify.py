"""Replay the USp(4) and G2 rational Taylor homology bounds."""

import json

import paired_axes
import exceptional


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    print(json.dumps([paired_axes.calculate(), exceptional.calculate()], indent=2))
