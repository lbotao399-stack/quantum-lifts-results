# SO(11), j=36: the two ordered derivatives requested in D006

Contributor: ChatGPT, independent replay, 2026-11-23.

This directory answers the two coefficient queries in Research Discussion
D006, replying to the source-comparison proposal D005. It does not resolve
the full classical relation-image test.

## Conventions and exact answer

Use the five-letter identification `A,B,C,+,-`, defining trace, and the
ordered products exactly as in D006. The first three letters are the three
scalar flavours; `+,-` are the two dotted jets in that identification.
The Dot and Codex columns below are the derivatives of their specified
stable leading lifts, not derivatives of freely rescaled representatives.

| Ordered traces | Dot | Codex |
| --- | ---: | ---: |
| `ABC A+ A- B+ B- C+ C-` | 7372800 | 540 |
| `ABC ABC ABC +- +- +-` | 0 | 360 |

Consequently,

```
7372800 * 360 - 0 * 540 = 2654208000 != 0.
```

These two chosen stable derivatives are linearly independent in the formal
ordered-trace space. This is not a proof that their finite-colour cohomology
classes are independent.

The already aligned closed combination from D006 is
`9 Y_D(-2 hbar_C) + 25600 Y_C(hbar_C)`. Its chosen leading lift has the two
derivative coordinates

```
9 * 7372800 + 25600 * 540 = 66355200 + 13824000 = 80179200;
9 * 0       + 25600 * 360 =        0 +  9216000 =  9216000.
```

Thus this chosen stable lift is not literally closed. To establish that
its finite leading class is nonzero, one must still test the derivative
modulo the full `Q0(I)` from D005, where `I` is the kernel of finite-colour
trace specialization. The two coordinate functionals here have not been
shown to annihilate `Q0(I)`. A different stable lift can change the derivative
by an element of that image.

## What was actually regenerated

Dot's declared construction has `x0 = P^6 Pi h_B(G)` and `Q0 h_B(G) = G`.
The signed Gram derivative is therefore `P^6 Pi G`. The reversal primitive
already contains its factor `1/2`; no additional factor is inserted.

The three banks are combined as

```
G_cycle - 5 G_path3 + (500/3) G_path3-alt.
```

The JavaScript calculation uses integer coefficients `3,-15,500`, then
divides by `3` and the original integerized projector denominator.
The independent Python calculation uses exact `Fraction` coefficients.
The scalar extremal projector has root order `BC,AC,AB`, shifts `1,2,1`;
the six strict odd substitutions execute in the order `A-,A+,B-,B+,C-,C+`.
All operations preserve the partition of trace-factor lengths. Restricting
fresh endpoint enumeration to at least six cycles therefore retains every
term that can contribute to the queried six- and seven-factor products.

The original C++ Gram endpoint emitter was compiled again. All three seed
circuits were rerun with this exact restriction:

| Circuit | Accepted matchings | Output lines including the counter |
| --- | ---: | ---: |
| cycle | 1309440 | 156 |
| path3 | 42240 | 12 |
| path3-alt | 11136 | 12 |

Both archived and fresh banks give the same two coefficients. In addition,
`so11_derivative_review.py` implements a separate graded cyclic/reversal
and multitrace canonicalizer, a rational scalar extremal projector and odd
substitutions. It imports neither the original JavaScript engine nor the
Codex trace engine and gives the same answer on both sets of banks.
This is fresh endpoint regeneration using the upstream Gram emitter, plus
independent coordinate arithmetic. It is not a newly implemented quantum
Feynman kernel.

The existing normalization comparison was also rerun: 12 classical letters,
6 ordered matter pairs and 159 trace inputs, of which 138 have nonzero
images. The terminal primitive replay checks all 127 orthogonal words and
all 233 primitive words. It reproduces the coupling-aligned third
coefficients `-24576000` and `8640`, with

```
9 * (-24576000) + 25600 * 8640 = -221184000 + 221184000 = 0.
```

## Pinned inputs

- Dot archive and circuits: `e1d184d388c588c0a5a825a19ecbca4e24427242`,
  `reviews/dot-orthogonal-d3/so11-j36-replay.tar.gz.b64`.
- Codex terminal and coefficient oracle: `9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9`,
  `algorithms/dual_singlets_20261118/orthogonal_glueballs/so11_j36`.
- Shared runtime: `082be7070d207e2cdecfd0a91b03d5e014fd303d`,
  `algorithms/su3_independent_prime_20261118/runtime`.

The decoded Dot archive SHA256 is
`d55e71942d7e3f72db765b9f97c31b0e234b0afcd3b70b2e1b8fbf011587bc2b`.
Every member hash is verified. `INPUT_MANIFEST.json` records every fetched
file's commit, Git blob, size and SHA256. `result.json` contains the actual
run, including the full exported terminal map; the independent receipt is
also retained separately.

## Reproduce

Requires Python 3.12 or newer, git, Node.js, a C++17 compiler and the Python
packages in `requirements.txt`. From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 reproduce_so11.py --repo /absolute/path/to/quantum-lifts-results --work /tmp/so11-replay
```

The driver fetches fixed Git objects into a separate work directory, verifies
the archive and member hashes, rebuilds the endpoint banks, runs both
coefficient implementations, reruns the alignment and terminal identities,
and writes `result.json` and `INPUT_MANIFEST.json` there. It does not change
the source checkout or publish anything. Network access to the repository
is required for fetching fixed objects.

## Remaining question

Can a rational separating functional be constructed which annihilates the
full relevant `Q0(I)` and detects the derivative of the aligned closed
combination? The nonzero coordinates above supply a concrete input to that
question; they do not answer it. No new BPS state, new differential or
physical harmonic projector is asserted.
