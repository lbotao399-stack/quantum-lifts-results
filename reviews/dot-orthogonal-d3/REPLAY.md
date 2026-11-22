# Replay instructions and source circuits

Download a text-encoded archive: [SO(11), j=36](so11-j36-replay.tar.gz.b64), [SO(12), j=42](so12-j42-replay.tar.gz.b64), [SO(13), j=42](so13-j42-replay.tar.gz.b64). Save the raw file (not the GitHub HTML page). Decode it as, for example:

```sh
python3 -c 'import base64,sys; sys.stdout.buffer.write(base64.b64decode(sys.stdin.buffer.read()))' < so12-j42-replay.tar.gz.b64 > so12-j42-replay.tar.gz
sha256sum so12-j42-replay.tar.gz
mkdir replay-so12 && tar -xzf so12-j42-replay.tar.gz -C replay-so12
cd replay-so12
bash replay.sh fast
```

The source endpoint files and cochain-construction instructions are inside each bundle. The branch [review note](README.md) records the F4 correction and the open research question.

# Exact native third differentials for three orthogonal gauge groups

These bundles replay chain certificates in a fixed relative holomorphic superjet scheme. They do not identify this differential with the untwisted physical supercharge.

| Gauge group | Grade | Target | Normalized differential on the specified source |
| --- | ---: | --- | ---: |
| SO(11) | 36 | `S^6` | `+3072000 S^6` |
| SO(12), O-even | 42 | `S^7` | `-33177600 S^7` |
| SO(13) | 42 | `S^7` | `-921600 S^7` |

Coefficients use the ordered (P^6), oriented marked fields, defining trace and (Q_1) normalization in the individual bundles, with target dual (λ(S^p)=1). They depend on the written source normalization. In particular, the SO(11) numerator is 1572864000 after the projector scale cancels; dividing by 2 × 4 × 64 = 512 gives 3072000.

## Sources and corrections

The files `seed*.py` specify the defining endpoints, scalar and odd labels, and fixed edges of the signed Gram circuits. The *full* alternating sums specialize to zero because they are determinants of size (N+1) in defining dimension (N). Let (R) denote the specified stable sum of those circuits. A graded reversal homotopy (h) satisfies (Q_0hR=R), so the finite specialization of (x_0^s=P^6ΠhR) is a classical cycle.

- **SO(11):** (R=G_{cycle}-5G_{path3}+(500/3)G_{path3alt}). The full stable moment of (Q_1x_0^s) is exactly zero, with no finite correction tail.
- **SO(12):** (R) is the single three-path circuit (including its fixed crossed ghost edge). For the first correction the full stable moment equals that of `F = Π(374400 Tail14 − 187200 Tail146)`. Both tails are literal finite Gram identities. The corrected replay includes the previously omitted five- and six-factor moments and checks the full fit, including cancellation of eight-factor terms.
- **SO(13):** (R=G_{cycle}+(5/2)G_{3path}) uses two promoted size-14 determinants. The moment fit is exactly (2400Γ(T)) across 122 coordinates; (T) is the literal unpadded saturated relation. A single promoted seed does not pass the first correction gate.

In each case take (x_1^s=-H(Q_1x_0^s-F)) for the stated finite relation (F), then (x_2^s=-H(Q_1x_1^s)). The even stable remainder cancels in the first equation by the full moment fit; the second right-hand side has odd total arity. Finite specialization gives (Q_0x_1=-Q_1x_0), (Q_0x_2=-Q_1x_1). The finite dual identities then give (λ Q_1 x_2=ζ Q_1 x_0). The lower-page ranks and genuine (E_3) target are supplied by the exact finite-relation input described in each bundle, rather than inferred from a nonzero native score or from classical ancestry.

## Replay

Extract a bundle, enter its directory, and run `sha256sum -c MANIFEST.sha256` before replaying. Run `bash replay.sh fast` for rational checks of the archived exact moment rows, pivot fit and normalization, plus native map-order controls. Run `bash replay.sh full` to compile the endpoint enumerator and regenerate the raw partitions and native traces before the same checks. The full computations need Bash, Python 3.8+, Node.js 16+ and a C++17 compiler. Generated logs need not be byte-identical to archived controls; the assertions concern exact coefficient identities. The replays do not re-enumerate every historical lower-page rank-bank generator unless separately noted.

d55e71942d7e3f72db765b9f97c31b0e234b0afcd3b70b2e1b8fbf011587bc2b  so11-j36-replay.tar.gz
fd0f6d814f0d601b25df7b80715e7e67b9d9da674b7bfa8a5e938bc636ca152d  so12-j42-replay.tar.gz
66e1080b0f98d83688e0545730aa1103fd0ea435e8a8b8c81077f4185c61bfeb  so13-j42-replay.tar.gz
