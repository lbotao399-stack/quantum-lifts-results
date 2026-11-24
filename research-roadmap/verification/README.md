# Public verification tier

This is the designated inexpensive tier for A08. It checks the unchanged
published coefficient dual, the complete independent orthogonal return
identity, the coupled-equation and lossless-packing controls, and the complete
new SU(3) generator certificate. A real copied input is corrupted and must be
rejected. Every input manifest, command, exit status and output hash is recorded.

From a **clean Linux checkout**, with Python 3.12:

```sh
python -m pip install -r research-roadmap/verification/requirements.txt
python -B research-roadmap/verification/run.py --output /tmp/public-verification-receipt
```

The runner refuses a dirty checkout unless explicitly invoked in development
mode. Development mode cannot produce a clean-checkout receipt. Historical
controls that write beside their source run from unchanged temporary copies.
No private paths, cluster credentials or scheduler impersonation are used.
No projection routine is invoked by these entry points.

## What is rebuilt

| Entry point | Exact scope |
|---|---|
| Roadmap validator | Weights, dependencies, pins and generated documentation |
| Sp(3) coefficient dual | The published candidate's exact pairing one; inherited closure and completeness |
| SO(11) global return | All 1000 marker monomials and the independent coefficient-operator bridge |
| Coupled equations | Complete small control, all homogeneous directions and restart |
| Map packing | All control entries, missing-row rejection and a corrupted matrix entry |
| SU(3) generator certificate | Complete rational ranks, relations and duals; archived frame entries reused |
| SU(3) independent polynomial check | The whole printed derivative and all 1588 native residual coordinates |

The corresponding expensive source rebuild is separately callable:

```sh
python -B algorithms/su3_generator_20261124/reproduce.py /tmp/su3-generator-rebuild
```

Historical production replays remain documented in their original bundles.
This inexpensive tier does not claim to reconstruct the unverified remaining
G2 source-minor entries, the dense SU(3) production correction maps, or the
complete F4 source. A08 additionally depends on completion of A07; introducing
this runner alone earns no point.
