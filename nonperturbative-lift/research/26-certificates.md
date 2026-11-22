<a id="np-section-24"></a>

## 26. Appendix: electronic certificates

[Download the electronic certificates (ZIP)](../artifacts/nonperturbative-certificates.zip), containing exact witnesses, verification records, and the source manifest.

**Table E1: machine-readable certificates for the pure-SYM calculations.**

Unlinked paths in the table refer to the original certificate ZIP. The four expanded polynomials have separate download links.

Its `manifest.json` records the complete provenance, byte counts, and SHA256 hashes.

PKL files contain exact rational polynomial dictionaries.

JSON files record normalization, equations, rank witnesses, or independent readback.

| Result | Full witness or coefficient data | Verification record |
|---|---|---|
| $C_0\star C_2=Q_{\mathrm{pert}}w$ | `rational_c4/C0_C2_boundary_exact.pkl`, `boundary` field | `rational_c4/C0_C2_boundary_readback.json` |
| $C_2\star C_2+3C_0\star C_4$ full primitive | `rational_c4/full_hidden_product_exact.pkl`, `boundary` field | `rational_c4/full_completion_summary.json` |
| $[C_0][C_4]\ne0$ | `rational_c4/C0C4_nonzero_exact.pkl` | `rational_c4/nonzero_summary.json` |
| SU(2) $C_2$ and first partner | `rational_c4/C2_complete_exact.json` | Source convention and component equations in the same file |
| SU(2) weight $(7,7)$ full partner class | `rational_c4/m4_lift/M4_full_exact.pkl` | `rational_c4/m4_lift/readback.json`, `uniqueness.md` |
| SU(2) all-degree first lift | `su2_transgression/all_degree_first_lift.pkl` | `su2_transgression/all_degree_first_lift.json` |
| SU(2) weight $(9,9)$ leading class | `c6_survival/cartan_target_certificate.json`, `cartan_witnesses_nb2_certificate.json` | `c6_survival/proofs.md` |
| SU(2) third correction | [Expanded third correction](https://github.com/lbotao399-stack/quantum-lifts-results/releases/download/nonperturbative-lift-v1/su2-third-correction.pkl.gz) | `rational_c4/m6_lift/compact_first/lift/level3_exact_summary.json`, `level3_independent_readback.json` |
| SU(3) first partner | `su3_rational/partner_Q.pkl` | `su3_rational/partner_Q_certificate.json`, `incoming_d1_Q_certificate.json` |
| SU(3) indeterminacy and $C_0$ annihilation | [Expanded product](../artifacts/witnesses/su3-c0-product.pkl.gz) and exact rank witnesses in `su3_rational/proofs.md` | `su3_rational/C0_squared_image_certificate.json`, `c0_product_Q_certificate.json` |
| SU(3) $C_1^2$ | `su3_rational/C1_star_Q.pkl`, `C1_square_Q.pkl` | `su3_rational/C1_star_Q_certificate.json`, `C1_square_Q_certificate.json` |
| SU(3) full $C_3$ | `su3_rational/C3_Q.pkl` and the trace appendix | `su3_rational/C3_Q_certificate.json`, `C3_scalar_complex_Q_certificate.json` |
| SU(3) full $C_4$ | `su3_rational/k4_full/C4_Q.pkl` and the trace appendix | `su3_rational/k4_full/independent_full_Q_readback.json`, `scalar_ambiguity_Q_certificate.json` |
| SU(3) $k=6$ second lift | `su3_rational/k6_two_loop/joint_two_loop_Q.pkl` | `su3_rational/k6_two_loop/independent_joint_two_loop_readback.json` |
| SU(3) all-degree $E_2$ | Full module and primitive formulas in `su3_rational/leading_module_and_radial_kernel.md` | `su3_rational/global_leading_E2_Q_certificate.json`, `su3_transgression/verify_universal_first_lift.json` |
| SU(4) tree module and Cartan kernel | Exact reduction data in the verification column | `su4_cartan/leading_tree_module.json`, `leading_radial.json`, `kernel_recurrence.json` |
| SU(4) full $k=2$ cocycle | `su4_cartan/k2_full.pkl` | `su4_cartan/k2_independent_readback.json` |
| SU(4) $k=4,5$ first lifts | [Degree four](../artifacts/witnesses/su4-k4-first-lift.pkl.gz) and [degree five](../artifacts/witnesses/su4-k5-first-lift.pkl.gz) | `su4_cartan/k4_first/first_lift.json`, `independent_readback.json`; `su4_cartan/k5_first/first_lift.json`, `independent_readback.json` |
| SU(4) remaining $k=6$ gate | Status of the open first-lift search | `su4_cartan/k6_first/bounded_status.json` |
| Hessian bridge | Exact determinant derivation in the body | `literature/verify_hessian_bridge.json`, `verify_hessian_bridge.py` |

The [historical file index](../artifacts/historical-files.json) identifies the original calculation programs, coefficient data, and saved intermediate results in the versioned release. The 17 portable calculations have their own [replay catalog](../checks/catalog.json).

The complete $934381$-term SU(2) third correction occupies $686635425$ uncompressed bytes. Its SHA-256 matches the digest in the original exact calculation. The [verification guide](../verification/README.md#external-witnesses) gives the download and byte-level checks for all four polynomials.

The product files use $C_{4,\mathrm{source}}$, while the article uses $C_4=-C_{4,\mathrm{source}}$.

The body gives the old/product normalization bridge for the complete $C_2$ export; retain it when comparing coefficients across files.
