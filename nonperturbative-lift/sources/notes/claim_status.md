# Research claim map at the end of the overnight continuation

Primary source: /Users/libotao/Desktop/研究项目_非Git/Nonperturbative_lifts_2026-11-11/article/main.tex. Source manuscripts and their normalization remain unchanged.

## Original mechanism and its inputs

The exact-product theorem assumes a compatible exact differential/product and the chain identity {Q,h}=L^N−nu Lambda^(3N), with [Q,L]=0. It then computes the instanton differential of the secondary class ha−L^(N−1)w modulo the specified indeterminacy. The assumed exact chain identity is not replaced by a separated two-point integral. The physical Ward/BV/composite matching is OPEN.

The iterated associator sum is necessary for N>2. Its corrected formula is in homotopy/iterated_multiplication.md. Unconditional claims of strict [L,h]=0 and h²=0 from the failed path-root homotopy were withdrawn; exact finite b-differential homotopy checks through order3 stand with their recorded scope.

## Exact results now supported

- SU2: the C2 partner is rationally verified. The full hidden-product relation C2 star C2+3 C0 star C4_article is a Qpert boundary, and C0 star C4 is nonzero. The complete M4 cocycle and the specified scalar-block uniqueness are certified. The M6-associated leading sector at weight(9,9), (Nb,Nc)=(5,5) has E_infinity dimension1 overQ; the explicit third correction now has an independent saved-state full-Q0 readback, while the fourth/fifth corrections remain uncomputed. This is not a theorem about the C6 generator itself. The all-degree leading two-ghost scalar E2 dimension is1 for even k and0 for odd k, with a universal first primitive from the order-two chain rule; higher-loop survival of the whole tower is not asserted.
- SU3: the original C1 partner is an exact full perturbative cocycle with E_infinity dimension1. Its C0 annihilation and the vanishing C0² indeterminacy image are proved. C1 star C1 is nonzero. Full C3 and C4 cocycles and their complete homogeneous scalar-block uniqueness are certified. All-degree leading scalar E2 dimension is0 for k=2 mod3 and1 otherwise. k6 now has a complete exact second lift and fresh full-Q readback: leading E3 dimension1; no E_infinity assertion.
- SU4: the complete leading tree module is free of rank6 over the Casimir ring. The all-degree Cartan Q1 kernel has dimension0 for k=3 mod4 and1 otherwise. k2 is a full cocycle; k4 and k5 have exact first lifts with independent readbacks. The final bounded k6 first-lift search ended with169 independent source directions and133 selected Q0-image columns, without an exact primitive; checkpoints are retained and no nonexistence is inferred. The all-degree full-E2 upgrade remains conditional on that final exact primitive; its universal chain-rule algebra has independently passed.
- Literature: the full leading SU3/SU4 Cartan representatives equal ordered ordinary Hessians of the known degree(k+2) invariant polynomial twisted periods, with factors +(k+1)²/3 and−(k+1)²/4 relative to the saved generating representatives. Feigin–Silantyev Thm4.8/Cor4.9/Section5 supplies the period classification. The bridge, finite full-Q0 transgressions, higher lifts and physical instanton matching are separate achievements/questions. No exhaustive novelty claim or all-N quantum-complex equivalence is made.

## Physical lane

The complete connection-only M1 block at background/Grassmann coefficient b_H² xi1 xi2 eta1 eta2 e1 e2, including b and ghost jets, is calculated with source symmetrization and exact angular/radial integration. It is nonzero and invariant under the checked h0 change and the specified eta1 homogeneous shift when its finite boundary term is retained. Root independently reran check_normal_connection_radial_readback.py successfully. Its value is in the saved qengine coefficient convention, with component arrays omitting the overall loop factor; it is not the physical instanton coefficient. The remaining minimal operator is the even two-fermion normal gaugino insertion; normal/BV, M0, action/measure and composite/contact contributions remain to be assembled.

## Methods that carry the argument

Finite trace invariant generation and Cayley–Hamilton reductions; genuine Cartan detection supplemented by complete source/target and incoming-map arguments; exact character counts; radial differential equations and coefficient recurrences; finite full-Q0 primitives plus the order-two chain rule; product forests and associator homotopies; modular column selection followed by rational coefficient reconstruction and full polynomial equalities; fresh source-descriptor readback; physical zero modes, regular Green solutions, exact sphere/beta moments and retained finite regulator boundary terms.

## Immediate exact ring consequence

In the original article convention, the certified identities $[C_0][C_2]=0$ and $[C_2]^2=-3[C_0][C_4]\ne0$ imply exact nilpotency index three for $[C_2]$ in the perturbative ring:

$$
[C_2]^3=-3[C_0][C_2][C_4]=0,\qquad [C_2]^2\ne0.
$$

This uses associativity on perturbative cohomology and the verified full quantum-product identity, without division by the loop parameter. The final formal-ring gate is now complete: `rational_c4/C0_C2_boundary_exact.pkl` saves the full primitive w2+hbar*w1+hbar²*w0, with component sizes24,221,105. The construction and independent saved-state full-product/full-Q0 readback both pass, including bare ghosts, gradings and rotations. Certificate SHA256:7c57131cd88155e3fc3e9019332dc5bac9f440d6ebf6af26f586db3907811ca6. The identity holds over Q[hbar] and hence Q[[hbar]]; no hbar division or inference from E_infinity-sector vanishing is used. Thus the displayed exact index-three conclusion is certified in the formal perturbative ring.
