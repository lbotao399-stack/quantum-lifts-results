# A neutral-sector test for the first instanton differential

This extends the repository's [one-orientation positivity observation](../../research/22-n4-sym.md#positivity-constraint-on-leading-fourier-sectors) to both orientations and makes explicit which Hamiltonian to test. It is a finite-dimensional conditional reduction, not an instanton calculation.

**Hypotheses.** Use a matched finite near-BPS radial window, with positive inner product and a nilpotent physical supercharge Q. Identify the radial inner products isometrically in a common perturbative frame which does not depend on the vacuum angle theta at order zero. Write t = exp(-8 pi²/g²). Suppose in that frame

    Q = Qp + t K(theta) + O(t²),
    K(theta) = exp(i theta) Q+ + exp(-i theta) Q-.

Here Qp is the completed perturbative differential, not just Q0. The expansion excludes other saddles through the indicated action; a neutral sector at order t would require another term in K. The kernel P of Hp = {Qp,Qp†} must be an isolated cluster. At fixed g it is enough that the supercharge admits an operator-norm first-order expansion as a function of the exponentially small parameters. For a weak-cusp asymptotic, also require that the resolvent and other normalization losses do not spoil exponential ordering (polynomial bounds in 1/g suffice); for angular averages require uniformity in theta. A formal twist computation does not provide these analytic matching conditions.

**Calculation.** Set A = P Q+ P, B = P Q- P acting on the fixed perturbative kernel. Let P_t be the spectral projector of H = {Q,Q†} onto the continued small cluster, and use a near-identity unitary U identifying P with P_t. At fixed gap U = 1 + O(t); for nonuniform gaps interpret each error with its bound in g. Because Q is nilpotent, [Q,H] = [Q†,H] = 0, so P_t reduces both. The exact small-cluster Hamiltonian transported to P is therefore the anticommutator of

    D = P U† Q U P = t [exp(i theta) A + exp(-i theta) B] + O(t²).

The first variation of U contributes no term to D since Qp P = P Qp = 0. Taking the anticommutator gives no term of order t. Its coefficient at order t² is

    H2(theta) = {A,A†} + {B,B†}
               + exp(2i theta) {A,B†} + exp(-2i theta) {B,A†}.

In particular the neutral coefficient is N = {A,A†} + {B,B†}. H2(theta) is positive for each theta, but can have additional zero directions at individual angles because the differential is a sum. Any fixed normalization relating H to the physical distance above the BPS bound must be applied before quoting a dimension shift.

**Trace-only criterion.** Taking the trace on the entire matched perturbative kernel gives

    tr N = 2 ||A||_HS² + 2 ||B||_HS²,

where the Hilbert–Schmidt norm is induced by the radial inner product. Thus a positive coefficient of t² in the neutral Fourier component of the *sum of the matched small-cluster gaps* is equivalent, under these hypotheses, to a nonzero first projected instanton or anti-instanton differential. A vanishing coefficient excludes both on that complete cluster. This is a trace of physical gaps, not a graded index, and individual eigenbranches need not admit angle-independent naming. A charge block need not contain all source and target spaces for a trace converse concerning both maps on a larger window.

**Why naive projection fails.** In even and odd two-dimensional spaces take Q to act only from even to odd, by the matrix

    M = [[t a, t b], [t c, 1]],     a = exp(i theta) - exp(-i theta).

The naive even compression of M†M to the first basis vector is t²(|a|²+|c|²), and the odd compression of MM† is t²(|a|²+|b|²). Both count perturbative-heavy admixtures. The true two small gaps have common leading coefficient |a|²: use det(M†M) = |t a - t² b c|² and the other eigenvalue tending to 1. At theta = 0 this t² coefficient is zero even though A and B are both nonzero; its angular average is 2 in each block. The same issue in an actual Ward computation is resolved by the radial spectral/Gram projection, not by interpreting raw two-point logarithms or unsandwiched variations as projected coefficients.

**Use for R-007.** A one-instanton two-point log of a perturbative BPS survivor need not show the effect even if the first instanton differential is present. A direct projected supercharge matrix element tests it one exponential order earlier than the dimension. Its entries can be defined by the Ward contour of the renormalized supercurrent between perturbatively harmonic radial endpoints; this is a test of finitely many entries, not a determination of the full instanton effective supercharge. An order-t change of identification adds a commutator with Qp to K, which the sandwich P K P kills. This does not allow dropping required contact terms in the Ward correlator. If a physical neutral gap is easier to compute, use the isolated-cluster trace rather than a raw compression and establish the saddle assumptions. This note supplies neither quantity for the SU(3) candidate.
