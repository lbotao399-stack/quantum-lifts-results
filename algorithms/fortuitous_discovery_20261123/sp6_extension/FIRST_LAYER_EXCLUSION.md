# First-layer exclusion for the native symplectic trace complex

## Scope and conventions

`Sp(r)` here is the repository's `USp(2r)`. Field length counts native
Taylor letters. `Q0` is the classical differential from
`{Q,Psi}=Psi^2`; it raises field length by one. All cyclic, reversal,
Grassmann and trace conventions are those of the pinned upstream code.

Stable abelianization is applied **before** a finite-rank quotient. It
sorts letters inside each trace with their Grassmann signs, annihilates
odd-length symplectic traces, and retains the original formal trace
generators. It is not evaluation on the finite-rank Cartan, which would
annihilate every finite-rank trace identity and provide no test.

## Statement

Stable abelianization is injective on the entire first nonzero
symplectic trace-identity layer, at field length `2r+2`.

Consequently, that layer has zero intersection with the image of the
universal `Q0`. The finite-rank and universal `Q0` cohomologies coincide
at every field length at most `2r+1`, in every native weight sector.

This is a consequence of the symplectic second fundamental theorem and
the native differential, not a proposed new general invariant-theory
theorem. It does **not** assert injectivity in subsequent identity
layers and does not apply to loop-corrected supercharges.

## Proof

1. The symplectic second fundamental theorem generates tensor relations
   by the Pfaffian antisymmetrizer on `2r+2` fundamental slots. An adjoint
   letter is represented by the symmetric two-index tensor `J X`.
   Marking both indices of one such tensor annihilates the contraction.
   Thus a nonzero trace relation requires at least `2r+2` letters.

2. At exactly `2r+2` letters, every letter has one marked index. The
   unmarked indices pair into `r+1` paths, each containing two letters.
   These paired recipes span the first trace-identity layer. Reversal
   inside a pair and exchange of pairs carry precisely the native
   Grassmann signs. A repeated odd pair is zero and is omitted.

3. For any canonical nonzero paired recipe, exactly one matching of the
   marked endpoints closes every path separately. Its contribution is
   the product of the `r+1` corresponding quadratic traces, with native
   coefficient `(-1)^(r+1)`. Every other matching joins at least two
   paths and has at most `r` trace factors.

4. Different canonical paired recipes give different nonzero products
   of quadratic traces. Therefore their highest trace-number terms are
   linearly independent. In particular, the spanning paired recipes are
   an independent basis of the first identity layer.

5. Stable abelianization preserves each of these quadratic-trace
   products and never increases trace number. It cannot cancel their
   independent highest trace-number terms. Its restriction to the
   first identity layer is consequently injective.

6. The native superfield `Psi` is odd. In the supercommutative letter
   algebra `Psi^2=0`, coefficient by coefficient. Hence stable
   abelianization annihilates every universal `Q0` image. Its
   injectivity on the first identity layer forces any universal
   `Q0` image in that layer to be zero.

7. Below field length `2r+2`, there are no finite-rank trace identities.
   The incoming and central cochain spaces at length at most `2r+1`
   therefore coincide with the universal spaces. At the last of these
   lengths, a purported extra finite-rank cycle would have a universal
   `Q0` variation in the first identity layer. Step 6 forces that
   variation to vanish universally as well. The closed and exact
   spaces, and hence the cohomologies, coincide.

## Consequence for the Sp(6) search

The searched Sp(6) operator has length 15 and native weight
`(3,3,11,0,0)`. At Sp(7), length 15 is `2r+1`, inside the exclusion
range. No fortuitous class at this length can be obtained by restricting
a class from Sp(7). In particular, a certified Sp(6) fortuitous class
at this length has no lift to a Sp(7) class under the trace covering.

The actual adjacent-rank calculation gives:

| Quantity at Sp(7), weight `(3,3,11,0,0)`, identity length 16 | Exact value |
|---|---:|
| Universal highest-weight dimension | 47346 |
| Finite-rank highest-weight dimension | 47342 |
| Complete identity dimension | 4 |
| Rank of `Q0` on the four identities | 1 |
| Rank of stable abelianization on the four identities | 4 |
| A four-by-four abelianization minor | 8 |
| Fortuitous dimension at operator length 15 | 0 |

The three closed relation directions must not be counted as fortuitous
classes: stable abelianization detects all of them. The independent
exact programs and coefficient probes are in
`next_rank_relations.py`, `next_rank_relations.json` and
`next_rank_character.json`.

## Primary mathematical source

Martin Rubey and Bruce W. Westbury, *Combinatorics of symplectic
invariant tensors*, arXiv:1504.02586, especially Definitions 2.10–2.11
and Theorem 2.13:
https://arxiv.org/pdf/1504.02586

The second fundamental theorem is the cited input. The highest
trace-number argument and its `Q0` exclusion consequence above are
spelled out so that the computational claim does not depend on a
numerical rank extrapolation.
