# A first supercharge differential in the SO(7) BMN sector at level 26

The source has BMN charges $(2J;q_1,q_2,q_3)=(3;5,4,4)$ and ten letters. Take the archived rational $SO(7)$ source $F$, replace one spacetime derivative of each spin by each of the three flavour derivatives in the saved order, and multiply by the closed quadratic graviton `w_gen((0,),1)`. The resulting highest-weight polynomial has 1,668 trace terms. Its classical variation is the same six substitutions applied to $Q_0 F$, multiplied by the graviton. The integer colour frame verifies $Q_0 F=0$. The first supercharge image has 94,926 terms of common denominator 21. [Source polynomial and substitutions](../so7_hair_first.py) · [Rational parent](../../so7_cartan_target_20261118/representatives.json.gz)

The required first correction has nine letters and Taylor weights $(1,1;4,3,3)$. Its finite-colour cochain dimension is 25,329. Job 5682281 found a minor of order 12,400 in the preceding classical image with residue 343790 modulo 1000003. Since $Q_0^2=0$, this gives

$$
\operatorname{rank}_{\mathbb Q}Q_0\le 25329-12400=12929.
$$

Job 5681838 evaluated the classical images of 26,633 spanning highest-weight trace polynomials and the proposed first-supercharge image on the same 19,472 coordinates. Their ranks over the finite field are 12,929 and 12,930 before and after adjoining the image. Reduction of the integral evaluations gives rational lower bounds. Adjoining one vector increases a rank by at most one. Hence

$$
\operatorname{rank}_{\mathbb Q}Q_0=12929,\qquad
\operatorname{rank}_{\mathbb Q}(Q_0,Q_1 O_0)=12930,\qquad d_1[O_0]\ne0 .
$$

This contributes one highest-weight pairing from a dressed source. [First-action ranks](first_action/first.json) · [Preceding minor](predecessors/previous.json.gz) · [Slurm accounting](jobs.json)

From the repository root, the first command verifies the recorded hashes, determinants and rank argument. The next two regenerate the evaluations with NumPy, SciPy, python-flint and a C++ compiler. The first-action reduction uses a high-memory node.

```sh
d=algorithms/dual_singlets_20261118
python3 "$d/so7_level26/verify.py"
python3 "$d/so7_hair_first.py" --source algorithms/so7_cartan_target_20261118/representatives.json.gz --out /tmp/so7-26-first --workers 32
python3 "$d/so7_hair_bound.py" --out /tmp/so7-26-bound --workers 32
```
