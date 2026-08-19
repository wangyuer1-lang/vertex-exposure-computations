# Computations for the vertex-exposure enumerator papers

Exact-arithmetic computations accompanying three manuscripts on the vertex-exposure
word enumerator

$$E_G=\sum_{\pi\in S_{V(G)}}\bigl[b_1^G(\pi),\dots,b_n^G(\pi)\bigr],
\qquad b_i^G(\pi)=\bigl|N_G(v_i)\cap\{v_1,\dots,v_{i-1}\}\bigr|.$$

1. *Vertex-exposure words of graphs: shuffle factorization, polymatroids, and reconstruction*
2. *Rooted cut laws, branch factorization, and Plücker syzygies for tree reconstruction*
3. *The vertex-exposure enumerator distinguishes all finite trees*

Only paper 2 uses a computation in the statement of a theorem (Theorem 5.2, the
rank certificate through order thirteen). Papers 1 and 3 contain no
computer-assisted proofs; the scripts here were used for development and as
independent checks of the identities proved symbolically in those papers.

## Requirements

Python 3 standard library only. No third-party packages, no build step.
Tested on CPython 3.12 and 3.13.

## The main certificate (paper 2, Theorem 5.2)

```
python root_product_rank_scan.py --max-n 13 \
  --maximum-roots-only --group-jdm \
  --left-excess 3 --connected-left-face --summary-only
```

This constructs, for every free tree of order at most thirteen and every
automorphism orbit of maximum-degree roots, the integer coefficient vector of
the connected cut jet $J_3^\circ(T,r)$, groups rows by order and directed joint
degree matrix, and performs sparse Gaussian elimination modulo
$p = 1{,}000{,}000{,}007$. Every stratum must come out with full row rank
(`deficient=0` on each line). Because the vectors have integer entries, full
row rank mod $p$ implies full row rank over $\mathbb{Q}$.

Expected output, matching Table 1 of paper 2:

```
ROOT-PRODUCT-RANK n=1  groups=1   rows=1    max_group=1  deficient=0
ROOT-PRODUCT-RANK n=2  groups=1   rows=1    max_group=1  deficient=0
ROOT-PRODUCT-RANK n=3  groups=1   rows=1    max_group=1  deficient=0
ROOT-PRODUCT-RANK n=4  groups=2   rows=2    max_group=1  deficient=0
ROOT-PRODUCT-RANK n=5  groups=3   rows=4    max_group=2  deficient=0
ROOT-PRODUCT-RANK n=6  groups=6   rows=7    max_group=2  deficient=0
ROOT-PRODUCT-RANK n=7  groups=11  rows=14   max_group=3  deficient=0
ROOT-PRODUCT-RANK n=8  groups=21  rows=29   max_group=3  deficient=0
ROOT-PRODUCT-RANK n=9  groups=40  rows=63   max_group=5  deficient=0
ROOT-PRODUCT-RANK n=10 groups=79  rows=145  max_group=9  deficient=0
ROOT-PRODUCT-RANK n=11 groups=152 rows=338  max_group=18 deficient=0
ROOT-PRODUCT-RANK n=12 groups=294 rows=811  max_group=39 deficient=0
ROOT-PRODUCT-RANK n=13 groups=570 rows=1970 max_group=78 deficient=0
```

Here `groups` is the number of joint-degree strata, `rows` the number of
maximum-root orbit representatives, and `max_group` the largest stratum.

## File map

**Base layer**

| File | Contents |
|---|---|
| `scan_edge_uncover_n7.py` | dynamic program `uncover_and_lower_dp` computing $E_G$ over vertex subsets |
| `scan_edge_uncover.py` | direct permutation-based enumerator, used to cross-check the DP |
| `tree_u_research.py` | unlabelled free-tree generation, AHU canonical forms, directed joint degree matrix |

**Rooted cut laws (paper 2, Sections 3 and 4)**

| File | Contents |
|---|---|
| `rooted_cut_law_scan.py` | construction of $K_{T,r}$ and scan over rooted types |
| `rooted_cut_recursion_check.py` | root-branch factorization identity |
| `unique_max_cut_factorization_check.py` | unique-maximum-degree slice |
| `max_root_sum_scan.py` | the maximum-root sum identity |

**Rank certificate and syzygies (paper 2, Sections 5 and 6)**

| File | Contents |
|---|---|
| `root_product_rank_scan.py` | JDM-stratified sparse rank certificate (Theorem 5.2) |
| `validate_bounded_root_product.py` | bounded product routine vs. full product then projection |
| `root_product_relation.py` | complete root-product relations |
| `root_product_extrema_scan.py` | counterexamples to the naive leading-term route |
| `cut_jet_relation_witness.py` | the $J_2$ deficient stratum and its depth-three witness |

**Auxiliary checks**

| File | Contents |
|---|---|
| `tree_pseudosimilar_scan.py` | pseudosimilar vertices with isomorphic deletion forests |
| `tree_jdm_obstruction.py` | joint degree matrix obstructions |
| `path_forest_binary_layer_check.py` | path and forest binary layers |

## Cross-checks reported in paper 2, Section 5.3

```
python validate_bounded_root_product.py --max-n 10
python rooted_cut_law_scan.py --max-n 10
python rooted_cut_recursion_check.py --max-n 9
python unique_max_cut_factorization_check.py --max-n 11
python cut_jet_relation_witness.py
```

All coefficients are constructed as exact integers throughout; the only modular
step is the Gaussian elimination in `root_product_rank_scan.py`, and its
conclusion is lifted to $\mathbb{Q}$ as described above.

## Integrity

`SHA256SUMS` lists the SHA-256 hash of every script in this repository. Verify with

```
sha256sum -c SHA256SUMS
```

## License

Not yet chosen. Add a `LICENSE` file before archiving.
