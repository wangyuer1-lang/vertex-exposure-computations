from __future__ import annotations

import argparse
from collections import Counter, defaultdict

from rooted_cut_law_scan import (
    prefix_suffix_dp,
    rooted_cut_law,
    rooted_orbit_representatives,
)
from rooted_cut_recursion_check import child_product
from tree_u_research import direct_joint_degree_matrix, free_trees
from unique_max_cut_factorization_check import components_after_root, shuffle_words


PRIME = 1_000_000_007


def bounded_child_product(
    adjacency: list[int],
    root: int,
    max_left: int,
    connected_left_face: bool = False,
):
    product = Counter({((), ()): 1})
    for branch, attachment in components_after_root(adjacency, root):
        prefix, suffix = prefix_suffix_dp(branch)
        factor = rooted_cut_law(attachment, prefix, suffix)
        if connected_left_face:
            factor = Counter(
                {
                    pair: coefficient
                    for pair, coefficient in factor.items()
                    if sum(pair[0]) == len(pair[0]) - 1
                }
            )
        updated = Counter()
        for (left_before, left_after), left_coefficient in product.items():
            for (right_before, right_after), right_coefficient in factor.items():
                if len(left_before) + len(right_before) > max_left:
                    continue
                for before, before_multiplicity in shuffle_words(
                    left_before, right_before
                ):
                    for after, after_multiplicity in shuffle_words(
                        left_after, right_after
                    ):
                        updated[(before, after)] += (
                            left_coefficient
                            * right_coefficient
                            * before_multiplicity
                            * after_multiplicity
                        )
        product = updated
    return product


def sparse_rank(rows: list[dict]) -> int:
    basis: dict[tuple, dict] = {}
    for source in rows:
        vector = {
            key: value % PRIME
            for key, value in source.items()
            if value % PRIME
        }
        while vector:
            pivot = max(vector)
            if pivot not in basis:
                inverse = pow(vector[pivot], PRIME - 2, PRIME)
                vector = {
                    key: (value * inverse) % PRIME
                    for key, value in vector.items()
                }
                basis[pivot] = vector
                break
            factor = vector[pivot]
            pivot_row = basis[pivot]
            for key, value in pivot_row.items():
                updated = (vector.get(key, 0) - factor * value) % PRIME
                if updated:
                    vector[key] = updated
                else:
                    vector.pop(key, None)
    return len(basis)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=1)
    parser.add_argument("--max-n", type=int, default=9)
    parser.add_argument("--maximum-roots-only", action="store_true")
    parser.add_argument("--group-degree-sequence", action="store_true")
    parser.add_argument("--group-jdm", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--minimum-left-slice", action="store_true")
    parser.add_argument("--left-excess", type=int)
    parser.add_argument("--connected-left-face", action="store_true")
    args = parser.parse_args()
    trees = free_trees(args.max_n)
    for n in range(args.min_n, args.max_n + 1):
        groups = defaultdict(list)
        for adjacency in trees[n].values():
            maximum = max(row.bit_count() for row in adjacency)
            for root in rooted_orbit_representatives(adjacency):
                degree = adjacency[root].bit_count()
                if args.maximum_roots_only and degree != maximum:
                    continue
                degree_sequence = tuple(
                    sorted((row.bit_count() for row in adjacency), reverse=True)
                )
                if args.group_jdm:
                    key = (
                        degree,
                        degree_sequence,
                        tuple(sorted(direct_joint_degree_matrix(adjacency).items())),
                    )
                elif args.group_degree_sequence:
                    key = (degree, degree_sequence)
                else:
                    key = degree
                if args.minimum_left_slice:
                    product = bounded_child_product(
                        adjacency,
                        root,
                        degree,
                        args.connected_left_face,
                    )
                elif args.left_excess is not None:
                    product = bounded_child_product(
                        adjacency,
                        root,
                        degree + args.left_excess,
                        args.connected_left_face,
                    )
                else:
                    product = child_product(adjacency, root)
                if args.minimum_left_slice:
                    left = (0,) * degree
                    vector = {
                        right: coefficient
                        for (before, right), coefficient in product.items()
                        if before == left
                    }
                elif args.left_excess is not None:
                    vector = {
                        pair: coefficient
                        for pair, coefficient in product.items()
                        if len(pair[0]) <= degree + args.left_excess
                    }
                else:
                    vector = dict(product)
                if args.connected_left_face:
                    vector = {
                        pair: coefficient
                        for pair, coefficient in vector.items()
                        if sum(pair[0]) == len(pair[0]) - degree
                    }
                groups[key].append(vector)
        fields = []
        deficient = []
        for degree, rows in sorted(groups.items()):
            rank = sparse_rank(rows)
            fields.append(f"d={degree}:{rank}/{len(rows)}")
            if rank != len(rows):
                deficient.append((degree, rank, len(rows)))
        if args.summary_only:
            maximum_group = max((len(rows) for rows in groups.values()), default=0)
            print(
                f"ROOT-PRODUCT-RANK n={n} groups={len(groups)} "
                f"rows={sum(map(len, groups.values()))} max_group={maximum_group} "
                f"deficient={len(deficient)}",
                flush=True,
            )
            for key, rank, row_count in deficient:
                print(
                    f"  DEFICIENT key={key!r} rank={rank}/{row_count}",
                    flush=True,
                )
        else:
            print(f"ROOT-PRODUCT-RANK n={n} " + " ".join(fields), flush=True)


if __name__ == "__main__":
    main()
