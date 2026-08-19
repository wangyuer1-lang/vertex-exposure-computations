from __future__ import annotations

import argparse

from rooted_cut_law_scan import rooted_orbit_representatives
from rooted_cut_recursion_check import child_product
from tree_u_research import direct_joint_degree_matrix, free_trees, rooted_ahu
from unique_max_cut_factorization_check import components_after_root


PRIME = 1_000_000_007


def clean_add(target: dict, source: dict, factor: int) -> None:
    for key, value in source.items():
        updated = (target.get(key, 0) + factor * value) % PRIME
        if updated:
            target[key] = updated
        else:
            target.pop(key, None)


def signed(value: int) -> int:
    return value if value <= PRIME // 2 else value - PRIME


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--degree-sequence", type=str)
    args = parser.parse_args()
    requested_sequence = (
        tuple(map(int, args.degree_sequence.split(",")))
        if args.degree_sequence
        else None
    )
    trees = free_trees(args.n)[args.n]
    rows = []
    metadata = []
    for tree_code, adjacency in trees.items():
        maximum = max(row.bit_count() for row in adjacency)
        degree_sequence = tuple(
            sorted((row.bit_count() for row in adjacency), reverse=True)
        )
        if maximum != args.degree:
            continue
        if requested_sequence is not None and degree_sequence != requested_sequence:
            continue
        for root in rooted_orbit_representatives(adjacency):
            if adjacency[root].bit_count() != args.degree:
                continue
            rooted_code = rooted_ahu(adjacency, root)
            orbit_size = sum(
                rooted_ahu(adjacency, vertex) == rooted_code
                for vertex in range(len(adjacency))
            )
            rows.append(dict(child_product(adjacency, root)))
            branches = tuple(
                sorted(
                    rooted_ahu(branch, attachment)
                    for branch, attachment in components_after_root(adjacency, root)
                )
            )
            joint_degree_matrix = tuple(
                sorted(direct_joint_degree_matrix(adjacency).items())
            )
            metadata.append(
                (
                    tree_code,
                    rooted_code,
                    orbit_size,
                    degree_sequence,
                    branches,
                    joint_degree_matrix,
                )
            )

    basis = {}
    dependencies = []
    for index, source in enumerate(rows):
        vector = {key: value % PRIME for key, value in source.items() if value % PRIME}
        combination = {index: 1}
        while vector:
            pivot = max(vector)
            if pivot not in basis:
                inverse = pow(vector[pivot], PRIME - 2, PRIME)
                vector = {key: value * inverse % PRIME for key, value in vector.items()}
                combination = {
                    key: value * inverse % PRIME
                    for key, value in combination.items()
                }
                basis[pivot] = (vector, combination)
                break
            factor = vector[pivot]
            pivot_vector, pivot_combination = basis[pivot]
            clean_add(vector, pivot_vector, -factor)
            clean_add(combination, pivot_combination, -factor)
        if not vector:
            dependencies.append(combination)

    print(f"TARGET rows={len(rows)} rank={len(basis)} dependencies={len(dependencies)}")
    for dependency_index, dependency in enumerate(dependencies, start=1):
        print(f"DEPENDENCY {dependency_index}")
        for row_index, coefficient in sorted(dependency.items()):
            print(
                f"  coefficient={signed(coefficient)} "
                f"tree-AHU={metadata[row_index][0]} "
                f"rooted-AHU={metadata[row_index][1]} "
                f"orbit-size={metadata[row_index][2]} "
                f"degrees={metadata[row_index][3]} "
                f"branches={metadata[row_index][4]} "
                f"jdm={metadata[row_index][5]}"
            )


if __name__ == "__main__":
    main()
