from __future__ import annotations

import argparse
from collections import defaultdict

from root_product_rank_scan import PRIME, bounded_child_product
from rooted_cut_law_scan import rooted_orbit_representatives
from tree_u_research import direct_joint_degree_matrix, free_trees, rooted_ahu


def clean_add(target: dict, source: dict, factor: int) -> None:
    for key, value in source.items():
        updated = (target.get(key, 0) + factor * value) % PRIME
        if updated:
            target[key] = updated
        else:
            target.pop(key, None)


def dependencies(rows: list[dict]) -> tuple[int, list[dict[int, int]]]:
    basis = {}
    result = []
    for index, source in enumerate(rows):
        vector = {
            key: value % PRIME
            for key, value in source.items()
            if value % PRIME
        }
        combination = {index: 1}
        while vector:
            pivot = max(vector)
            if pivot not in basis:
                inverse = pow(vector[pivot], PRIME - 2, PRIME)
                vector = {
                    key: value * inverse % PRIME
                    for key, value in vector.items()
                }
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
            result.append(combination)
    return len(basis), result


def signed(value: int) -> int:
    value %= PRIME
    return value if value <= PRIME // 2 else value - PRIME


def normalize(relation: dict[int, int]) -> dict[int, int]:
    first = relation[min(relation)]
    inverse = pow(first, PRIME - 2, PRIME)
    return {index: coefficient * inverse % PRIME for index, coefficient in relation.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=11)
    parser.add_argument("--base-excess", type=int, default=2)
    parser.add_argument("--witness-excess", type=int, default=3)
    args = parser.parse_args()

    groups = defaultdict(list)
    for tree_code, adjacency in free_trees(args.n)[args.n].items():
        degrees = tuple(row.bit_count() for row in adjacency)
        maximum = max(degrees)
        degree_sequence = tuple(sorted(degrees, reverse=True))
        jdm = tuple(sorted(direct_joint_degree_matrix(adjacency).items()))
        key = (maximum, degree_sequence, jdm)
        for root in rooted_orbit_representatives(adjacency):
            if degrees[root] != maximum:
                continue
            rooted_code = rooted_ahu(adjacency, root)
            orbit_size = sum(
                rooted_ahu(adjacency, vertex) == rooted_code
                for vertex in range(args.n)
            )
            base = dict(
                bounded_child_product(
                    adjacency, root, maximum + args.base_excess
                )
            )
            witness = dict(
                bounded_child_product(
                    adjacency, root, maximum + args.witness_excess
                )
            )
            groups[key].append(
                (base, witness, tree_code, rooted_code, orbit_size)
            )

    bad_group_count = 0
    for key, records in sorted(groups.items()):
        rank, relations = dependencies([record[0] for record in records])
        if not relations:
            continue
        bad_group_count += 1
        print(
            f"BAD-GROUP key={key!r} rank={rank}/{len(records)} "
            f"dependencies={len(relations)}"
        )
        for relation_number, raw_relation in enumerate(relations, start=1):
            relation = normalize(raw_relation)
            print(f"  RELATION {relation_number}")
            for index, coefficient in sorted(relation.items()):
                record = records[index]
                print(
                    f"    row={index} coefficient={signed(coefficient)} "
                    f"tree-AHU={record[2]} rooted-AHU={record[3]} "
                    f"orbit-size={record[4]}"
                )

            residual = {}
            for index, coefficient in relation.items():
                clean_add(residual, records[index][1], coefficient)
            assert residual
            witness_pair = min(
                residual,
                key=lambda pair: (len(pair[0]), pair[0], pair[1]),
            )
            print(
                f"    FIRST-WITNESS left-length={len(witness_pair[0])} "
                f"pair={witness_pair!r} residual={signed(residual[witness_pair])}"
            )
            for index, coefficient in sorted(relation.items()):
                entry = records[index][1].get(witness_pair, 0)
                if entry:
                    print(
                        f"      row={index} relation-coefficient={signed(coefficient)} "
                        f"entry={entry}"
                    )
    print(f"SUMMARY bad-groups={bad_group_count}")


if __name__ == "__main__":
    main()
