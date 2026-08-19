from __future__ import annotations

import argparse
from collections import defaultdict

from rooted_cut_law_scan import rooted_orbit_representatives
from rooted_cut_recursion_check import child_product
from tree_u_research import direct_joint_degree_matrix, free_trees, rooted_ahu


def extrema(law):
    pairs = list(law)
    keys = {
        "raw-min": min(pairs),
        "raw-max": max(pairs),
        "leftlen-min": min(pairs, key=lambda p: (len(p[0]), p[0], p[1])),
        "leftlen-max": max(pairs, key=lambda p: (len(p[0]), p[0], p[1])),
        "rightlen-min": min(pairs, key=lambda p: (len(p[1]), p[1], p[0])),
        "rightlen-max": max(pairs, key=lambda p: (len(p[1]), p[1], p[0])),
        "flat-lowsep-min": min(pairs, key=lambda p: p[0] + (-1,) + p[1]),
        "flat-lowsep-max": max(pairs, key=lambda p: p[0] + (-1,) + p[1]),
        "flat-highsep-min": min(pairs, key=lambda p: p[0] + (99,) + p[1]),
        "flat-highsep-max": max(pairs, key=lambda p: p[0] + (99,) + p[1]),
    }
    return {
        name: (pair, law[pair])
        for name, pair in keys.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=9)
    parser.add_argument("--maximum-roots-only", action="store_true")
    parser.add_argument("--group-jdm", action="store_true")
    args = parser.parse_args()
    trees = free_trees(args.max_n)
    for n in range(1, args.max_n + 1):
        records = []
        for tree_code, adjacency in trees[n].items():
            maximum = max(row.bit_count() for row in adjacency)
            degree_sequence = tuple(
                sorted((row.bit_count() for row in adjacency), reverse=True)
            )
            jdm = tuple(sorted(direct_joint_degree_matrix(adjacency).items()))
            for root in rooted_orbit_representatives(adjacency):
                degree = adjacency[root].bit_count()
                if args.maximum_roots_only and degree != maximum:
                    continue
                law = child_product(adjacency, root)
                records.append(
                    (
                        degree,
                        degree_sequence,
                        jdm,
                        tree_code,
                        rooted_ahu(adjacency, root),
                        extrema(law),
                    )
                )
        names = list(records[0][5])
        summary = []
        for name in names:
            buckets = defaultdict(list)
            for degree, degree_sequence, jdm, tree_code, rooted_code, values in records:
                invariant = (degree, degree_sequence, jdm) if args.group_jdm else degree
                buckets[(invariant, values[name])].append((tree_code, rooted_code))
            collisions = [bucket for bucket in buckets.values() if len(bucket) > 1]
            summary.append(
                f"{name}={len(collisions)}/{sum(map(len, collisions))}"
            )
        print(f"EXTREMA n={n} rooted={len(records)} " + " ".join(summary))


if __name__ == "__main__":
    main()
