from __future__ import annotations

import argparse
from collections import Counter, defaultdict

from scan_edge_uncover_n7 import uncover_and_lower_dp
from tree_u_research import free_trees


def maximum_root_sum(
    law: Counter[tuple[int, ...]], maximum: int
) -> Counter[tuple[tuple[int, ...], tuple[int, ...]]]:
    result: Counter[tuple[tuple[int, ...], tuple[int, ...]]] = Counter()
    for word, coefficient in law.items():
        for position, letter in enumerate(word):
            if letter == maximum:
                result[(word[:position], word[position + 1 :])] += coefficient
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=10)
    args = parser.parse_args()
    trees = free_trees(args.max_n)
    for n in range(1, args.max_n + 1):
        bare_buckets = defaultdict(list)
        enriched_buckets = defaultdict(list)
        for tree_code, adjacency in trees[n].items():
            degrees = [row.bit_count() for row in adjacency]
            degree_counts = tuple(sorted(Counter(degrees).items()))
            maximum = max(degrees)
            law, _ = uncover_and_lower_dp(adjacency)
            root_sum = tuple(sorted(maximum_root_sum(law, maximum).items()))
            bare_buckets[(maximum, root_sum)].append(tree_code)
            enriched_buckets[(degree_counts, root_sum)].append(tree_code)
        bare = [bucket for bucket in bare_buckets.values() if len(bucket) > 1]
        enriched = [
            bucket for bucket in enriched_buckets.values() if len(bucket) > 1
        ]
        print(
            f"MAX-ROOT-SUM n={n} trees={len(trees[n])} "
            f"bare-collision-buckets={len(bare)} "
            f"degree-enriched-collision-buckets={len(enriched)}"
        )
        for collision in enriched[:3]:
            print(f"COLLISION {collision}")


if __name__ == "__main__":
    main()
