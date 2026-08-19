from __future__ import annotations

import argparse
from collections import Counter, defaultdict

from tree_u_research import free_trees, rooted_ahu


CutLaw = Counter[tuple[tuple[int, ...], tuple[int, ...]]]


def prefix_suffix_dp(
    adjacency: list[int],
) -> tuple[list[Counter[tuple[int, ...]]], list[Counter[tuple[int, ...]]]]:
    n = len(adjacency)
    full = (1 << n) - 1
    prefix: list[Counter[tuple[int, ...]]] = [Counter() for _ in range(1 << n)]
    suffix: list[Counter[tuple[int, ...]]] = [Counter() for _ in range(1 << n)]
    prefix[0][()] = 1
    for subset in range(1, 1 << n):
        remaining = subset
        while remaining:
            low_bit = remaining & -remaining
            vertex = low_bit.bit_length() - 1
            previous = subset ^ low_bit
            letter = (adjacency[vertex] & previous).bit_count()
            for word, coefficient in prefix[previous].items():
                prefix[subset][word + (letter,)] += coefficient
            remaining ^= low_bit

    suffix[full][()] = 1
    for subset in range(full - 1, -1, -1):
        available = full ^ subset
        while available:
            low_bit = available & -available
            vertex = low_bit.bit_length() - 1
            following = subset | low_bit
            letter = (adjacency[vertex] & subset).bit_count()
            for word, coefficient in suffix[following].items():
                suffix[subset][(letter,) + word] += coefficient
            available ^= low_bit
    return prefix, suffix


def rooted_cut_law(
    root: int,
    prefix: list[Counter[tuple[int, ...]]],
    suffix: list[Counter[tuple[int, ...]]],
) -> CutLaw:
    law: CutLaw = Counter()
    root_bit = 1 << root
    for subset in range(1, len(prefix)):
        if not subset & root_bit:
            continue
        for left, left_coefficient in prefix[subset].items():
            for right, right_coefficient in suffix[subset].items():
                law[(left, right)] += left_coefficient * right_coefficient
    return law


def rooted_orbit_representatives(adjacency: list[int]) -> list[int]:
    representatives: dict[str, int] = {}
    for root in range(len(adjacency)):
        representatives.setdefault(rooted_ahu(adjacency, root), root)
    return list(representatives.values())


def scan(
    max_n: int,
    connected_left_face: bool = False,
    root_first_slice: bool = False,
) -> None:
    trees = free_trees(max_n)
    for n in range(1, max_n + 1):
        buckets: dict[tuple, list[tuple[str, int]]] = defaultdict(list)
        rooted_count = 0
        for tree_code, adjacency in trees[n].items():
            prefix, suffix = prefix_suffix_dp(adjacency)
            for root in rooted_orbit_representatives(adjacency):
                law = rooted_cut_law(root, prefix, suffix)
                if connected_left_face:
                    law = Counter(
                        {
                            pair: coefficient
                            for pair, coefficient in law.items()
                            if sum(pair[0]) == len(pair[0]) - 1
                        }
                    )
                if root_first_slice:
                    law = Counter(
                        {
                            right: coefficient
                            for (left, right), coefficient in law.items()
                            if left == (0,)
                        }
                    )
                signature = tuple(sorted(law.items()))
                buckets[signature].append((tree_code, root))
                rooted_count += 1
        collisions = [bucket for bucket in buckets.values() if len(bucket) > 1]
        print(
            f"ROOTED-CUT n={n} rooted-orbits={rooted_count} "
            f"collision-buckets={len(collisions)}"
        )
        for collision in collisions[:3]:
            print(f"COLLISION {collision}")


def pseudosimilar_target() -> None:
    adjacency = [0] * 11
    edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (8, 10),
    ]
    for u, v in edges:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    prefix, suffix = prefix_suffix_dp(adjacency)
    left = rooted_cut_law(3, prefix, suffix)
    right = rooted_cut_law(7, prefix, suffix)
    differences = sorted(
        (pair, left[pair], right[pair])
        for pair in set(left) | set(right)
        if left[pair] != right[pair]
    )
    print(
        f"PSEUDOSIMILAR-CUT support-left={len(left)} "
        f"support-right={len(right)} differing-pairs={len(differences)}"
    )
    print(f"PSEUDOSIMILAR-CUT witness={differences[0]}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--pseudosimilar-target", action="store_true")
    parser.add_argument("--connected-left-face", action="store_true")
    parser.add_argument("--root-first-slice", action="store_true")
    args = parser.parse_args()
    if args.pseudosimilar_target:
        pseudosimilar_target()
    else:
        scan(
            args.max_n,
            args.connected_left_face,
            args.root_first_slice,
        )


if __name__ == "__main__":
    main()
