from __future__ import annotations

import argparse
import itertools
from collections import Counter, defaultdict


def edge_pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def relabel_mask(mask: int, pairs: list[tuple[int, int]], pair_index: dict[tuple[int, int], int], permutation: tuple[int, ...]) -> int:
    out = 0
    for bit, (u, v) in enumerate(pairs):
        if mask & (1 << bit):
            a, b = permutation[u], permutation[v]
            if a > b:
                a, b = b, a
            out |= 1 << pair_index[(a, b)]
    return out


def unlabeled_masks(n: int) -> list[int]:
    pairs = edge_pairs(n)
    pair_index = {edge: i for i, edge in enumerate(pairs)}
    permutations = list(itertools.permutations(range(n)))
    representatives: set[int] = set()
    for mask in range(1 << len(pairs)):
        canonical = min(relabel_mask(mask, pairs, pair_index, p) for p in permutations)
        representatives.add(canonical)
    return sorted(representatives)


def adjacency_masks(n: int, edge_mask: int) -> list[int]:
    adj = [0] * n
    for bit, (u, v) in enumerate(edge_pairs(n)):
        if edge_mask & (1 << bit):
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return adj


def uncover(adj: list[int]) -> Counter[tuple[int, ...]]:
    exact: Counter[tuple[int, ...]] = Counter()
    for order in itertools.permutations(range(len(adj))):
        prefix = 0
        word = []
        for v in order:
            word.append((adj[v] & prefix).bit_count())
            prefix |= 1 << v
        exact[tuple(word)] += 1
    return exact


def lower(uncover_law: Counter[tuple[int, ...]]) -> tuple[tuple[tuple[int, ...], int], ...]:
    lowered: Counter[tuple[int, ...]] = Counter()
    for word, coefficient in uncover_law.items():
        for i, value in enumerate(word):
            if value:
                child = list(word)
                child[i] -= 1
                lowered[tuple(child)] += coefficient * value
    return tuple(sorted(lowered.items()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=6)
    args = parser.parse_args()

    for n in range(1, args.max_n + 1):
        representatives = unlabeled_masks(n)
        buckets: dict[tuple[int, tuple], list[int]] = defaultdict(list)
        for mask in representatives:
            law = uncover(adjacency_masks(n, mask))
            buckets[(mask.bit_count(), lower(law))].append(mask)
        collisions = [
            (key, masks)
            for key, masks in buckets.items()
            if len(masks) > 1 and key[0] >= 4
        ]
        print(
            f"n={n} unlabeled={len(representatives)} "
            f"edge-uncover-collision-buckets(m>=4)={len(collisions)}"
        )
        for (m, _), masks in collisions[:10]:
            print(f"  m={m}: {masks}")


if __name__ == "__main__":
    main()
