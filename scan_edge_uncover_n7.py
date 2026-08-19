from __future__ import annotations

import argparse
import itertools
import math
import time
from collections import Counter, defaultdict


def edge_pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def adjacency_from_mask(n: int, mask: int) -> list[int]:
    adjacency = [0] * n
    for bit, (u, v) in enumerate(edge_pairs(n)):
        if mask & (1 << bit):
            adjacency[u] |= 1 << v
            adjacency[v] |= 1 << u
    return adjacency


def stable_colors(adjacency: list[int]) -> list[int]:
    n = len(adjacency)
    degrees = [row.bit_count() for row in adjacency]
    degree_values = {value: i for i, value in enumerate(sorted(set(degrees)))}
    colors = [degree_values[value] for value in degrees]

    while True:
        number_of_colors = max(colors) + 1
        signatures = []
        for vertex in range(n):
            neighbor_color_counts = [0] * number_of_colors
            neighbors = adjacency[vertex]
            while neighbors:
                low_bit = neighbors & -neighbors
                neighbor = low_bit.bit_length() - 1
                neighbor_color_counts[colors[neighbor]] += 1
                neighbors ^= low_bit
            signatures.append((colors[vertex], tuple(neighbor_color_counts)))
        signature_values = {
            value: i for i, value in enumerate(sorted(set(signatures)))
        }
        refined = [signature_values[value] for value in signatures]
        if refined == colors:
            return colors
        colors = refined


def canonical_mask(adjacency: list[int]) -> int:
    """Exact canonical form, using stable color refinement only as safe pruning."""

    n = len(adjacency)
    colors = stable_colors(adjacency)
    color_classes = [
        tuple(vertex for vertex, color in enumerate(colors) if color == class_color)
        for class_color in range(max(colors) + 1)
    ]
    pairs = edge_pairs(n)
    best: int | None = None

    for choices in itertools.product(
        *(itertools.permutations(group) for group in color_classes)
    ):
        # New labels are assigned class-by-class; order[new_label] = old_label.
        order = tuple(vertex for choice in choices for vertex in choice)
        candidate = 0
        for bit, (new_u, new_v) in enumerate(pairs):
            old_u, old_v = order[new_u], order[new_v]
            if adjacency[old_u] & (1 << old_v):
                candidate |= 1 << bit
        if best is None or candidate < best:
            best = candidate

    assert best is not None
    return best


def unlabeled_masks_by_augmentation(max_n: int) -> dict[int, list[int]]:
    """Generate every unlabeled graph by adding one vertex, then canonicalize."""

    representatives: dict[int, list[int]] = {1: [0]}
    for n in range(2, max_n + 1):
        start = time.perf_counter()
        current: set[int] = set()
        for previous_mask in representatives[n - 1]:
            previous_adjacency = adjacency_from_mask(n - 1, previous_mask)
            for neighborhood in range(1 << (n - 1)):
                adjacency = previous_adjacency[:] + [0]
                new_vertex = n - 1
                for vertex in range(n - 1):
                    if neighborhood & (1 << vertex):
                        adjacency[vertex] |= 1 << new_vertex
                        adjacency[new_vertex] |= 1 << vertex
                current.add(canonical_mask(adjacency))
        representatives[n] = sorted(current)
        elapsed = time.perf_counter() - start
        print(
            f"generated n={n}: unlabeled={len(current)} "
            f"expected-known={'yes' if len(current) in {2, 4, 11, 34, 156, 1044, 12346} else 'check'} "
            f"seconds={elapsed:.3f}",
            flush=True,
        )
    return representatives


def uncover(
    adjacency: list[int], orders: list[tuple[int, ...]]
) -> Counter[tuple[int, ...]]:
    exact: Counter[tuple[int, ...]] = Counter()
    for order in orders:
        prefix = 0
        word = []
        for vertex in order:
            word.append((adjacency[vertex] & prefix).bit_count())
            prefix |= 1 << vertex
        exact[tuple(word)] += 1
    return exact


def lower_signature(
    uncover_law: Counter[tuple[int, ...]]
) -> tuple[tuple[tuple[int, ...], int], ...]:
    lowered: Counter[tuple[int, ...]] = Counter()
    for word, coefficient in uncover_law.items():
        for position, value in enumerate(word):
            if value:
                child = list(word)
                child[position] -= 1
                lowered[tuple(child)] += coefficient * value
    return tuple(sorted(lowered.items()))


def uncover_and_lower_dp(
    adjacency: list[int],
) -> tuple[
    Counter[tuple[int, ...]],
    Counter[tuple[int, ...]],
]:
    """Compute U and ∂U together from the last-vertex recurrence on subsets."""

    n = len(adjacency)
    laws: list[Counter[tuple[int, ...]]] = [Counter() for _ in range(1 << n)]
    lowers: list[Counter[tuple[int, ...]]] = [Counter() for _ in range(1 << n)]
    laws[0][()] = 1

    for subset in range(1, 1 << n):
        remaining = subset
        while remaining:
            vertex_bit = remaining & -remaining
            vertex = vertex_bit.bit_length() - 1
            previous = subset ^ vertex_bit
            degree = (adjacency[vertex] & previous).bit_count()

            for word, coefficient in laws[previous].items():
                laws[subset][word + (degree,)] += coefficient
                if degree:
                    lowers[subset][word + (degree - 1,)] += coefficient * degree
            for word, coefficient in lowers[previous].items():
                lowers[subset][word + (degree,)] += coefficient

            remaining ^= vertex_bit

    return laws[-1], lowers[-1]


def recover_uncover_no_isolates(
    lowered: Counter[tuple[int, ...]], n: int
) -> Counter[tuple[int, ...]]:
    """Invert ∂ on a graph-derived law when the graph has no isolated vertices."""

    suffix_slices: list[Counter[tuple[int, ...]]] = [
        Counter() for _ in range(n)
    ]
    for word, coefficient in lowered.items():
        suffix_slices[word[-1]][word[:-1]] += coefficient

    degree_card_sums: list[Counter[tuple[int, ...]]] = [
        Counter() for _ in range(n)
    ]
    # B_0 = 0 because the graph has no isolated vertices.  The identity
    # S_k = ∂B_k + (k+1)B_(k+1) then gives every B_(k+1).
    for degree in range(n - 1):
        derivative = Counter(dict(lower_signature(degree_card_sums[degree])))
        numerator = suffix_slices[degree].copy()
        numerator.subtract(derivative)
        numerator = Counter(
            {word: value for word, value in numerator.items() if value}
        )
        divisor = degree + 1
        next_sum: Counter[tuple[int, ...]] = Counter()
        for word, value in numerator.items():
            assert value % divisor == 0
            quotient = value // divisor
            assert quotient >= 0
            if quotient:
                next_sum[word] = quotient
        degree_card_sums[degree + 1] = next_sum

    recovered: Counter[tuple[int, ...]] = Counter()
    for degree, card_sum in enumerate(degree_card_sums):
        for word, coefficient in card_sum.items():
            recovered[word + (degree,)] += coefficient
    return recovered


def graph6(n: int, mask: int) -> str:
    """Graph6 for n <= 62, useful as a reproducible collision certificate."""

    bits = []
    # graph6 uses upper triangle by columns: (0,1),(0,2),(1,2),...
    pair_to_bit = {pair: i for i, pair in enumerate(edge_pairs(n))}
    for column in range(1, n):
        for row in range(column):
            bits.append((mask >> pair_to_bit[(row, column)]) & 1)
    while len(bits) % 6:
        bits.append(0)
    payload = []
    for start in range(0, len(bits), 6):
        value = 0
        for bit in bits[start : start + 6]:
            value = (value << 1) | bit
        payload.append(chr(value + 63))
    return chr(n + 63) + "".join(payload)


def scan(
    representatives: dict[int, list[int]], max_n: int, use_dp: bool
) -> None:
    for n in range(1, max_n + 1):
        start = time.perf_counter()
        orders = [] if use_dp else list(itertools.permutations(range(n)))
        buckets: dict[tuple[int, tuple], list[int]] = defaultdict(list)
        u_buckets: dict[tuple[int, tuple], list[int]] = defaultdict(list)
        for index, mask in enumerate(representatives[n], start=1):
            adjacency = adjacency_from_mask(n, mask)
            if use_dp:
                law, lowered = uncover_and_lower_dp(adjacency)
                lowered_signature = tuple(sorted(lowered.items()))
            else:
                law = uncover(adjacency, orders)
                lowered_signature = lower_signature(law)
            law_signature = tuple(sorted(law.items()))
            buckets[(mask.bit_count(), lowered_signature)].append(mask)
            u_buckets[(mask.bit_count(), law_signature)].append(mask)
            if index % 100 == 0:
                print(
                    f"scan n={n}: {index}/{len(representatives[n])}",
                    flush=True,
                )

        collisions = [
            (key, masks)
            for key, masks in buckets.items()
            if len(masks) > 1 and key[0] >= 4
        ]
        u_collisions = [
            (key, masks)
            for key, masks in u_buckets.items()
            if len(masks) > 1
        ]
        elapsed = time.perf_counter() - start
        print(
            f"RESULT n={n} unlabeled={len(representatives[n])} "
            f"permutations={math.factorial(n)} "
            f"method={'subset-dp' if use_dp else 'permutations'} "
            f"U-collision-buckets={len(u_collisions)} "
            f"edge-uncover-collision-buckets(m>=4)={len(collisions)} "
            f"seconds={elapsed:.3f}",
            flush=True,
        )
        for (edge_count, _), masks in collisions[:20]:
            certificates = [graph6(n, mask) for mask in masks]
            print(
                f"COLLISION n={n} m={edge_count} graph6={certificates}",
                flush=True,
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=7)
    parser.add_argument(
        "--permutations",
        action="store_true",
        help="Use direct permutation enumeration instead of subset DP.",
    )
    args = parser.parse_args()

    representatives = unlabeled_masks_by_augmentation(args.max_n)
    scan(representatives, args.max_n, use_dp=not args.permutations)


if __name__ == "__main__":
    main()
