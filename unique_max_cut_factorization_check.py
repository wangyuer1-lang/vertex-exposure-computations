from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache

from rooted_cut_law_scan import prefix_suffix_dp, rooted_cut_law
from scan_edge_uncover_n7 import uncover_and_lower_dp
from tree_u_research import free_trees


PairLaw = Counter[tuple[tuple[int, ...], tuple[int, ...]]]


@lru_cache(maxsize=None)
def shuffle_words(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[tuple[tuple[int, ...], int], ...]:
    if not left:
        return ((right, 1),)
    if not right:
        return ((left, 1),)
    result: Counter[tuple[int, ...]] = Counter()
    for suffix, coefficient in shuffle_words(left[1:], right):
        result[(left[0],) + suffix] += coefficient
    for suffix, coefficient in shuffle_words(left, right[1:]):
        result[(right[0],) + suffix] += coefficient
    return tuple(sorted(result.items()))


def tensor_shuffle(left: PairLaw, right: PairLaw) -> PairLaw:
    result: PairLaw = Counter()
    for (left_before, left_after), left_coefficient in left.items():
        for (right_before, right_after), right_coefficient in right.items():
            for before, before_multiplicity in shuffle_words(
                left_before, right_before
            ):
                for after, after_multiplicity in shuffle_words(
                    left_after, right_after
                ):
                    result[(before, after)] += (
                        left_coefficient
                        * right_coefficient
                        * before_multiplicity
                        * after_multiplicity
                    )
    return result


def components_after_root(
    adjacency: list[int], root: int
) -> list[tuple[list[int], int]]:
    n = len(adjacency)
    unseen = set(range(n)) - {root}
    components: list[tuple[list[int], int]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        vertices: list[int] = []
        unseen.remove(start)
        while stack:
            u = stack.pop()
            vertices.append(u)
            for v in list(unseen):
                if adjacency[u] & (1 << v):
                    unseen.remove(v)
                    stack.append(v)
        old_to_new = {old: new for new, old in enumerate(vertices)}
        local = [0] * len(vertices)
        attachment = -1
        for old_u in vertices:
            new_u = old_to_new[old_u]
            if adjacency[root] & (1 << old_u):
                attachment = new_u
            for old_v in vertices:
                if adjacency[old_u] & (1 << old_v):
                    local[new_u] |= 1 << old_to_new[old_v]
        assert attachment >= 0
        components.append((local, attachment))
    return components


def extracted_max_slice(adjacency: list[int], maximum: int) -> PairLaw:
    law, _ = uncover_and_lower_dp(adjacency)
    result: PairLaw = Counter()
    for word, coefficient in law.items():
        positions = [index for index, letter in enumerate(word) if letter == maximum]
        if not positions:
            continue
        assert len(positions) == 1
        position = positions[0]
        result[(word[:position], word[position + 1 :])] += coefficient
    return result


def branch_product(adjacency: list[int], root: int) -> PairLaw:
    product: PairLaw = Counter({((), ()): 1})
    for branch, attachment in components_after_root(adjacency, root):
        prefix, suffix = prefix_suffix_dp(branch)
        factor = rooted_cut_law(attachment, prefix, suffix)
        product = tensor_shuffle(product, factor)
    return product


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=9)
    args = parser.parse_args()
    trees = free_trees(args.max_n)
    checked = 0
    for n in range(1, args.max_n + 1):
        count = 0
        for adjacency in trees[n].values():
            degrees = [row.bit_count() for row in adjacency]
            maximum = max(degrees)
            roots = [v for v, degree in enumerate(degrees) if degree == maximum]
            if len(roots) != 1:
                continue
            root = roots[0]
            actual = extracted_max_slice(adjacency, maximum)
            expected = branch_product(adjacency, root)
            assert actual == expected
            count += 1
            checked += 1
        print(f"UNIQUE-MAX n={n} trees={count} factorization-passed={count}")
    print(f"SUMMARY unique-max-trees-checked={checked} failures=0")


if __name__ == "__main__":
    main()
