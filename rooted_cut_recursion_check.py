from __future__ import annotations

import argparse
from collections import Counter

from rooted_cut_law_scan import (
    prefix_suffix_dp,
    rooted_cut_law,
    rooted_orbit_representatives,
)
from tree_u_research import free_trees
from unique_max_cut_factorization_check import (
    PairLaw,
    components_after_root,
    tensor_shuffle,
)


def extract_full_root_slice(law: PairLaw, vertex_count: int) -> PairLaw:
    cumulative: list[Counter[tuple[int, ...]]] = [
        Counter() for _ in range(vertex_count + 1)
    ]
    for (left, right), coefficient in law.items():
        cumulative[len(left)][left + right] += coefficient

    exact: list[Counter[tuple[int, ...]]] = [
        Counter() for _ in range(vertex_count + 1)
    ]
    for position in range(1, vertex_count + 1):
        exact[position] = cumulative[position] - cumulative[position - 1]

    last_words = exact[vertex_count]
    root_degrees = {word[-1] for word in last_words}
    assert len(root_degrees) == 1
    root_degree = next(iter(root_degrees))

    result: PairLaw = Counter()
    for position in range(1, vertex_count + 1):
        for word, coefficient in exact[position].items():
            if word[position - 1] != root_degree:
                continue
            result[
                (word[: position - 1], word[position:])
            ] += coefficient
    return result


def child_product(adjacency: list[int], root: int) -> PairLaw:
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
            prefix, suffix = prefix_suffix_dp(adjacency)
            for root in rooted_orbit_representatives(adjacency):
                law = rooted_cut_law(root, prefix, suffix)
                actual = extract_full_root_slice(law, n)
                expected = child_product(adjacency, root)
                assert actual == expected
                count += 1
                checked += 1
        print(f"ROOTED-RECURSION n={n} rooted-orbits={count} passed={count}")
    print(f"SUMMARY rooted-orbits-checked={checked} failures=0")


if __name__ == "__main__":
    main()
