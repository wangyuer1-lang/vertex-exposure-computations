from __future__ import annotations

import argparse
import math
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, str(__file__.rsplit("\\", 1)[0]))
from scan_edge_uncover_n7 import uncover_and_lower_dp


RootedTree = tuple["RootedTree", ...]


def rooted_tree_representations(max_n: int) -> dict[int, list[RootedTree]]:
    """Generate canonical rooted-tree tuples as multisets of child tuples."""

    by_size: dict[int, list[RootedTree]] = {1: [()]}
    for n in range(2, max_n + 1):
        pool = [
            (tree, size)
            for size in range(1, n)
            for tree in by_size[size]
        ]
        pool.sort(key=lambda item: (item[1], repr(item[0])))
        generated: set[RootedTree] = set()
        children: list[RootedTree] = []

        def extend(start: int, remaining: int) -> None:
            if remaining == 0:
                generated.add(tuple(children))
                return
            for index in range(start, len(pool)):
                child, size = pool[index]
                if size > remaining:
                    break
                children.append(child)
                extend(index, remaining - size)
                children.pop()

        extend(0, n - 1)
        by_size[n] = sorted(generated, key=repr)
    return by_size


def adjacency_from_rooted(tree: RootedTree) -> list[int]:
    adjacency: list[int] = []

    def add_subtree(node: RootedTree, parent: int | None) -> int:
        vertex = len(adjacency)
        adjacency.append(0)
        if parent is not None:
            adjacency[vertex] |= 1 << parent
            adjacency[parent] |= 1 << vertex
        for child in node:
            add_subtree(child, vertex)
        return vertex

    add_subtree(tree, None)
    return adjacency


def rooted_ahu(adjacency: list[int], root: int, parent: int = -1) -> str:
    child_codes = [
        rooted_ahu(adjacency, child, root)
        for child in range(len(adjacency))
        if adjacency[root] & (1 << child) and child != parent
    ]
    return "(" + "".join(sorted(child_codes)) + ")"


def free_tree_code(adjacency: list[int]) -> str:
    return min(rooted_ahu(adjacency, root) for root in range(len(adjacency)))


def free_trees(max_n: int) -> dict[int, dict[str, list[int]]]:
    rooted = rooted_tree_representations(max_n)
    free: dict[int, dict[str, list[int]]] = {}
    expected = {
        1: 1,
        2: 1,
        3: 1,
        4: 2,
        5: 3,
        6: 6,
        7: 11,
        8: 23,
        9: 47,
        10: 106,
        11: 235,
        12: 551,
        13: 1301,
        14: 3159,
        15: 7741,
        16: 19320,
    }
    for n in range(1, max_n + 1):
        representatives: dict[str, list[int]] = {}
        for tree in rooted[n]:
            adjacency = adjacency_from_rooted(tree)
            representatives.setdefault(free_tree_code(adjacency), adjacency)
        assert len(representatives) == expected[n], (
            n,
            len(representatives),
            expected[n],
        )
        free[n] = representatives
        print(
            f"generated free trees n={n}: {len(representatives)}",
            flush=True,
        )
    return free


def comm_signature(
    law: Counter[tuple[int, ...]]
) -> tuple[tuple[tuple[int, ...], int], ...]:
    commuted: Counter[tuple[int, ...]] = Counter()
    for word, coefficient in law.items():
        commuted[tuple(sorted(word))] += coefficient
    return tuple(sorted(commuted.items()))


def direct_joint_degree_matrix(
    adjacency: list[int],
) -> Counter[tuple[int, int]]:
    degrees = [neighbors.bit_count() for neighbors in adjacency]
    matrix: Counter[tuple[int, int]] = Counter()
    for u, neighbors in enumerate(adjacency):
        remaining = neighbors
        while remaining:
            low_bit = remaining & -remaining
            v = low_bit.bit_length() - 1
            matrix[(degrees[u], degrees[v])] += 1
            remaining ^= low_bit
    return matrix


def joint_degree_matrix_from_uncover(
    law: Counter[tuple[int, ...]],
) -> tuple[Counter[int], Counter[tuple[int, int]]]:
    n = len(next(iter(law)))
    degree_counts: Counter[int] = Counter()
    last_factor = math.factorial(n - 1)
    for word, coefficient in law.items():
        degree_counts[word[-1]] += coefficient
    for degree in list(degree_counts):
        assert degree_counts[degree] % last_factor == 0
        degree_counts[degree] //= last_factor

    q_counts: Counter[tuple[int, int]] = Counter()
    pair_factor = math.factorial(n - 2) if n >= 2 else 1
    if n >= 2:
        for word, coefficient in law.items():
            q_counts[(word[-2], word[-1])] += coefficient
        for key in list(q_counts):
            assert q_counts[key] % pair_factor == 0
            q_counts[key] //= pair_factor

    maximum_degree = max(degree_counts, default=0)
    edge_matrix: Counter[tuple[int, int]] = Counter()
    for second_degree in range(maximum_degree + 1):
        next_value = 0
        for first_degree in range(maximum_degree, -1, -1):
            total_pairs = (
                degree_counts[first_degree] * degree_counts[second_degree]
                - (
                    degree_counts[first_degree]
                    if first_degree == second_degree
                    else 0
                )
            )
            value = (
                total_pairs
                + next_value
                - q_counts[(first_degree, second_degree)]
            )
            assert value >= 0
            if value:
                edge_matrix[(first_degree, second_degree)] = value
            next_value = value
    return degree_counts, edge_matrix


def has_distinct_internal_degrees(adjacency: list[int]) -> bool:
    internal = [
        neighbors.bit_count()
        for neighbors in adjacency
        if neighbors.bit_count() >= 2
    ]
    return len(internal) == len(set(internal))


def is_spider(adjacency: list[int]) -> bool:
    return sum(neighbors.bit_count() >= 3 for neighbors in adjacency) <= 1


def scan(max_n: int) -> None:
    trees = free_trees(max_n)
    total = 0
    distinct_internal_total = 0
    spider_total = 0
    for n in range(1, max_n + 1):
        start = time.perf_counter()
        u_buckets: dict[tuple, list[str]] = defaultdict(list)
        comm_buckets: dict[tuple, list[str]] = defaultdict(list)
        for code, adjacency in trees[n].items():
            law, _ = uncover_and_lower_dp(adjacency)
            u_signature = tuple(sorted(law.items()))
            u_buckets[u_signature].append(code)
            comm_buckets[comm_signature(law)].append(code)

            _, recovered_jdm = joint_degree_matrix_from_uncover(law)
            assert recovered_jdm == direct_joint_degree_matrix(adjacency)

            distinct_internal_total += int(
                has_distinct_internal_degrees(adjacency)
            )
            spider_total += int(is_spider(adjacency))
            total += 1

        u_collisions = [values for values in u_buckets.values() if len(values) > 1]
        comm_collisions = [
            values for values in comm_buckets.values() if len(values) > 1
        ]
        print(
            f"RESULT n={n} trees={len(trees[n])} "
            f"U-collision-buckets={len(u_collisions)} "
            f"comm-collision-buckets={len(comm_collisions)} "
            f"seconds={time.perf_counter() - start:.3f}",
            flush=True,
        )
        for collision in u_collisions[:5]:
            print(f"U-COLLISION {collision}", flush=True)

    print(
        f"SUMMARY trees={total} "
        f"distinct-internal-degree={distinct_internal_total} "
        f"spiders={spider_total}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=10)
    args = parser.parse_args()
    scan(args.max_n)


if __name__ == "__main__":
    main()
