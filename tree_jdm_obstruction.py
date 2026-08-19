from __future__ import annotations

from collections import defaultdict

from tree_u_research import direct_joint_degree_matrix, free_trees


def edge_list(adjacency: list[int]) -> list[tuple[int, int]]:
    return [
        (u, v)
        for u in range(len(adjacency))
        for v in range(u + 1, len(adjacency))
        if adjacency[u] & (1 << v)
    ]


def main() -> None:
    trees = free_trees(11)
    first_reported = False
    for n in range(1, 12):
        buckets: dict[tuple, list[tuple[str, list[int]]]] = defaultdict(list)
        for code, adjacency in trees[n].items():
            signature = tuple(sorted(direct_joint_degree_matrix(adjacency).items()))
            buckets[signature].append((code, adjacency))

        collisions = [bucket for bucket in buckets.values() if len(bucket) > 1]
        print(
            f"JDM n={n} collision-buckets={len(collisions)} "
            f"trees-in-buckets={sum(map(len, collisions))} "
            f"largest-bucket={max(map(len, collisions), default=0)}"
        )
        if collisions and not first_reported:
            first_reported = True
            print("FIRST-JDM-COLLISION")
            for index, (code, adjacency) in enumerate(collisions[0], start=1):
                print(f"tree-{index} AHU={code}")
                print(
                    f"tree-{index} degrees="
                    f"{sorted((row.bit_count() for row in adjacency), reverse=True)}"
                )
                print(f"tree-{index} edges={edge_list(adjacency)}")


if __name__ == "__main__":
    main()
