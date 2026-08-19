from __future__ import annotations

from collections import defaultdict

from tree_u_research import free_tree_code, free_trees, rooted_ahu


def delete_vertex(adjacency: list[int], deleted: int) -> list[list[int]]:
    n = len(adjacency)
    remaining = [v for v in range(n) if v != deleted]
    old_to_new = {old: new for new, old in enumerate(remaining)}
    reduced = [0] * (n - 1)
    for old_u in remaining:
        new_u = old_to_new[old_u]
        for old_v in remaining:
            if adjacency[old_u] & (1 << old_v):
                reduced[new_u] |= 1 << old_to_new[old_v]

    components: list[list[int]] = []
    unseen = set(range(n - 1))
    while unseen:
        start = min(unseen)
        stack = [start]
        vertices: list[int] = []
        unseen.remove(start)
        while stack:
            u = stack.pop()
            vertices.append(u)
            for v in list(unseen):
                if reduced[u] & (1 << v):
                    unseen.remove(v)
                    stack.append(v)
        local_index = {old: new for new, old in enumerate(vertices)}
        component = [0] * len(vertices)
        for old_u in vertices:
            for old_v in vertices:
                if reduced[old_u] & (1 << old_v):
                    component[local_index[old_u]] |= 1 << local_index[old_v]
        components.append(component)
    return components


def forest_code(components: list[list[int]]) -> tuple[str, ...]:
    return tuple(sorted(free_tree_code(component) for component in components))


def edge_list(adjacency: list[int]) -> list[tuple[int, int]]:
    return [
        (u, v)
        for u in range(len(adjacency))
        for v in range(u + 1, len(adjacency))
        if adjacency[u] & (1 << v)
    ]


def main() -> None:
    trees = free_trees(11)
    first = False
    for n in range(1, 12):
        examples = []
        pair_count = 0
        for tree_code, adjacency in trees[n].items():
            by_deleted_forest: dict[tuple[str, ...], list[int]] = defaultdict(list)
            for root in range(n):
                by_deleted_forest[forest_code(delete_vertex(adjacency, root))].append(root)
            for vertices in by_deleted_forest.values():
                rooted_classes: dict[str, list[int]] = defaultdict(list)
                for root in vertices:
                    rooted_classes[rooted_ahu(adjacency, root)].append(root)
                if len(rooted_classes) > 1:
                    classes = list(rooted_classes.values())
                    pair_count += sum(
                        len(classes[i]) * len(classes[j])
                        for i in range(len(classes))
                        for j in range(i + 1, len(classes))
                    )
                    examples.append((tree_code, adjacency, classes))
        print(
            f"PSEUDOSIMILAR n={n} trees-with-examples={len(examples)} "
            f"cross-orbit-pairs={pair_count}"
        )
        if examples and not first:
            first = True
            tree_code, adjacency, classes = examples[0]
            print(f"FIRST tree-AHU={tree_code}")
            print(f"FIRST edges={edge_list(adjacency)}")
            print(f"FIRST rooted-orbit-classes={classes}")
            for root_class in classes:
                root = root_class[0]
                print(
                    f"ROOT v={root} degree={adjacency[root].bit_count()} "
                    f"rooted-AHU={rooted_ahu(adjacency, root)} "
                    f"deleted={forest_code(delete_vertex(adjacency, root))}"
                )


if __name__ == "__main__":
    main()
