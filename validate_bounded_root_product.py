from __future__ import annotations

import argparse

from rooted_cut_law_scan import rooted_orbit_representatives
from rooted_cut_recursion_check import child_product
from root_product_rank_scan import bounded_child_product
from tree_u_research import free_trees


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=10)
    parser.add_argument("--max-excess", type=int, default=3)
    args = parser.parse_args()

    trees = free_trees(args.max_n)
    checks = 0
    for n in range(1, args.max_n + 1):
        level_checks = 0
        for adjacency in trees[n].values():
            degrees = [row.bit_count() for row in adjacency]
            maximum = max(degrees)
            for root in rooted_orbit_representatives(adjacency):
                degree = degrees[root]
                if degree != maximum:
                    continue
                full = child_product(adjacency, root)
                for excess in range(args.max_excess + 1):
                    bound = degree + excess
                    expected = type(full)(
                        {
                            pair: coefficient
                            for pair, coefficient in full.items()
                            if len(pair[0]) <= bound
                        }
                    )
                    actual = bounded_child_product(adjacency, root, bound)
                    assert actual == expected
                    connected_expected = type(full)(
                        {
                            pair: coefficient
                            for pair, coefficient in expected.items()
                            if sum(pair[0]) == len(pair[0]) - degree
                        }
                    )
                    connected_actual = bounded_child_product(
                        adjacency,
                        root,
                        bound,
                        connected_left_face=True,
                    )
                    assert connected_actual == connected_expected
                    checks += 2
                    level_checks += 2
        print(
            f"BOUNDED-PRODUCT n={n} checks={level_checks} failures=0",
            flush=True,
        )
    print(f"SUMMARY checks={checks} failures=0")


if __name__ == "__main__":
    main()
