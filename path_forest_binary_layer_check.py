from __future__ import annotations

from collections import Counter
from functools import lru_cache

from scan_edge_uncover_n7 import uncover_and_lower_dp


def partitions(total: int, minimum: int = 1) -> list[tuple[int, ...]]:
    if total == 0:
        return [()]
    result: list[tuple[int, ...]] = []
    for first in range(minimum, total + 1):
        for rest in partitions(total - first, first):
            result.append((first,) + rest)
    return result


def path_forest(lengths: tuple[int, ...]) -> list[int]:
    adjacency = [0] * sum(lengths)
    offset = 0
    for length in lengths:
        for local in range(length - 1):
            u = offset + local
            v = u + 1
            adjacency[u] |= 1 << v
            adjacency[v] |= 1 << u
        offset += length
    return adjacency


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


def shuffle_law(
    left: Counter[tuple[int, ...]], right: Counter[tuple[int, ...]]
) -> Counter[tuple[int, ...]]:
    result: Counter[tuple[int, ...]] = Counter()
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            for word, multiplicity in shuffle_words(left_word, right_word):
                result[word] += (
                    left_coefficient * right_coefficient * multiplicity
                )
    return result


def expected_layer(lengths: tuple[int, ...]) -> Counter[tuple[int, ...]]:
    result: Counter[tuple[int, ...]] = Counter({(): 1})
    for length in lengths:
        factor = Counter({(0,) + (1,) * (length - 1): 2 ** (length - 1)})
        result = shuffle_law(result, factor)
    return result


def main() -> None:
    checked = 0
    for n in range(1, 12):
        count = 0
        for lengths in partitions(n):
            component_count = len(lengths)
            law, _ = uncover_and_lower_dp(path_forest(lengths))
            actual = Counter(
                {
                    word: coefficient
                    for word, coefficient in law.items()
                    if word.count(0) == component_count
                    and all(letter <= 1 for letter in word)
                }
            )
            assert actual == expected_layer(lengths), lengths
            count += 1
            checked += 1
        print(f"PATH-FOREST n={n} integer-partitions={count} passed={count}")
    print(f"SUMMARY path-forests-checked={checked} failures=0")


if __name__ == "__main__":
    main()
