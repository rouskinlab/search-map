from itertools import combinations

items = ['AS1', 'PS1', 'PS2', 'PS3', 'LS2b', 'LS3', 'LS4']
exclusive = {("AS1", "PS1"), ("PS1", "LS4"), ("LS4", "PS2"), ("LS3", "PS3"), ("LS2b", "PS2"), ("PS2", "LS3")}


def is_valid(result):
    for item1, item2 in combinations(result, 2):
        if (item1, item2) in exclusive or (item2, item1) in exclusive:
            # This combination of items is forbidden.
            return False
    return True


results = list()
for size in range(1, len(items)):
    for result in combinations(items, size):
        if is_valid(result):
            results.append(set(result))


def is_best(result: set[str], results: list[set[str]]):
    for other in results:
        if other - result and not result - other:
            # result is a proper subset of another result, so not the best
            return False
    return True

best_results = [sorted(result) for result in results if is_best(result, results)]
for result in best_results:
    print(" + ".join(result))
