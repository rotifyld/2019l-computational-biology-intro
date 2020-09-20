from typing import Dict, List, Collection, Tuple

V = str
Tree = Dict[V, List[V]]

g: Tree = {"abcde": ["abcd", "e"], "abcd": ["ab", "cd"], "e": [], "ab": ["a", "b"], "cd": ["c", "d"], "a": [], "b": [],
           "c": [], "d": []}
s: Tree = {"abcde": ["abcd", "e"], "abcd": ["abd", "c"], "abd": ["a", "bd"], "bd": ["b", "d"], "a": [], "b": [],
           "c": [], "d": [], "e": []}


def v_le(v: V, u: V) -> bool:
    return {c for c in v}.issubset({c for c in u})


def v_min(vs: Collection[V]) -> V:
    if len(vs) == 0:
        raise RuntimeError("Minimal element of empty collection")

    ret = next(iter(vs))
    for v in vs:
        if v_le(v, ret):
            ret = v
    return ret


def t_root(t: Tree) -> V:
    return max(((len(v), v) for v in t.keys()))[1]


def t_parent(t: Tree, v: V) -> V:
    if t_root(t) == v:
        raise RuntimeError("Parent of a root")

    for parent, children in t.items():
        for child in children:
            if child == v:
                return parent

    raise RuntimeError("No parent found")


def t_distance(t: Tree, v: V, ancestor: V) -> int:
    dist = 0
    while v != ancestor:
        dist += 1
        v = t_parent(t, v)

    return dist


def lca(g: Tree, s: Tree) -> Dict[V, V]:
    return {v: v_min(list(filter(lambda u: v_le(v, u), (u for u in s.keys())))) for v in g.keys()}


def dc_d(g: Tree, s: Tree) -> Tuple[int, int]:
    mapping = lca(g, s)
    total_dc = 0
    total_d = 0

    for v in g.keys():
        if v != t_root(g):
            dist = t_distance(s, mapping[v], mapping[t_parent(g, v)])
            if dist == 0:
                total_d += 1
            total_dc += dist

    return total_dc, total_d


def dc(g: Tree, s: Tree) -> int:
    return dc_d(g, s)[0]


def d(g: Tree, s: Tree) -> int:
    return dc_d(g, s)[1]


def l(g: Tree, s: Tree) -> int:
    dc, d = dc_d(g, s)
    tree_size = len(g)
    return dc + 2 * d - tree_size + 1


def ex():
    print(lca(g, s))
    print(dc_d(g, s))
    print(l(g, s))


def test_tree():
    assert t_root(g) == "abcde"
    assert t_root(s) == "abcde"
    assert t_parent(g, "a") == "ab"
    assert t_parent(s, "a") == "abd"
    assert v_le("a", "ab")
    assert v_le("a", "qvdfirfeulirhflkejrhba")
    assert v_min(["a", "as", "asdf"]) == "a"
    assert v_min(["z", "az", "abcdz"]) == "z"
    print("pronto")


if __name__ == '__main__':
    ex()
    # test_tree()
