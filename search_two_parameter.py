import itertools

from two_parameter_lp import evaluate_two_parameter_bound, make_direction_family, standard_family


def plus_minus_family(m, shifts_u, shifts_v, shifts_uv, shifts_umv):
    return (
        make_direction_family((1, 0), shifts_u, prefix="u")
        + make_direction_family((0, 1), shifts_v, prefix="v")
        + make_direction_family((1, 1), shifts_uv, dilation=m, prefix="uv")
        + make_direction_family((1, -1), shifts_umv, dilation=m, prefix="umv")
    )


def search_standard_family(m, max_shift=2, max_size=2):
    best_bound = 1.0
    best_case = None

    options = []
    for size in range(1, max_size + 1):
        options.extend(itertools.combinations(range(max_shift + 1), size))

    for shifts_u, shifts_v, shifts_uv in itertools.product(options, repeat=3):
        forms = standard_family(m, shifts_u, shifts_v, shifts_uv)
        bound = evaluate_two_parameter_bound(m, forms)
        if bound < best_bound - 1e-9:
            best_bound = bound
            best_case = (shifts_u, shifts_v, shifts_uv)
            print(f"standard m={m} new best {best_bound:.6f} with {best_case}", flush=True)

    return best_bound, best_case


def search_plus_minus_family(m):
    best_bound = 1.0
    best_case = None
    options = [(0,), (1,), (0, 1)]

    for shifts_u, shifts_v, shifts_uv, shifts_umv in itertools.product(options, repeat=4):
        forms = plus_minus_family(m, shifts_u, shifts_v, shifts_uv, shifts_umv)
        bound = evaluate_two_parameter_bound(m, forms)
        if bound < best_bound - 1e-9:
            best_bound = bound
            best_case = (shifts_u, shifts_v, shifts_uv, shifts_umv)
            print(f"plus-minus m={m} new best {best_bound:.6f} with {best_case}", flush=True)

    return best_bound, best_case


def main():
    for m in (3, 4, 5):
        print(f"\nSearching standard family for m={m}")
        bound, witness = search_standard_family(m)
        print(f"Best standard bound for m={m}: {bound:.6f} with {witness}")

        print(f"\nSearching plus-minus family for m={m}")
        bound, witness = search_plus_minus_family(m)
        print(f"Best plus-minus bound for m={m}: {bound:.6f} with {witness}")


if __name__ == "__main__":
    main()
