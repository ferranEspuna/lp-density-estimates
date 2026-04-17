import itertools

from two_parameter_lp import (
    evaluate_two_parameter_bound,
    make_direction_family,
    standard_family,
)


def plus_minus_family(m, shifts_u, shifts_v, shifts_uv, shifts_umv):
    return (
        make_direction_family((1, 0), shifts_u, prefix="u")
        + make_direction_family((0, 1), shifts_v, prefix="v")
        + make_direction_family((1, 1), shifts_uv, dilation=m, prefix="uv")
        + make_direction_family((1, -1), shifts_umv, dilation=m, prefix="umv")
    )


def multi_direction_family(m, direction_specs):
    forms = []
    for direction, shifts, dilation, prefix in direction_specs:
        forms.extend(make_direction_family(direction, shifts, dilation=dilation, prefix=prefix))
    return forms


def _shift_options(max_shift=2, max_size=2):
    options = []
    for size in range(1, max_size + 1):
        options.extend(itertools.combinations(range(max_shift + 1), size))
    return options


def search_standard_family(m, max_shift=2, max_size=2):
    best_bound = 1.0
    best_case = None

    options = _shift_options(max_shift=max_shift, max_size=max_size)

    for shifts_u, shifts_v, shifts_uv in itertools.product(options, repeat=3):
        forms = standard_family(m, shifts_u, shifts_v, shifts_uv)
        bound = evaluate_two_parameter_bound(m, forms)
        if bound is None:
            continue
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
        if bound is None:
            continue
        if bound < best_bound - 1e-9:
            best_bound = bound
            best_case = (shifts_u, shifts_v, shifts_uv, shifts_umv)
            print(f"plus-minus m={m} new best {best_bound:.6f} with {best_case}", flush=True)

    return best_bound, best_case


def search_extended_directions(m, max_shift=2, max_size=2):
    best_bound = 1.0
    best_case = None
    options = _shift_options(max_shift=max_shift, max_size=max_size)

    direction_pool = [
        ((1, 0), 1, "u"),
        ((0, 1), 1, "v"),
        ((1, 1), m, "uv"),
        ((1, -1), m, "umv"),
        ((2, 1), m, "2u_v"),
        ((1, 2), m, "u_2v"),
        ((2, -1), m, "2u_mv"),
        ((1, -2), m, "u_m2v"),
    ]

    for extra_directions in itertools.combinations(direction_pool[2:], 2):
        chosen = direction_pool[:2] + list(extra_directions)
        for shift_choice in itertools.product(options, repeat=len(chosen)):
            direction_specs = []
            for (direction, dilation, prefix), shifts in zip(chosen, shift_choice):
                direction_specs.append((direction, shifts, dilation, prefix))

            forms = multi_direction_family(m, direction_specs)
            bound = evaluate_two_parameter_bound(m, forms)
            if bound is None:
                continue
            if bound < best_bound - 1e-9:
                best_bound = bound
                best_case = direction_specs
                print(
                    f"extended m={m} new best {best_bound:.6f} with {best_case}",
                    flush=True,
                )

    return best_bound, best_case


def main():
    for m in (3,):
        print(f"\nSearching standard family for m={m}")
        bound, witness = search_standard_family(m)
        print(f"Best standard bound for m={m}: {bound:.6f} with {witness}")

        print(f"\nSearching plus-minus family for m={m}")
        bound, witness = search_plus_minus_family(m)
        print(f"Best plus-minus bound for m={m}: {bound:.6f} with {witness}")

        print(f"\nSearching extended-direction family for m={m}")
        bound, witness = search_extended_directions(m, max_shift=5, max_size=2)
        print(f"Best extended bound for m={m}: {bound:.6f} with {witness}")


if __name__ == "__main__":
    main()
