from math import pi, cos
from decimal import Decimal
from solution import Solution


def main():
    a = Decimal('0')
    b = Decimal('1.5')
    n = 16
    epsilon = Decimal('1e-6')

    solution1 = Solution(
        a=a, b=b, n=n, epsilon=epsilon,
        qn_func=lambda x, i: -1 * (4 * i + 1) * (Decimal(pi) / 2) ** 2 * (x ** 4) / (
                (2 * i + 1) * (2 * i + 2) * (4 * i + 5)),
        a0_func=lambda x: x,
        fi_from_t=lambda t: Decimal(cos((Decimal.from_float(pi) * t * t) / 2)),
        num_points=2,
        type_nodes=1,
        type_interpolation=1)
    # solution1 = Solution(
    #     a=a, b=b, n=n, epsilon=epsilon,
    #     qn_func=lambda x, i: (-1) * ((x*x*(2*i+1))/((2*i+2)*((2*i+3)*(2*i+3)))),
    #     a0_func=lambda x: x,
    #     fi_from_t=lambda t: Decimal(cos((Decimal.from_float(pi) * t * t) / 2)),
    #     num_points=1,
    #     type_nodes=1,
    #     type_interpolation=2)

    # solution1.show_console_tabulation()
    # print("Полином")
    # print(*map(str, solution1.calculate_interpolation()), sep="\n")

    # print("x")
    # print(*solution1.interpolate_nodes(), sep="\n")
    # print("fi")
    # print(*[f_from_x(x, solution1.epsilon, solution1.a0_func, solution1.qn_func) for x in solution1.interpolate_nodes()], sep="\n")
    # print(*map(str, solution1.inaccuracy()), sep="\n")
    # print(*solution1.calculate_interpolation(), sep="\n")

    print(*solution1.inaccuracy())
    # solution1.show_graph_inaccuracy(only_normal_er=True, epsilon=Decimal('1e-4'))
    # solution1.show_graph_inaccuracy()

    # print(" n")
    # print(f"{0:>3}")
    # erf_xes = [Decimal('0.0')] + solution1.tabulate_erf_x(formula='gauss')
    # print("      C(x)          erf(x)         err")
    # for erf_x, f in zip(erf_xes, solution1.f_from_x_values):
    #     print(f"{f:2.12f} {erf_x:2.12f} {abs(f - erf_x):2.10f}")

    # solution1.show_graph_erf_x_error(methods=["center_rect", "trapezoid", "simpson", "gauss"])
    # print(*solution1.f_from_x_values)
    # print(solution1.tabulate_erf_x('simpson', Decimal('1e-4')))

    # print(f"\n-----------------------------------\nЗависимость погрешности от количества узлов")
    # dict_for_max_inaccuracy = {
    #     'num_points': 2,
    #     'stop_on_large_er': True,
    #     'max_n': 107,
    #     'start_n': 6,
    #     'stop_n': 100,
    # }
    # max_inaccuracy_list = solution1.error_vs_nodes_dependency(**dict_for_max_inaccuracy)
    # print(max_inaccuracy_list)
    # print(len(max_inaccuracy_list))
    #
    # print("n    max_er")
    # for i in range(len(max_inaccuracy_list)):
    #     print(f'{int(dict_for_max_inaccuracy.get("start_n", n + 1)) + i*10:2}, '
    #           f'{str(max_inaccuracy_list[i])}')

    # solution1.show_graph_error_vs_nodes_dependency(**dict_for_max_inaccuracy)
    #
    # solution1.n = 1000

    # solution1.show_graph_inaccuracy(True, rate=Decimal('10'), epsilon=Decimal('1e-2'))
    # solution1.show_graph_tabulation()
    # solution1.show_graph_inaccuracy(only_normal_er=False, num_points=1, epsilon=Decimal('1e-3'))

    # print(erf_x_with_formula(
    #     fi_from_t=lambda t: Decimal(cos((Decimal.from_float(pi) * t * t) / 2)),
    #     x1=Decimal('0'),
    #     x2=Decimal('0.3'),
    #     formula="simpson",
    #     epsilon=Decimal('1e-4')
    # ))

    # def fi_from_t(t: Decimal):
    #     return Decimal(cos((Decimal.from_float(pi) * t * t) / 2))


if __name__ == "__main__":
    main()
