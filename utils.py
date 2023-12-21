from typing import Callable, Union
from decimal import Decimal
from math import cos, pi, sqrt


def chebyshev_node(a: Decimal, b: Decimal, k: int, n: int) -> Decimal:
    """
    Вычисляет значение узла Чебышева x_k для заданных параметров.

    :param a: Начальное значение отрезка
    :param b: Конечное значение отрезка
    :param k: Порядковый номер узла (целое число)
    :param n: Количество узлов полинома

    :return: Значение узла Чебышева x в узле k
    """
    return Decimal('0.5') * (a + b) + \
        Decimal('0.5') * (b - a) * Decimal(cos(Decimal((2 * k + 1) * pi) / Decimal((2 * n))))


def f_from_x(x: Decimal,
             epsilon: Decimal,
             a0_func: Callable[[Decimal], Decimal],
             qn_func: Callable[[Decimal, int], Decimal]) -> Decimal:
    """
    Вычисляет значение функции f(x) с использованием метода численного анализа.

    :param x: Значение x
    :param epsilon: Погрешность для вычисления значения
    :param a0_func: Функция для расчета значения a0 при i = 0. Принимает в качестве параметра значение x
    :param qn_func: Функция вычисления множителя qn для шага i.
    Первым параметром принимает значение x, вторым - номер шага i

    :return: Результат вычисления f(x).

    :raises Exception: Бросает исключение при недостаточной сходимости.

    """
    res = 0
    an = None
    an_1 = a0_func(x)
    i = 0
    while an is None or abs(an_1 - an) > epsilon:
        res = res + an_1
        qn = qn_func(x, i)
        an = an_1
        an_1 *= qn

        i += 1
        if i > 1000:
            raise Exception("Недостаточная сходимость")
    return res


def interpolate_nodes(nodes: list[Decimal], num_points: int = 1) -> list[Decimal]:
    """
    Функция равномерно распределяет заданное количество точек `num_points` между всеми узлами в списке `nodes`

    Если в функии кол-во узлов n = 5, [a, b] = [0, 12] и они распределены равномерно,
    узлы, будут иметь ззначение: 0, 3, 6, 9, 12.
    При вызове функции interpolate_nodes(nodes, num_points=1) получим полином с n = 9,
    значение узлов: 0, 1.5, 3, 4.5, 7, 8.5, 9, 10.5, 12

    :param nodes: Список узлов
    :param num_points: Количество точек, которое необходимо вставить между узлами, по умолчанию 1
    :return: Новый список узлов с дополнительными точкаами между узлами.
    """
    new_x_values = [nodes[0]]
    for i in range(len(nodes) - 1):
        x1 = nodes[i]
        x2 = nodes[i + 1]
        step = (x2 - x1) / (num_points + 1)
        next_x = x1 + step
        for _ in range(num_points):
            new_x_values.append(next_x)
            next_x = next_x + step
        new_x_values.append(nodes[i + 1])
    return new_x_values


def newton_interpolation(x: Decimal, nodes: list[Decimal], f_values: list[Decimal]) -> Decimal:
    """
    Вычисляет приближенное значение функции для заданного значения x
    на основе узлов полинома и значений в этих узлах
    :param x: Значение x.
    :param nodes: Список узлов
    :param f_values: Список значений в заданных узлах
    :return: Приближенное значение f(x) - Запрограммировано по формуле Логранжа.
    :raise Exaption: Возникает если количество узлов не совпадает с количеством значчений
    """
    n = len(nodes)
    if len(f_values) != n:
        raise Exception("Количество узлов должно совпадать с количеством значений")

    divided_differences = [f_values[:]]
    for i in range(1, n):
        divided_differences.append([])
        for j in range(n - i):
            divided_difference = (divided_differences[i - 1][j + 1] - divided_differences[i - 1][j]) / (
                    nodes[j + i] - nodes[j])
            divided_differences[i].append(divided_difference)

    result = f_values[0]
    term = 1
    for i in range(1, n):
        term *= (x - nodes[i - 1])
        result += term * divided_differences[i][0]

    return result


def lagrange_interpolation(x: Decimal, nodes: list[Decimal], f_values: list[Decimal]) -> Decimal:
    """
    Вычисляет приближенное значение функции для заданного значения x
    на основе узлов полинома и значений в этих узлах
    :param x: Значение x.
    :param nodes: Список узлов
    :param f_values: Список значений в заданных узлах
    :return: Приближенное значение f(x) - Запрограммировано по формуле Логранжа.
    :raise Exaption: Возникает если количество узлов не совпадает с количеством значчений
    """
    if len(nodes) != len(f_values):
        raise Exception('Количество узлов в списке x_values должно совпадать с количеством значений f_from_x_values')

    result = Decimal('0.0')
    for i in range(len(nodes)):
        # term = f(xi)*Li(x)
        li_x = 1
        for j in range(len(nodes)):
            if i != j:
                li_x *= (x - nodes[j]) / (nodes[i] - nodes[j])

        term = f_values[i] * li_x
        result += term
    return result


def calculate_interpolation(
        nodes: list[Decimal],
        f_values: list[Decimal],
        num_points: int = 1,
        method_interpolation: Callable = lagrange_interpolation) -> list[Decimal]:
    """
    Вычисляет приближенное значение функции для всех узлов проинтерполированного списка nodes
    на основе узлов полинома и значений в этих узлах.

    Примечение: для интерполирования узлов используется функция `interpolate_nodes()`,
    параметр `num_points` будет использоваться для этой функции

    :param nodes: Список узлов
    :param f_values: Список значений в заданных узлах
    :param num_points: Количество точек, которое необходимо вставить между узлами, по умолчанию 1
    :param method_interpolation: Функция для интерполяции, по умолчанию использует полином Лагранжа
    :return: Приближенное значение f(x) - Запрограммировано по формуле Логранжа.
    :raise Exaption: Возникает если количество узлов не совпадает с количеством значчений

    :return: Список, представляющий собой интерполяционный полином.
    """
    if len(nodes) != len(f_values):
        raise Exception('Количество узлов в списке x_values должно совпадать с количеством значений f_from_x_values')

    interpolated_nodes = interpolate_nodes(nodes, num_points)
    return [method_interpolation(x, nodes, f_values) if i % (num_points + 1) != 0 else f_values[i // (num_points + 1)]
            for i, x in enumerate(interpolated_nodes)]


def inaccuracy(nodes: list[Decimal],
               f_values: list[Decimal],
               a0_func: Callable[[Decimal], Decimal],
               qn_func: Callable[[Decimal, int], Decimal],
               epsilon: Decimal,
               num_points: int = 1) -> list[Decimal]:
    """
    Вычисляет погрешности для всех узлов проинтерполированного полинома
    :param nodes: Список узлов
    :param f_values: Список значений в заданных узлах
    :param epsilon: Погрешность для вычисления значения
    :param a0_func: Функция для расчета значения a0 при i = 0. Принимает в качестве параметра значение x
    :param qn_func: Функция вычисления множителя qn для шага i.
    :param num_points: Количество точек, которое необходимо вставить между узлами
    :return:
    """

    result = []

    x_values_er = interpolate_nodes(nodes, num_points)
    f_from_x_values_er = calculate_interpolation(nodes, f_values, num_points)

    i = 0
    j = 0
    while i < len(nodes):
        er_x = x_values_er[j]
        if nodes[i] == er_x:
            result.append(Decimal('0'))
            i += 1
            j += 1
        else:
            f_from_x_er = f_from_x_values_er[j]
            f_from_x_normal = f_from_x(er_x, epsilon, a0_func, qn_func)
            result.append(abs(f_from_x_normal - f_from_x_er))
            j += 1
    return result


def erf_x_with_formula(fi_from_t: Callable[[Decimal], Decimal],
                       x1: Decimal,
                       x2: Decimal,
                       formula: str,
                       epsilon: Decimal,
                       max_n: int = None) -> Decimal:
    """
    Вычисление приближенного значения значения erf(x), используя составную квадратурную формулу Гаусса с двумя узлами
    :param max_n:
    :param formula: Формула
    :param fi_from_t: Функция ф(t)
    :param x1: Первое значение (a)
    :param x2: Второе значение (b)
    :param epsilon: Точность вычисления интегралла
    :return: Приближенное значение интегралла
    """

    formula_dict = {
        "gauss": gauss_si,
        "left_rect": left_rect_si,
        "right_rect": right_rect_si,
        "center_rect": center_rect_si,
        "trapezoid": trapezoid_si,
        "simpson": simpson_si
    }

    formula_si = formula_dict.get(formula)
    if formula_si is None:
        raise Exception(
            f'Недопустимая формула: {formula}.\nДопустимые формулы:\n{", ".join([f for f in formula_dict])}')

    n = 2
    sn: Union[None, Decimal] = None
    s2n: Union[None, Decimal] = None
    while s2n is None or abs(s2n - sn) > epsilon:
        hn = (x2 - x1) / n
        s = Decimal('0')
        for i in range(1, n + 1):
            si = formula_si(fi_from_t, x1, i, hn)
            s += si

        if sn is None:
            sn = s
        else:
            sn = s2n if s2n is not None else sn
            s2n = s

        n *= 2
        if n > 500000:
            e_string = f'Не получается получить интеграл на отрезке [{x1}, {x2}] с заданой точностью ({epsilon}). \n' \
                       f'N = {n}.\nТекущая точность ( |Sn(ф) - S2n(ф)| ) = {s2n - sn}.\nS2n(ф) = {s2n}.'
            raise Exception(e_string)
        if max_n is not None and n > max_n:
            print(f"{n}")
            return s2n
    print(f"{n}")
    return s2n


def gauss_si(fi_from_t, x1, i, hn):
    zi_1 = x1 + (i - 1) * hn
    t1 = zi_1 + (hn / 2) * (1 - (1 / Decimal(sqrt(3))))
    t2 = zi_1 + (hn / 2) * (1 + (1 / Decimal(sqrt(3))))

    si = (hn / 2) * (fi_from_t(t1) + fi_from_t(t2))
    return si


def left_rect_si(fi_from_t, x1, i, hn):
    b = x1 + i * hn
    a = b - hn
    si = (b - a) * fi_from_t(a)
    return si


def right_rect_si(fi_from_t, x1, i, hn):
    b = x1 + i * hn
    a = b - hn
    si = (b - a) * fi_from_t(b)
    return si


def center_rect_si(fi_from_t, x1, i, hn):
    b = x1 + i * hn
    a = b - hn
    si = (b - a) * fi_from_t((b + a) / 2)
    return si


def trapezoid_si(fi_from_t, x1, i, hn):
    b = x1 + i * hn
    a = b - hn
    si = (b - a) / 2 * (fi_from_t(a) + fi_from_t(b))
    return si


def simpson_si(fi_from_t, x1, i, hn):
    b = x1 + i * hn
    a = b - hn
    c = (a + b) / 2
    si = (b - a) / 6 * (fi_from_t(a) + 4 * fi_from_t(c) + fi_from_t(b))
    return si


