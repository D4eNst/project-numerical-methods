from decimal import Decimal, getcontext
from typing import Callable

import matplotlib.pyplot as plt

from utils import f_from_x, chebyshev_node, calculate_interpolation, interpolate_nodes, lagrange_interpolation, \
    newton_interpolation, erf_x_with_formula, inaccuracy

getcontext().prec = 10


class Solution:
    """
    Класс для решения задачи
    """

    def __init__(self,
                 a: Decimal,
                 b: Decimal,
                 n: int,
                 epsilon: Decimal = Decimal('1e-6'),
                 qn_func: Callable = None,
                 a0_func: Callable = None,
                 fi_from_t: Callable = None,
                 **kwargs
                 ):
        """
        Инициализация класса Solution.
        :param a: Начальное значение интервала.
        :param b: Конечное значение интервала.
        :param n: Количество узлов.
        :param epsilon: Точность вычислений.
        :param qn_func: Функция для вычисления qn. (1 задание)
        :param a0_func: Функция для вычисления a0. (1 задание)
        :param fi_from_t: Функция ф(t). (3 задание)
        :param kwargs: Дополнительные аргументы.
        :keyword num_points: количество точек между узлами для интерполяции (чтобы использовать по умолчанию для функций в классе)
        :keyword type_nodes: int: Тип распределения узлов (1 - равномерно, 2 - полином Чебышева).
        :keyword type_interpolation: int: Параметр для выбора функции интерполяции (1 - Логранжа, 2 - Ньютона)

        """
        valid_kwargs = {'type_nodes', 'type_interpolation', 'num_points'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")

        self.a = a
        self.b = b
        self.n = n
        self.h = Decimal((b - a) / (n - 1))
        self.epsilon = epsilon
        self.qn_func = qn_func
        self.a0_func = a0_func
        self.fi_from_t = fi_from_t

        ti = {
            1: lagrange_interpolation,
            2: newton_interpolation,
        }

        self.num_points = kwargs.get('num_points', 1)
        self.type_nodes = kwargs.get('type_nodes', 1)
        self.type_interpolation = ti[kwargs.get('type_interpolation', 1)]

        self._x_values_cache = None
        self._f_from_x_values_cache = None

    def calculate_uniform_nodes(self) -> list[Decimal]:
        """
        Возвращает список равномерно распределенных узлов для заданных параметров у объекта класса.
        Не имеет отношение к свойству x_values и всегда считается отдельно
        :return: Список значений x.
        """
        return [self.a + self.h * i for i in range(self.n)]

    def calculate_chebyshev_nodes(self):
        """
        Возвращает списо узлов полинома Чебышева для заданных параметров у объекта класса.
        Не имеет отношение к свойству `x_values` и всегда считается отдельно
        :return: Список значений x.
        """
        return [chebyshev_node(self.a, self.b, k, self.n) for k in range(self.n - 1, -1, -1)]

    @property
    def x_values(self) -> list[Decimal]:
        """
        Возвращает список значений x в интервале [a, b] с шагом h. Если

        :return: Список значений x.
        :raises Exception: Вызывает исключение при некорректом указании типа распределения узлов
        """
        node_calculations = {
            1: self.calculate_uniform_nodes,
            2: self.calculate_chebyshev_nodes
        }
        if self._x_values_cache is None:
            if self.type_nodes in node_calculations:
                self._x_values_cache = node_calculations[self.type_nodes]()
            else:
                raise Exception(
                    'Укажите один из доступных способов распределения узлов: \n'
                    ' 1) Равномерно\n 2) По полиному Чебышева')
        return self._x_values_cache

    @property
    def f_from_x_values(self) -> list[Decimal]:
        """
        Вычисляет значения функции f(x) для всех x_values.

        :return: Список результатов f(x).
        """
        if self._f_from_x_values_cache is None:
            self._f_from_x_values_cache = [f_from_x(x, self.epsilon, self.a0_func, self.qn_func) for x in self.x_values]
        return self._f_from_x_values_cache

    def interpolate_nodes(self, **kwargs) -> list[Decimal]:
        """
        Функция равномерно распределяет заданное количество точек `num_points` между всеми узлами

        :keyword num_points: Количество точек, которое необходимо вставить между узлами
        :return: Новый список узлов с дополнительными точкаами между узлами.
        """
        valid_kwargs = {'num_points'}
        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        num_points = kwargs.get('num_points', self.num_points)

        return interpolate_nodes(self.x_values, num_points)

    def calculate_interpolation(self, **kwargs) -> list[Decimal]:
        """
        Вычисляет приближенное значение функции для всех узлов проинтерполированного списка nodes
        на основе узлов полинома и значений в этих узлах.

        Примечение: для интерполирования узлов используется функция `interpolate_nodes()`,
        параметр `num_points` будет использоваться для этой функции

        :keyword num_points: Количество точек, которое необходимо вставить между узлами
        :keyword type_interpolation: Функция для интерполяции
        :return: Список, представляющий собой интерполяционный полином.
        :raise Exaption: Возникает если количество узлов не совпадает с количеством значчений
        """
        valid_kwargs = {'num_points', 'type_interpolation'}
        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        num_points = kwargs.get('num_points', self.num_points)
        type_interpolation = kwargs.get('type_interpolation', self.type_interpolation)

        return calculate_interpolation(self.x_values, self.f_from_x_values, num_points, type_interpolation)

    def inaccuracy(self, **kwargs) -> list[Decimal]:
        """
        Вычисляет погрешности для всех узлов проинтерполированного полинома
        :keyword num_points: Количество точек, которое необходимо вставить между узлами
        :return: список погрешностей
        """
        valid_kwargs = {'num_points'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        num_points = kwargs.get('num_points', self.num_points)

        return inaccuracy(
            self.x_values,
            self.f_from_x_values,
            self.a0_func,
            self.qn_func,
            self.epsilon,
            num_points
        )

    def error_vs_nodes_dependency(self,
                                  stop_on_large_er: bool = True,
                                  max_n: int = 100,
                                  **kwargs) -> list[Decimal]:
        """
        Вычисляет список значений максимальной погрешости в зависимости от количества узлов
        :param stop_on_large_er: Автоматически завершить расчет, Если погрешность превышает единицу
        :param max_n: Максимальное количество узлов (если stop_on_large_er задан, может игнорироваться)
        :keyword num_points: Количество точек, которое необходимо вставить между узлами
        :keyword start_n: По умолчанию текущее значение количества узлов, но может быть задано другое число > 2
        :keyword stop_n: Указывает при каком количестве узлов завершить вычисления
        :keyword step: Шаг для увеличения узлов
        :return:
        """
        valid_kwargs = {'num_points', 'start_n', 'stop_n', 'step'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")

        self_n = self.n
        num_points = kwargs.get('num_points', self.num_points)
        start_n = kwargs.get('start_n', self_n + 1)
        stop_n = kwargs.get('stop_n', max_n)
        step = kwargs.get('step', 2)

        max_inaccuracy_list = []
        while True:
            self.n = start_n
            inaccuracy_list = self.inaccuracy(num_points=num_points)
            max_inaccuracy = max(inaccuracy_list)
            max_inaccuracy_list.append(max_inaccuracy)
            start_n = start_n + step
            if abs(max_inaccuracy) > Decimal('1') and stop_on_large_er:
                print('Погрешность больше единицы, далее считать бессмысленно')
                print(f"n = {int(start_n) - 1}, Погрешность: {str(max_inaccuracy)}")
                break
            if start_n > max_n or start_n > stop_n:
                print(f'Погрешность мала, но превышено максиамльное число узлов: {start_n - 1}')
                break

        self.n = self_n
        return max_inaccuracy_list

    def get_normal_er(self, f_values: list[Decimal] = None, epsilon=Decimal('1e-3'), **kwargs) -> tuple[int, int]:
        """
        Вычисляет диапазон значений погрешности, где она оптимальна (меньше заданной)
        :param f_values: Значения погрешности (если не передан, вычисляется)
        :param epsilon: Допустимая погрешность

        :return: Первый и последний (включительно) индекс диапазона.
        В случае, если такого не существует, возвращает (-1, -1)
        """
        valid_kwargs = {'num_points'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        num_points = kwargs.get('num_points', self.num_points)

        f_val = self.inaccuracy(num_points=num_points) if f_values is None else f_values

        first = 0
        last = len(f_val) - 1

        while first <= last:
            if f_val[first] != Decimal('0') and abs(f_val[first]) < epsilon:
                break
            first += 1

        while last >= first:
            if f_val[last] != Decimal('0') and abs(f_val[last]) < epsilon:
                break
            last -= 1

        if last < first:
            return -1, -1
        return first, last

    def tabulate_erf_x(self, formula: str, epsilon: Decimal = Decimal('1e-4'), **kwargs) -> list[Decimal]:
        """
        Вычисляет значения интеграллов для всех диапазонов

        [ self.x_values[0], self.x_values[i] ], i = 1, 2, ... , self.n
        :param formula: Формула, по который вычисляется значние интегралла
        :param epsilon: Точность вычисления интеграла
        :keyword max_n: Максимальное число разбиений отрезка при вычислении интегралов
        :return: Список значений интеграллов
        """
        valid_kwargs = {'max_n'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        max_n = kwargs.get('max_n')

        return [erf_x_with_formula(self.fi_from_t, self.x_values[0], self.x_values[i], formula, epsilon, max_n)
                for i in range(1, self.n)]

    def show_console_tabulation(self):
        """
        Выводит протабуллированную функцию на атрезке [a, b] с шагом в консоль
        """
        print(f"{str('xi'):12} {str('fi'):12}")
        for xi, fi in zip(self.x_values, self.f_from_x_values):
            print(f"{str(xi):12} {str(fi):12}")

    def show_graph_tabulation(self):
        """
        Строит график протабуллированной функции на атрезке [a, b] с шагом в консоль
        """
        plt.plot([float(k) for k in self.x_values], [float(k) for k in self.f_from_x_values], marker='o', linestyle='-')

        plt.title('График функции')
        plt.xlabel('x\u1d62')
        plt.ylabel('f\u1d62')

        plt.grid(True)
        plt.figtext(0.17, 0.82, f'a = {self.a}', ha='left')
        plt.figtext(0.17, 0.77, f'b = {self.b}', ha='left')
        plt.figtext(0.17, 0.72, f'h = {self.h}', ha='left')
        plt.show()

    def show_graph_inaccuracy(self, only_normal_er: bool = False, **kwargs) -> None:
        """
        Построение графика погрешностей
        :param only_normal_er: Получение только диапазона ошибок, меньше epsilon

        """
        valid_kwargs = {'num_points', 'epsilon'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        num_points = kwargs.get('num_points', 1)
        epsilon = kwargs.get('epsilon', Decimal('1e-3'))

        f_val = self.inaccuracy(num_points=num_points)
        x_val = interpolate_nodes(self.x_values, num_points=num_points)

        if only_normal_er:
            first, last = self.get_normal_er(f_values=f_val, epsilon=epsilon, num_points=num_points)
            if last >= 0:
                x_val = x_val[first:last + 1]
                f_val = f_val[first:last + 1]
            else:
                raise Exception(f'Не удается получить диапазон погрешностей, меньше чем epsilon: {str(epsilon)}\n'
                                f'Наименьшая погрешность: {str(min(f_val))}')

        plt.plot([float(x) for x in x_val], [float(f) for f in f_val])

        plt.title('График погрешностей')
        plt.xlabel('x\u1d62')
        plt.ylabel('Δ\u1d62')
        plt.grid(True)
        plt.figtext(0.88, 0.92, f'n = {self.n}', ha='left')
        plt.figtext(0.88, 0.82, f'a = {round(x_val[0], 4)}', ha='left')
        plt.figtext(0.88, 0.72, f'b = {round(x_val[-1], 4)}', ha='left')
        plt.show()

    def show_graph_error_vs_nodes_dependency(self, **kwargs):
        """
        Строит график зависимости максимальной погрешнотси он количества узлов
        :param kwargs: Параметры для функции error_vs_nodes_dependency
        :return:
        """

        valid_kwargs = {'step', 'stop_on_large_er', 'max_n', 'num_points', 'start_n', 'stop_n'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        step = kwargs.get('step', 2)

        max_inaccuracy_list = self.error_vs_nodes_dependency(**kwargs)
        start_n = kwargs.get('start_n', Decimal(f'{self.n + 1}'))
        n_list = [float(start_n) + i * step for i in range(len(max_inaccuracy_list))]

        plt.plot(n_list, [float(er) for er in max_inaccuracy_list])

        plt.title('График погрешностей от кол-ва узлов')
        plt.xlabel('n\u1d62')
        plt.ylabel('Δ\u1d62')
        plt.grid(True)
        plt.show()

    def show_graph_erf_x_error(self, *args: str, methods: str | list[str] = 'all', **kwargs) -> None:
        """
        :param args: Формулы для вычилений интегралов. Для каждой формулы строится свой график
        :param methods: Формулы можно передать в виде списка или один метод в виде строки. Для указания
        всех методов можно использовать строку "all". В случае, если формулы были переданы в args, будут
        использоваться формулы из обоих параметров (если methods == "all", и args заданы,
        будут использованы только args)
        :keyword max_n: Максимальное число разбиений отрезка при вычислении интегралов
        :keyword epsilon: Точность вычисления интеграла
        :return:
        """
        valid_kwargs = {'max_n', 'epsilon'}

        for key in kwargs:
            if key not in valid_kwargs:
                raise ValueError(f"Недопустимый аргумент в kwargs: {key}")
        max_n = kwargs.get('max_n')
        epsilon = kwargs.get('epsilon', Decimal('1e-4'))

        if methods == 'all':
            if not args:
                methods = ["left_rect", "right_rect", "center_rect", "trapezoid", "simpson", "gauss"]
            else:
                methods = args
        else:
            if isinstance(methods, str):
                methods = [methods]
            methods = set(methods + list(args))

        print(methods)
        arrays = [[Decimal(0.0)] for _ in range(len(methods))]
        for i, method in enumerate(methods):
            arrays[i] += self.tabulate_erf_x(method, epsilon=epsilon, max_n=max_n)

        arrays = [list(map(float, [abs(f2 - f1) for f2, f1 in zip(self.f_from_x_values, array)])) for array in arrays]
        x_val = [float(x) for x in self.x_values]

        for array in arrays:
            plt.plot(x_val, array)

        plt.title('График погрешностей при вычислении интеграла')
        plt.xlabel('x')
        plt.ylabel('Δ\u1d62')
        plt.grid(True)
        plt.show()

    def __setattr__(self, name, value):
        if name in vars(self).keys() and name[0] != "_":
            if getattr(self, f"{name}", None) != value:
                super().__setattr__(name, value)
                self._reset_cache()
                self.h = Decimal((self.b - self.a) / (self.n - 1))
        else:
            super().__setattr__(name, value)

    def _reset_cache(self):
        self._f_from_x_values_cache = None
        self._x_values_cache = None
