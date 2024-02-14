Описание проекта

В рамках учебного проекта, я занялся исследованием численных методов, а именно интерполированиием трансцендентных функций и вычислением интеграла с использованием квадратурной формы. 
Этот проект написан на языке программирования Python, при этом акцент сделан на использование базовых возможностей языка и объекта Decimal для точных вычислений.

По самому проекту есть подробный отчет, который находится в этом репозитории.

Основныые технологии:
- Python / Git
- Decimal
- Matplotlib

# Пример использования

```
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
   qn_func=lambda x, i: -1 * (4 * i + 1) * (Decimal(pi) / 2) ** 2 * (x ** 4) /
                        ((2 * i + 1) * (2 * i + 2) * (4 * i + 5)),
   a0_func=lambda x: x,
   fi_from_t=lambda t: Decimal(cos((Decimal.from_float(pi) * t * t) / 2)),
   num_points=2,
   type_nodes=1,
   type_interpolation=1)
if __name__ == "__main__":
   main()
```

Вывести значения узлов и посчитанные по первому заданию значения в этих узлах:
```
solution1.show_console_tabulation()
```

Вот так можно вывести в столбик все значения интерполяционного полинома:
```
print("Полином")
print(*map(str, solution1.calculate_interpolation()), sep="\n")
```

При этом также можно вывести иксы этого полинома:
```
print("x")
print(*solution1.interpolate_nodes(), sep="\n")
```

Список погрешностей можно вывести так:
```
print(*solution1.inaccuracy())
А для того, чтобы вывести график, как в отчете, можно использовать метод:
solution1.show_graph_inaccuracy()
```
Или
```
solution1.show_graph_inaccuracy(only_normal_er=True, epsilon=Decimal('1e-4'))
```

Таблицы из части 4 я составлял с помощью такого кода:
```
erf_xes = [Decimal('0.0')] + solution1.tabulate_erf_x(formula='gauss')
print(" C(x) erf(x) err")
for erf_x, f in zip(erf_xes, solution1.f_from_x_values):
 print(f"{f:2.12f} {erf_x:2.12f} {abs(f - erf_x):2.10f}")
```

А графики можно вывести так:
```
solution1.show_graph_erf_x_error(methods=["center_rect", "trapezoid", "simpson", 
"gauss"])
```
Или
```
solution1.show_graph_erf_x_error(methods= "all")
```

Вывод списка максимальных погрешностей до тех пор, пока она не превысит единицу, с шагом 2:
```
print(f"\n-----------------------------------\nЗависимость погрешности от количества 
узлов")
dict_for_max_inaccuracy = {
 'num_points': 2,
 'stop_on_large_er': True,
 'max_n': 107,
 'start_n': 6,
 'stop_n': 100,
 'step': 2
}
max_inaccuracy_list = solution1.error_vs_nodes_dependency(**dict_for_max_inaccuracy)
print("n max_er")
for i in range(len(max_inaccuracy_list)):
 print(f'{int(dict_for_max_inaccuracy.get("start_n", n + 1)) + i*10:2}, '
 f'{str(max_inaccuracy_list[i])}')
```

Теперь график:
```
solution1.show_graph_error_vs_nodes_dependency(**dict_for_max_inaccuracy)
```

