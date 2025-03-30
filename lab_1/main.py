import matplotlib.pyplot as plt
import math
from methods import *

def func_3(x, u):
    return x + u ** 3


def analit_1(x):
    return math.cos(math.sqrt(x))

def teilor_1(x):
    return 1 - 1 / 2 * x + 1 / 24 * x ** 2 - 1 / 720 * x ** 3


def analit_2(u):
    return math.exp(u * u) - u * u / 2 - 0.5

def pikar_first_2(u):
    return 0.5 + u ** 2 / 2 + u ** 4 / 4

def pikar_second_2(u):
    return 0.5 + u ** 2 / 2 + u ** 4 / 2 + u ** 6 / 12

def pikar_third_2(u):
    return 0.5 + u ** 2 / 2 + u ** 4 / 2 + u ** 6 / 6 + u ** 8 / 48

def pikar_forth_2(u):
    return 0.5 + u ** 2 / 2 + u ** 4 / 2 + u ** 6 / 6 + u ** 8 / 24 + u ** 10 / 240


def pikar_first_3(x):
    return x ** 2 / 2

def pikar_second_3(x):
    return x ** 7 / 56 + x ** 2 / 2

def pikar_third_3(x):
    return (x ** 22 / 3863552
            + 3 * x ** 17 / 106624
            + x ** 12 / 896
            + x ** 7 / 56
            + x ** 2 / 2)

def pikar_forth_3(x):
    return (x ** 67 / 3863981943017739124736
            + 9 * x ** 62 / 98677964914244452352
            + 1081 * x ** 57 / 73440052228804050944
            + 76623 * x ** 52 / 53388985337387155456
            + 310503 * x ** 47 / 3244062457475366912
            + (1069 * x ** 42) / 224099368435712
            + (790667 * x ** 37) / 4145838316060672
            + (11871 * x ** 32) / 1883187970048
            + (361 * x ** 27) / 2101772288
            + (47 * x ** 22) / 11941888
            + (33 * x ** 17) / 426496
            + (x ** 12) / 896
            + (x ** 7) / 56
            + (x ** 2) / 2)



def table_print(x, y1, y2, y3, y4, y5):
    print("------------------------------------------------------------------------------------------")
    print("|      x      |    Пикар 1   |    Пикар 2   |    Пикар 3   |    Пикар 4   |     Эйлер    |")
    print("------------------------------------------------------------------------------------------")
    for i in range(len(x)):
        print("|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|{:^14.5f}|".format(
            x[i], y1[i], y2[i], y3[i], y4[i], y5[i]))
    print("------------------------------------------------------------------------------------------")

def task1():
    x0 = 0
    x_max = 20
    N = 1000
    h = (x_max - x0) / N
    Xvals = []
    Yvals_analit = []
    Yvals_teilor = []
    Yvals_euler = euler_task1(0, 1, -0.5, h, N)
    for i in range(N):
        x = x0 + i * h
        Xvals.append(x)
        Yvals_analit.append(analit_1(x))
        Yvals_teilor.append(teilor_1(x))
    plt.plot(Xvals, Yvals_analit, label='Аналитическое')
    plt.plot(Xvals, Yvals_teilor, label='Тейлор')
    plt.plot(Xvals, Yvals_euler, label='Эйлер')
    plt.legend()
    plt.show()

def task2():
    x0 = 0
    x_max = 1
    N = 1000
    h = (x_max - x0) / N
    Xvals = []
    Yvals_analit = []
    Yvals_first = []
    Yvals_second = []
    Yvals_third = []
    Yvals_forth = []
    for i in range(N):
        x = x0 + i * h
        Xvals.append(x)
        Yvals_analit.append(analit_2(x))
        Yvals_first.append(pikar_first_2(x))
        Yvals_second.append(pikar_second_2(x))
        Yvals_third.append(pikar_third_2(x))
        Yvals_forth.append(pikar_forth_2(x))
    plt.plot(Yvals_analit, Xvals, label='Аналитическое')
    plt.plot(Yvals_first, Xvals, label='1 приближение')
    plt.plot(Yvals_second, Xvals, label='2 приближение')
    plt.plot(Yvals_third, Xvals, label='3 приближение')
    plt.plot(Yvals_forth, Xvals, label='4 приближение')
    plt.legend()
    plt.show()

def task3(x_max):
    x0 = 0
    print("Максимальное значение x_max: ", x_max)
    N = 100
    h = (x_max - x0) / N
    Xvals = []
    Yvals_first = []
    Yvals_second = []
    Yvals_third = []
    Yvals_forth = []
    Yvals_euler = euler_base(0, 0, h, N, func_3)
    for i in range(N):
        x = x0 + i * h
        Xvals.append(x)
        Yvals_first.append(pikar_first_3(x))
        Yvals_second.append(pikar_second_3(x))
        Yvals_third.append(pikar_third_3(x))
        Yvals_forth.append(pikar_forth_3(x))
    table_print(Xvals, Yvals_first, Yvals_second, Yvals_third, Yvals_forth, Yvals_euler)
    plt.plot(Xvals, Yvals_euler, label='Эйлер')
    plt.plot(Xvals, Yvals_first, label='1 приближение')
    plt.plot(Xvals, Yvals_second, label='2 приближение')
    plt.plot(Xvals, Yvals_third, label='3 приближение')
    plt.plot(Xvals, Yvals_forth, label='4 приближение')
    plt.legend()
    plt.show()

def calc_xmax():
    lhs, rhs = 0, 2
    accur = 1e-15
    prev = lhs
    cur = (rhs - lhs) / 2
    while abs(prev - cur) > accur:
        try:
            N = 100
            h = cur / N
            euler_base(0, 0, h, N, func_3)
            lhs = cur
        except OverflowError:
            rhs = cur
        prev = cur
        cur = (rhs + lhs) / 2
    return lhs


while True:
    action = int(input("Введите номер задания или 0, если хотите завершить работу: "))
    while action < 0 or action > 3:
        action = int(input("Введите корректный номер задания (1, 2 или 3) или 0, если хотите завершить работу: "))

    if action == 1:
        task1()
    elif action == 2:
        task2()
    elif action == 3:
        task3(calc_xmax())
    elif action == 0:
        break
