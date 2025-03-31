import math as m
import matplotlib.pyplot as plt
from runge import rungeKuttaSystem
from decimal import Decimal, getcontext

# Установка точности для Decimal
getcontext().prec = 37

variant = int(input("Введите вариант задания (1-2): "))
while variant != 1 and variant != 2:
    variant = int(input("Введите вариант задания (1-2): "))

if variant == 1:
    table_in = [ [ Decimal(2000),  Decimal(8.2e-3)    ],
                  [ Decimal(3000),  Decimal(2.768e-2) ],
                  [ Decimal(4000),  Decimal(6.56e-2)  ],
                  [ Decimal(5000),  Decimal(1.281e-1) ],
                  [ Decimal(6000),  Decimal(2.214e-1) ],
                  [ Decimal(7000),  Decimal(3.516e-1) ],
                  [ Decimal(8000),  Decimal(5.248e-1) ],
                  [ Decimal(9000),  Decimal(7.472e-1) ],
                  [ Decimal(10000), Decimal(1.025)    ] ]
else:
    table_in = [ [ Decimal(2000),  Decimal(1.6)    ],
                  [ Decimal(3000),  Decimal(5.4)   ],
                  [ Decimal(4000),  Decimal(12.8)  ],
                  [ Decimal(5000),  Decimal(25)    ],
                  [ Decimal(6000),  Decimal(43.2)  ],
                  [ Decimal(7000),  Decimal(68.6)  ],
                  [ Decimal(8000),  Decimal(102.4) ],
                  [ Decimal(9000),  Decimal(145.8) ],
                  [ Decimal(10000), Decimal(200)   ] ]

table_with_replaces = [[Decimal(m.log(float(measure))) for measure in row] for row in table_in]

R0 = Decimal(0)
R = Decimal(0.35)
Tw = Decimal(2000)
T0 = Decimal(10000)
p = Decimal(4)
EPS = Decimal(1e-36)
c = Decimal(3 * 10 ** 10)

def calcT(r):
    return (Tw - T0) * (r / R) ** p + T0

def calcPlanca(r):
    return Decimal(3.084e-4) / (Decimal(m.exp(Decimal(4.799 * 10 ** 4) / calcT(r))) - Decimal(1))

def calcK(table, r):
    Tlog = Decimal(m.log(calcT(r)))
    klog = table[0][1] + (Tlog - table[0][0]) / (table[1][0] - table[0][0]) * (table[1][1] - table[0][1])
    return Decimal(m.exp(klog))

def fFunc(r, u, f):
    return (-Decimal(3)) * calcK(table_with_replaces, r) * f / c

def phiFunc(r, u, f):
    if r == 0:
        return c * calcK(table_with_replaces, r) * (calcPlanca(r) - u) / Decimal(2)
    return -(f / r) + c * calcK(table_with_replaces, r) * (calcPlanca(r) - u)

def calcSecondCond(u, v):
    return v[-1] - Decimal(0.39) * c * u[-1]

def sign(val):
    if abs(val) < EPS:
        return 0
    elif val < 0:
        return -1
    return 1

def shotMethodSolve():
    n = 100
    ksi_l, ksi_r = Decimal(0), Decimal(1)
    u0_l, u0_r = ksi_l * calcPlanca(R0), ksi_r * calcPlanca(R0)
    F0 = Decimal(0)

    x_solve, u_l, f_l = rungeKuttaSystem(Decimal(0), u0_l, F0, (R - R0) / n, n, fFunc, phiFunc)
    u_r, f_r = rungeKuttaSystem(Decimal(0), u0_r, F0, (R - R0) / n, n, fFunc, phiFunc)[1:]
    l_sign = sign(calcSecondCond(u_l, f_l))
    r_sign = sign(calcSecondCond(u_r, f_r))

    u_solve, f_solve = [], []
    iter = 1
    while abs((ksi_l - ksi_r) / ((ksi_l + ksi_r) / Decimal(2))) > EPS:
        ksi_mid = (ksi_l + ksi_r) / Decimal(2)
        u0_mid = ksi_mid * calcPlanca(Decimal(0))
        u_solve, f_solve = rungeKuttaSystem(Decimal(0), u0_mid, F0, (R - R0) / n, n, fFunc, phiFunc)[1:]
        mid_sign = sign(calcSecondCond(u_solve, f_solve))
        print(f'Итерация {iter}; Знак {mid_sign}; значение ksi = {ksi_mid}')
        iter += 1
        if mid_sign == l_sign:
            ksi_l = ksi_mid
        else:
            ksi_r = ksi_mid

    u_p_solve = list(map(calcPlanca, x_solve))
    log_u_p_solve = [m.log(i) for i in u_p_solve]
    log_u_solve = [m.log(i) for i in u_solve]
    fig, ax = plt.subplots(1, 2)
    ax[0].plot(x_solve, log_u_solve, label="u")
    ax[0].plot(x_solve, log_u_p_solve, label="u_p")
    ax[0].legend()
    ax[1].plot(x_solve, f_solve, label="F")
    ax[1].legend()
    plt.show()

shotMethodSolve()