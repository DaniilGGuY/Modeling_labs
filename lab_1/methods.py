def euler_base(x0, y0, h, n, func):
    res = []
    x, y = x0, y0
    for i in range(n):
        y += h * func(x, y)
        x += h
        res.append(y)
    return res


def euler_task1(x0, y0, dy0, h, n):
    res = []
    x, y, dy = x0, y0, dy0
    for i in range(n):
        try:
            ddy = -(2 * dy + y) / (4 * x)
        except ZeroDivisionError:
            ddy = 1 / 12
        y += h * dy
        dy += h * ddy
        x += h
        res.append(y)
    return res

# alpha = 1 / 2
def runge_kutta4(x0, y0, h, n, func):
    res = []
    x, y = x0, y0
    for i in range(n):
        y += h * (func(x, y) / 2 + func(x + h, y + h * func(x, y)) / 2)
        x += h
        res.append(y)
    return res
