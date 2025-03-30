def rungeKuttaSystem(x0, u0, v0, h, n, f, phi):
    res_x, res_u, res_v = [], [], []
    x, u, v = x0, u0, v0
    for i in range(n):
        k1 = f(x, u, v)
        q1 = phi(x, u, v)
        k2 = f(x + h / 2, u + h / 2 * k1, v + h / 2 * q1)
        q2 = phi(x + h / 2, u + h / 2 * k1, v + h / 2 * q1)
        k3 = f(x + h / 2, u + h / 2 * k2, v + h / 2 *  q2)
        q3 = phi(x + h / 2, u + h / 2 * k2, v + h / 2 *  q2)
        k4 = f(x + h, u + h * k3, v + h * q3)
        q4 = phi(x + h, u + h * k3, v + h * q3)

        x += h
        u = u + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        v = v + h / 6 * (q1 + 2 * q2 + 2 * q3 + q4)
        res_u.append(u)
        res_v.append(v)
        res_x.append(x)

    return res_x, res_u, res_v