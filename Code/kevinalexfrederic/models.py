from scipy.integrate import odeint
from numpy import linspace, pi as npi
from sympy import diff, lambdify, symbols, log, sqrt


class SimpleGalaxy:
    _g_const = 4.302e-6  # in kpc . (km/s)^2 . M_sol^-1

    def __init__(self, v, a, b, x0, y0, z0, vx0, vy0, vz0, **kwargs):
        self.u = v**2
        self.a = a
        self.b = b
        self.qinit = [x0, vx0, y0, vy0, z0, vz0]
        self._dvx = None
        self._dvy = None
        self._dvz = None
        self._dvx2 = None
        self._dvy2 = None
        self._dvz2 = None
        self.pot = None

        self.update()

    def update(self):
        x, y, z = symbols('x y z')

        p = self._potexpr(x, y, z)

        self._dvx = lambdify((x, y, z), diff(p, x), modules="numpy")
        self._dvy = lambdify((x, y, z), diff(p, y), modules="numpy")
        self._dvz = lambdify((x, y, z), diff(p, z), modules="numpy")

        self._dvx2 = lambdify((x, y, z), diff(p, x, 2), modules="numpy")
        self._dvy2 = lambdify((x, y, z), diff(p, y, 2), modules="numpy")
        self._dvz2 = lambdify((x, y, z), diff(p, z, 2), modules="numpy")

        self.pot = lambdify((x, y, z), p, modules="numpy")

    def _potexpr(self, x, y, z):
        return (self.u / 2) * log(1 + x**2 + (y / self.a)**2 + (z / self.b)**2)

    def rho(self, x, y, z):
        return (self._dvx2(x, y, z) + self._dvy2(x, y, z) + self._dvz2(x, y, z)) / (4 * self._g_const * npi)

    def deriv(self, q, t):
        nq = [0, 0, 0, 0, 0, 0]  # new array to return
        nq[0] = q[1]
        nq[2] = q[3]
        nq[4] = q[5]

        nq[1] = -self._dvx(q[0], q[2], q[4])
        nq[3] = -self._dvy(q[0], q[2], q[4])
        nq[5] = -self._dvz(q[0], q[2], q[4])

        return nq

    def integrate(self, period, steps):
        time = linspace(start=0, stop=period, num=steps)
        return odeint(self.deriv, self.qinit, time)

    def energy(self, x, y, z, vx, vy, vz):
        return self.pot(x, y, z) + (vx**2 + vy**2 + vz**2) / 2

    def energy0(self):
        return self.energy(self.qinit[0], self.qinit[2], self.qinit[4], self.qinit[1], self.qinit[3], self.qinit[5])


class BlackHoleGalaxy(SimpleGalaxy):
    _g_const = 4.302e-6  # in kpc . (km/s)^2 . M_sol^-1

    def __init__(self, v, a, b, x0, y0, z0, vx0, vy0, vz0, **kwargs):
        self.m = kwargs.pop("m", 0)
        self.e = kwargs.pop("e", 0)
        SimpleGalaxy.__init__(self, v, a, b, x0, y0, z0, vx0, vy0, vz0)

    def _potexpr(self, x, y, z):
        r2 = x**2 + y**2 + z**2
        p = SimpleGalaxy._potexpr(self, x, y, z)
        p += -self._g_const * self.m / sqrt(self.e**2 + r2)
        return p
