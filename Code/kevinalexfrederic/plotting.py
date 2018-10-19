import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import linspace, meshgrid
import math as mt
# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D


class ModelPlotter:
    def __init__(self, model):
        self.model = model
        self._q = []
        self.fig = plt.figure(figsize=plt.figaspect(1)*3)
        self._u = []
        self._v = []
        self._minval = 0,
        self._maxval = 0

    def prepare(self, period, steps):
        self._q = self.model.integrate(period, steps)
        self._minval = min(min(self._q[:, 0]), min(self._q[:, 2]), min(self._q[:, 4]))
        self._maxval = max(max(self._q[:, 0]), max(self._q[:, 2]), max(self._q[:, 4]))

        self._minval *= 1.4  # add some extra space to the axes
        self._maxval *= 1.4

        [self._u, self._v] = meshgrid(linspace(self._minval, self._maxval, num=800),
                                      linspace(self._minval, self._maxval, num=800))

    def show(self):
        plt.show()

    def save(self, name):
        plt.savefig(name)

    def close(self):
        plt.close()


class OrbitPlot(ModelPlotter):
    def __init__(self, model, **kwargs):
        ModelPlotter.__init__(self, model)

    def plot(self, title, period, steps):
        self.prepare(period, steps)

        print "\tXY-vlak"
        self.subplot(221, "XY")
        print "\tYZ-vlak"
        self.subplot(222, "YZ")
        print "\tXZ-vlak"
        self.subplot(223, "XZ")
        print "\t3D plot"
        self.subplot3d(224)

        plt.suptitle(title, fontsize=34)
        plt.tick_params(labelsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2)

    def subplot(self, which, axes):
        plt.subplot(which)
        ax = self.fig.gca()
        ax.set_aspect('equal')

        if axes == "YZ":
            plt.xlabel('y [kpc]', fontsize=18)
            plt.ylabel('z [kpc]', fontsize=18)
            x = 0
            y = self._u
            z = self._v
        elif axes == "XZ":
            plt.xlabel('x [kpc]', fontsize=18)
            plt.ylabel('z [kpc]', fontsize=18)
            x = self._u
            y = 0
            z = self._v
        else:
            plt.xlabel('x [kpc]', fontsize=18)
            plt.ylabel('y [kpc]', fontsize=18)
            x = self._u
            y = self._v
            z = 0

        t = self.model.pot(x, y, z)

        plt.contourf(self._u, self._v, t, 24, cmap=cm.Spectral)
        plt.contour(self._u, self._v, t, 8, colors='w', zorder=10, alpha=1)

        e0 = self.model.energy0()
        e = self.model.energy(x, y, z, 0, 0, 0) - e0

        if axes == "YZ":
            plt.plot(self._q[:, 2], self._q[:, 4], color='k', alpha=0.8)
        elif axes == "XZ":
            plt.plot(self._q[:, 0], self._q[:, 4], color='k', alpha=0.8)
        else:
            plt.plot(self._q[:, 0], self._q[:, 2], color='k', alpha=0.8)

        plt.contour(self._u, self._v, e, [0], colors='k', linewidths='3', zorder=15, alpha=0.8)

        plt.grid(True)
        plt.title(axes + "-vlak", fontsize=26)

    def subplot3d(self, which):
        plt.subplot(which, projection='3d')
        ax = self.fig.gca()
        ax.set_aspect('equal')
        ax.auto_scale_xyz([self._minval, self._maxval], [self._minval, self._maxval], [self._minval, self._maxval])
        ax.set_xlim(self._minval, self._maxval)
        ax.set_ylim(self._minval, self._maxval)
        ax.set_zlim(self._minval, self._maxval)
        ax.set_xlabel('x [kpc]', fontsize=18)
        ax.set_ylabel('y [kpc]', fontsize=18)
        ax.set_zlabel('z [kpc]', fontsize=18)

        e0 = self.model.energy0()
        e1 = self.model.energy(self._u, self._v, 0, 0, 0, 0) - e0
        e2 = self.model.energy(self._u, 0, self._v, 0, 0, 0) - e0
        e3 = self.model.energy(0, self._u, self._v, 0, 0, 0) - e0

        plt.contour(self._u, self._v, e1, [0], zdir='z', offset=self._minval, colors='k', linewidths='2', alpha=0.8)
        plt.contour(self._u, e2, self._v, [0], zdir='y', offset=self._maxval, colors='k', linewidths='2', alpha=0.8)
        plt.contour(e3, self._u, self._v, [0], zdir='x', offset=self._minval, colors='k', linewidths='2', alpha=0.8)

        plt.plot(self._q[:, 2], self._q[:, 4], zs=self._minval, zdir='x', lw=.5, c='.5')
        plt.plot(self._q[:, 0], self._q[:, 4], zs=self._maxval, zdir='y', lw=.5, c='.5')
        plt.plot(self._q[:, 0], self._q[:, 2], zs=self._minval, zdir='z', lw=.5, c='.5')
        plt.plot(self._q[:, 0], self._q[:, 2], self._q[:, 4], lw=.5, c='b')

        plt.title("3D plot", fontsize=26)
        ax.patch.set_visible(False)


class DensityPlot(ModelPlotter):
    def __init__(self, model, **kwargs):
        ModelPlotter.__init__(self, model)

    def plot(self, title, period, steps):
        self.prepare(period, steps)
        ax = self.fig.gca()
        ax.set_aspect('equal')

        plt.xlabel('x [kpc]', fontsize=26)
        plt.ylabel('y [kpc]', fontsize=26)

        m = self.model.rho(self._u, self._v, 0)

        plt.contourf(self._u, self._v, m, 24, cmap=cm.Spectral)
        plt.contour(self._u, self._v, m, 8, colors='w', zorder=10, alpha=1)

        plt.plot(self._q[:, 0], self._q[:, 2], color='k', alpha=0.8)

        plt.grid(True)
        plt.suptitle(title, fontsize=34)
        plt.tick_params(labelsize=22)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])


class BlackHolePlot(ModelPlotter):
    def __init__(self, model, **kwargs):
        ModelPlotter.__init__(self, model)
        self.masses = kwargs.pop("masses", [0])
        self.period = 0
        self.steps = 0
        self.columns = mt.floor(mt.sqrt(len(self.masses)))
        self.rows = round(len(self.masses) / self.columns)
        self.fig.set_size_inches(plt.figaspect(self.rows/self.columns)*3)

    def plot(self, title, period, steps):
        self.period = period
        self.steps = steps
        i = 1
        for mass in self.masses:
            print "\tMassa " + str(i) + " van " + str(len(self.masses)) + ": " + "%.2e M_sol" % mass
            which = 100 * self.rows + 10 * self.columns + i
            self.massplot(which, mass)
            i += 1

        plt.suptitle(title, fontsize=34)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    def massplot(self, which, mass):
        self.model.m = mass
        self.model.update()

        self.prepare(self.period, self.steps)
        plt.subplot(which)
        ax = self.fig.gca()
        ax.set_aspect('equal')

        plt.xlabel('y [kpc]', fontsize=18)
        plt.ylabel('x [kpc]', fontsize=18)

        t = self.model.pot(self._u, self._v, 0)

        plt.contourf(self._u, self._v, t, 24, cmap=cm.Spectral)
        plt.contour(self._u, self._v, t, 8, colors='w', zorder=10, alpha=1)

        e0 = self.model.energy0()
        e = self.model.energy(self._u, self._v, 0, 0, 0, 0) - e0

        plt.plot(self._q[:, 0], self._q[:, 2], color='k', alpha=0.8)

        plt.contour(self._u, self._v, e, [0], colors='k', linewidths='3', zorder=15, alpha=0.8)

        plt.grid(True)
        plt.title("BH mass = " + '%.1e' % mass + " $M_{\odot}$", fontsize=26)
