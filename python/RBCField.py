import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import h5py as h5


class FieldType:

    def __init__(self, Prandtl_number=1.0, Rayleigh_number=1000.0, wave_number=1.0):
        self.Pr = Prandtl_number
        self.Ra = Rayleigh_number
        self.a = wave_number
        self.Nx = None
        self.Ny = None
        self.l = None
        self.f = None
        self.fig = None
        self.ax = None
        self._x = None
        self._y = None
        self._xx = None
        self._yy = None
        self._xxx = None
        self._yyy = None
        self._nnn = None
        self._xxxx = None
        self._yyyy = None
        self._mmmm = None
        self._nnnn = None
        self._phi = None
        self._theta = None
        self._psi = None
        self._temperature = None
        self._vorticity = None
        self._u = None
        self._v = None
        self._speed = None
        return

    def loadFile(self, inf, h5grp_name):
        self.Pr = inf[h5grp_name].attrs["Prandtl_number"]
        self.Ra = inf[h5grp_name].attrs["Rayleigh_number"]
        self.a = inf[h5grp_name].attrs["wave_number"]
        self.Nx = inf[h5grp_name].attrs["Nx"]
        self.Ny = inf[h5grp_name].attrs["Ny"]
        self.l = inf[h5grp_name].attrs["length"]
        self.f = np.asarray(
            inf["{0:s}/states/{1:d}".format(h5grp_name, int(self.l-1))], dtype=np.float)
        return

    def initialize(self, x, y):
        self._x = x
        self._y = y
        self._xx, self._yy = np.meshgrid(self._x, self._y)
        self._xxx, self._yyy, self._nnn = np.meshgrid(
            self._x, self._y, np.arange(2.0, 2.0*(self.Ny+1.0), 2.0, dtype=np.float), indexing="ij")
        self._xxxx, self._yyyy, self._mmmm, self._nnnn = np.meshgrid(x, y, np.empty(
            shape=(self.Nx,), dtype=np.float), np.empty(shape=(self.Ny,), dtype=np.float), indexing="ij")
        mm = np.load("mm.npy")[:self.Nx, :self.Ny]
        nn = np.load("nn.npy")[:self.Nx, :self.Ny]
        self._mmmm = np.tile(mm, (len(self._x), len(self._y), 1, 1))
        self._nnnn = np.tile(nn, (len(self._x), len(self._y), 1, 1))
        return

    def evaluate(self):
        amps = np.tile(self.f[:self.Ny],
                       (len(self._x), len(self._y), 1))
        self._phi = np.transpose(
            np.sum(amps*np.sin(np.pi*self._nnn*self._yyy), axis=2))
        amps = np.tile(np.reshape(
            self.f[self.Ny:int((self.Nx+1)*self.Ny)], newshape=(self.Nx, self.Ny)), (len(self._x), len(self._y), 1, 1))
        self._theta = np.transpose(
            np.sum(amps*np.cos(2.0*np.pi*self.a*self._mmmm*self._xxxx)*np.sin(np.pi*self._nnnn*self._yyyy), axis=(2, 3)))
        amps = np.tile(np.reshape(self.f[
            int((self.Nx+1)*self.Ny):], newshape=(self.Nx, self.Ny)), (len(self._x), len(self._y), 1, 1))
        self._psi = np.transpose(np.sum(
            amps*np.sin(2.0*np.pi*self.a*self._mmmm*self._xxxx)*np.sin(np.pi*self._nnnn*self._yyyy), axis=(2, 3)))
        self._vorticity = np.transpose(np.sum(((-4.0*np.pi*np.pi*self.a*self.a)*self._mmmm*self._mmmm+(-np.pi*np.pi)*self._nnnn*self._nnnn)
                                              * amps*np.sin(2.0*np.pi*self.a*self._mmmm*self._xxxx)*np.sin(np.pi*self._nnnn*self._yyyy), axis=(2, 3)))
        self._u = np.transpose(np.sum(-amps*np.pi*self._nnnn*np.sin(2.0*np.pi*self.a *
                                                                    self._mmmm*self._xxxx)*np.cos(np.pi*self._nnnn*self._yyyy), axis=(2, 3)))
        self._v = np.transpose(np.sum(amps*2.0*np.pi*self._mmmm*np.cos(2.0*np.pi*self.a *
                                                                       self._mmmm*self._xxxx)*np.sin(np.pi*self._nnnn*self._yyyy), axis=(2, 3)))
        self._speed = np.sqrt(self._u*self._u + self._v*self._v)
        self._temperature = self.Ra*(1-self._yy)+self._phi+self._theta
        return

    def plotFigures(self, component=None):
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.rc('font', size=16)
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim([0.0, 1.0/self.a])
        self.ax.set_ylim([0.0, 1.0])
        self.ax.set_xlabel(r"\(x\)")
        self.ax.set_ylabel(r"\(y\)")
        self.ax.set_aspect("equal")
        # self.fig.tight_layout()
        if component == "psi":
            self.ax.set_title(r"\(\psi\)")
            levels = MaxNLocator(nbins=25).tick_values(
                self._psi.min(), self._psi.max())
            cntr = self.ax.contourf(
                self._xx, self._yy, self._psi, levels=levels, cmap="binary")
            crb = self.fig.colorbar(cntr, ax=self.ax)
            crb.formatter.set_scientific(True)
            crb.formatter.set_powerlimits((-1, 1))
            crb.update_ticks()
        elif component == "theta":
            self.ax.set_title(r"\(\theta\)")
            levels = MaxNLocator(nbins=25).tick_values(
                self._theta.min(), self._theta.max())
            cntr = self.ax.contourf(
                self._xx, self._yy, self._theta, levels=levels, cmap="jet")
            crb = self.fig.colorbar(cntr, ax=self.ax)
            crb.formatter.set_scientific(True)
            crb.formatter.set_powerlimits((-1, 1))
            crb.update_ticks()
        elif component == "phi":
            self.ax.set_title(r"\(\phi\)")
            levels = MaxNLocator(nbins=25).tick_values(
                self._phi.min(), self._phi.max())
            cntr = self.ax.contourf(
                self._xx, self._yy, self._phi, levels=levels, cmap="jet")
            crb = self.fig.colorbar(cntr, ax=self.ax)
            crb.formatter.set_scientific(True)
            crb.formatter.set_powerlimits((-1, 1))
            crb.update_ticks()
        elif component == "velocity":
            self.ax.set_title(r"Velocity Field")
            cd = 5*self._speed/self._speed.max()
            strm = self.ax.streamplot(self._xx, self._yy, self._u, self._v, density=1.5,
                                      linewidth=1, color=cd, cmap="viridis")
            self.fig.colorbar(strm.lines, ax=self.ax)
        elif component == "vorticity":
            self.ax.set_title(r"Vorticity Field")
            levels = MaxNLocator(nbins=25).tick_values(
                self._vorticity.min(), self._vorticity.max())
            cntr = self.ax.contourf(
                self._xx, self._yy, self._vorticity, levels=levels, cmap="bwr")
            crb = self.fig.colorbar(cntr, ax=self.ax)
            crb.formatter.set_scientific(True)
            crb.formatter.set_powerlimits((-1, 1))
            crb.update_ticks()
        elif component == "temperature":
            self.ax.set_title(r"Temperature Field")
            levels = MaxNLocator(nbins=25).tick_values(
                self._temperature.min(), self._temperature.max())
            cntr = self.ax.contourf(
                self._xx, self._yy, self._temperature, levels=levels, cmap="jet")
            crb = self.fig.colorbar(cntr, ax=self.ax)
            crb.formatter.set_scientific(True)
            crb.formatter.set_powerlimits((-1, 1))
            crb.update_ticks()
        return

    def savefig(self, file_name):
        self.fig.savefig(file_name)
        return


def main():
    x = np.linspace(0.0, 2.0, 201)
    y = np.linspace(0.0, 1.0, 101)
    components = ["psi", "velocity", "vorticity",
                  "temperature"]
    s = FieldType()
    with h5.File("data.h5", "r") as inf:
        s.loadFile(inf, "/group")
    s.initialize(x, y)
    s.evaluate()
    for component in components:
        s.plotFigures("{0:s}".format(component))
        s.savefig("{0:s}.pdf".format(component))
        plt.close()
    return 0


if __name__ == "__main__":
    main()
