import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from CurvePlot import CurveType, PloterType


def main():
    with h5.File("current_state.h5", "r") as inf:
        Pr = inf["/group"].attrs["Prandtl_number"]
        Ra = inf["/group"].attrs["Rayleigh_number"]
        a = inf["/group"].attrs["wave_number"]
        Nx = inf["/group"].attrs["Nx"]
        Ny = inf["/group"].attrs["Ny"]
        l = inf["/group"].attrs["length"]
        label = r"Pr={0:.1f}, Ra={1:.0e}, a={2:.2f}".format(Pr, Ra, a)
        x = np.linspace(0.0, 1.0, 1001)
        n = np.arange(2.0, 2*(Ny+1), 2.0)
        amps = np.tile(np.asarray(
            inf["/group/states/1"][:Ny], dtype=np.float), (len(x), 1))
        xx, nn = np.meshgrid(x, n, indexing="ij")
    y = 1-x+np.transpose(
        np.sum(amps*np.sin(np.pi*nn*xx), axis=1))
    ploter = PloterType(title=r"Temperature Profile",
                        xlabel=r"\(y\)", ylabel=r"\(T\)")
    curve = CurveType(
        label=r"Pr={0:.1f}, Ra={1:.0e}, a={2:.1f}".format(Pr, Ra, a), x=x, y=y)
    ploter.appendCurve(curve)
    ploter.plotFigure()
    ploter.ax.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))
    ploter.savefig("temperature_profile.pdf")
    plt.close()
    return 0


if __name__ == "__main__":
    main()
