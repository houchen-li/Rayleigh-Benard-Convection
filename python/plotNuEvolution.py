import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from CurvePlot import CurveType, PloterType


def main():
    with h5.File("data.h5", "r") as inf:
        Pr = inf["/group"].attrs["Prandtl_number"]
        Ra = inf["/group"].attrs["Rayleigh_number"]
        a = inf["/group"].attrs["wave_number"]
        label = r"Pr={0:.1f}, Ra={1:.0e}, a={2:.2f}".format(Pr, Ra, a)
        x = np.asarray(inf["/group/t"], dtype=np.float64)
        y = np.asarray(inf["/group/Nu"], dtype=np.float64)
    curve = CurveType(label, x, y)
    ploter = PloterType(title=r"Nu(t)", xlabel=r"t", ylabel=r"Nu")
    ploter.appendCurve(curve)
    ploter.plotFigure()
    ploter.ax.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))
    ploter.savefig("Nu.pdf")
    plt.close()
    return 0


if __name__ == "__main__":
    main()
