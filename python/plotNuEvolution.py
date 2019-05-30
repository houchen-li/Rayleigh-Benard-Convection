import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from CurvePlot import CurveType, PloterType


def main():
    with open("Nu.dat", "r") as f:
        f.readline()
        Pr = np.float(f.readline().split()[1])
        Ra = np.float(f.readline().split()[1])
        a = np.float(f.readline().split()[1])
    data = np.genfromtxt("Nu.dat", dtype=np.float,
                         skip_header=7, skip_footer=1)
    label = r"Pr={0:.1f}, Ra={1:.0e}, a={2:.1f}".format(Pr, Ra, a)
    curve = CurveType(label, data[:, 0], data[:, 1])
    ploter = PloterType(title=r"Nu(t)", xlabel=r"t", ylabel=r"Nu")
    ploter.appendCurve(curve)
    ploter.plotFigure()
    ploter.ax.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))
    ploter.savefig("Nu.pdf")
    plt.close()
    return 0


if __name__ == "__main__":
    main()
