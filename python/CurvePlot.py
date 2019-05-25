import numpy as np
import matplotlib.pyplot as plt
import h5py


plt.rc('text', usetex=True)
plt.rc('font', family='serif')


class CurveType:
    def __init__(self, label, x, y):
        self.label = label
        self.data = np.column_stack((x, y))

    @property
    def x(self):
        return self.data[:, 0]

    @property
    def y(self):
        return self.data[:, 1]

    def saveFile(self, outf, dset_name):
        dset = outf.create_dataset(dset_name, (len(
            self.x), 2), compression="gzip", dtype=np.float64, data=self.data)
        dset.attrs["label"] = self.label

    def loadFile(self, inf, dset_name):
        self.data = np.asarray(inf[dset_name], dtype=np.float64)
        self.label = inf[dset_name].attrs["label"]


class PloterType:
    def __init__(self, title="Figure", xlabel=r"\(x\)", ylabel=r"\(y\)"):
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cuvs = []
        self.fig = None
        self.ax = None

    def appendCurve(self, curve):
        self.cuvs.append(curve)

    def saveCurvesToFile(self, outf, grp_name):
        grp = outf.require_group(grp_name)
        for i in range(len(self.cuvs)):
            self.cuvs[i].saveFile(outf, str(i))

    def loadCurvesFromFile(self, inf, grp_name):
        grp = inf[grp_name]
        for key in grp.keys:
            cuv = CurveType(grp[key].attrs["label"], np.asarray(
                grp[key][:, 0], dtype=np.float64), np.asarray(grp[key][:, 1], dtype=np.float64))
            self.cuvs.append(cuv)

    def plotFigure(self):
        self.fig, self.ax = plt.subplots()
        for cuv in self.cuvs:
            self.ax.plot(cuv.x, cuv.y,
                         label=cuv.label)
        self.ax.grid()
        self.ax.legend()
        self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.fig.tight_layout()
        plt.close()

    def savefig(self, file_name):
        self.fig.savefig(file_name)
