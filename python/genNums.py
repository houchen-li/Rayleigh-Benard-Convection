import numpy as np


def genNums(Nx, Ny):
    m = np.arange(1, Nx+1, 1.0, dtype=np.float)
    n = np.arange(1, 2*Ny+1, 2.0, dtype=np.float)
    mm, nn = np.meshgrid(m, n, indexing="ij")
    for i in range(1, Nx, 2):
        nn[i, :] += 1.0
    return mm, nn


def main():
    Nx = 2000
    Ny = 1000
    mm, nn = genNums(Nx, Ny)
    np.save("mm.npy", mm)
    np.save("nn.npy", nn)
    return 0


if __name__ == "__main__":
    main()
