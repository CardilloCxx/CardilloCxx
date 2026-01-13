import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    # TODO: This needs to be adapted!
    file = "/home/jonas/gitprojects/CardilloMPI/vtk_out/wilberforce_tracked.csv"
    usecols = (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    data = np.loadtxt(file, delimiter=",", skiprows=1, usecols=usecols)
    print(f"data.shape: {data.shape}")
    t, x, y, z, ux, uy, uz, omegax, omegay, omegaz = data.T

    dt = t[1] - t[0]
    alpha = np.cumsum(omegax) * dt * 180 / np.pi
    beta = np.cumsum(omegay) * dt * 180 / np.pi
    gamma = np.cumsum(omegaz) * dt * 180 / np.pi

    fig, ax = plt.subplots(3, 1)

    ax[0].plot(t, x, label="x")
    ax[0].plot(t, y, label="y")
    ax[0].plot(t, z, label="z")
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(t, alpha, label="alpha")
    ax[1].plot(t, beta, label="beta")
    ax[1].plot(t, gamma, label="gamma")
    ax[1].grid()
    ax[1].legend()

    ax[2].plot(t, omegax, label="omega x")
    ax[2].plot(t, omegay, label="omega y")
    ax[2].plot(t, omegaz, label="omega z")
    ax[2].grid()
    ax[2].legend()

    plt.show()