import numpy as np
from geofound.fadums_chart import calc_fadums_from_m_and_n



def create():
    import matplotlib.pyplot as plt

    ms1 = np.array([0, 0.01, 0.02, 0.04, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 10])
    ms = ms1[:, np.newaxis] * np.ones((len(ms1), len(ms1)))
    ns = ms.T

    i = calc_fadums_from_m_and_n(ms, ns)
    plt.semilogx(ms1, i)
    plt.show()


if __name__ == '__main__':
    create()

