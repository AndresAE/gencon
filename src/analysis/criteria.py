from control import bode
from matplotlib import pyplot as plt
from numpy import rad2deg, log10


def nichols(sys, title=''):
    mag, phase, omega = bode(sys, plot=False)

    plt.plot(rad2deg(phase), 20*log10(mag))
    plt.plot([-240, -180, -120, -180, -240], [0, 6, 0, -6, 0], 'y')
    plt.plot([-210, -180, -150, -180, -210], [0, 3, 0, -3, 0], 'r')
    plt.xlim([-315, 0])
    plt.ylim([-50, 50])
    plt.xlabel('Phase [deg]')
    plt.xlabel('Magnitude [dB]')
    plt.title('Stability Margins' + title)
    plt.grid()
    plt.show()
