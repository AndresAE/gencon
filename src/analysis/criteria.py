from control import bode
from matplotlib import pyplot as plt
from numpy import rad2deg, log10


def bode_alt(sys, title=''):
    mag, phase, omega = bode(sys, plot=False)

    plt.subplot(2, 1, 1)
    plt.plot(omega, 20 * log10(mag))
    plt.xscale('log')
    plt.ylabel('Magnitude [dB]')
    plt.grid()
    plt.title('Bode' + title)
    plt.subplot(2, 1, 2)
    plt.plot(omega, rad2deg(phase))
    plt.xscale('log')
    plt.ylabel('Phase [deg]')
    plt.xlabel('frequency [rps]')
    plt.grid()
    plt.show()


def nichols(sys, title=''):
    mag, phase, omega = bode(sys, plot=False)

    plt.plot(rad2deg(phase), 20*log10(mag))
    plt.plot([-240, -180, -120, -180, -240], [0, 6, 0, -6, 0], 'y')
    plt.plot([-210, -180, -150, -180, -210], [0, 3, 0, -3, 0], 'r')
    plt.xlabel('Phase [deg]')
    plt.xlabel('Magnitude [dB]')
    plt.title('Stability Margins' + title)
    plt.grid()
    plt.show()
