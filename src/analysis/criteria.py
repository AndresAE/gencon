from control import tf
from matplotlib import pyplot as plt
from numpy import abs, arccos, arctan2, array, count_nonzero, cos, delete, imag, linspace, log10, meshgrid, polymul, \
    polyval, rad2deg, real, roots, sqrt, tan
from scipy.signal import tf2zpk


def bode_alt(sys, title=''):
    mag, phase, omega = frequency_response(sys)

    plt.subplot(2, 1, 1)
    plt.plot(omega, mag)
    plt.xscale('log')
    plt.ylabel('Magnitude [dB]')
    plt.grid()
    plt.title('Bode' + title)
    plt.subplot(2, 1, 2)
    plt.plot(omega, phase)
    plt.xscale('log')
    plt.ylabel('Phase [deg]')
    plt.xlabel('frequency [rps]')
    plt.grid()
    plt.show()


def frequency_response(tf, w_in=[0.01, 100], plot=0):
    w = linspace(w_in[0], w_in[1], 10000)
    wi = w*1j
    g = 20 * log10(abs(polyval(tf.num[0][0], wi))) - 20 * log10(abs(polyval(tf.den[0][0], wi)))
    zero, pole, gain = tf2zpk(tf.num[0][0], tf.den[0][0])

    m = count_nonzero(zero == 0)
    n = count_nonzero(pole == 0)
    zero = zero[zero != 0]
    pole = pole[pole != 0]

    real_zero = zero[imag(zero) == 0]
    real_pole = pole[imag(pole) == 0]
    p_real_zero = 0
    p_real_pole = 0
    for zi in real_zero:
        p_real_zero += rad2deg(arctan2(w, -real(zi)))
    for ip in real_pole:
        p_real_pole += rad2deg(arctan2(w, -real(ip)))

    imag_zero = zero[imag(zero) != 0]
    imag_pole = pole[imag(pole) != 0]
    n_pairs_zero = len(imag_zero) / 2
    n_pairs_pole = len(imag_pole) / 2
    p_imag_zero = 0
    p_imag_pole = 0
    for ii in range(0, int(n_pairs_zero)):
        char_equ = real(polymul([1, -imag_zero[0]], [1, -imag_zero[1]]))
        w_n = sqrt(char_equ[2])
        zeta = char_equ[1] / (2 * w_n)
        p_imag_zero += rad2deg(arctan2((2 * zeta * w / w_n), (1 - (w / w_n) ** 2)))
        imag_zero = delete(imag_zero, [0, 1])

    for jj in range(0, int(n_pairs_pole)):
        char_equ = real(polymul([1, -imag_pole[0]], [1, -imag_pole[1]]))
        w_n = sqrt(char_equ[2])
        zeta = char_equ[1] / (2 * w_n)
        p_imag_pole += rad2deg(arctan2((2 * zeta * w / w_n), (1 - (w / w_n) ** 2)))
        imag_pole = delete(imag_pole, [0, 1])

    p = 90*m - 90*n + p_real_zero - p_real_pole + p_imag_zero - p_imag_pole

    if plot:
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(w, g)
        plt.xscale('log')
        plt.grid()
        plt.subplot(2, 1, 2)
        plt.plot(w, p)
        plt.xscale('log')
        plt.grid()
        plt.show()

    return g, p, w


def nichols(sys, title=''):
    mag, phase, omega = frequency_response(sys)

    plt.plot(phase, mag)
    plt.plot([-240, -180, -120, -180, -240], [0, 6, 0, -6, 0], 'y')
    plt.plot([-210, -180, -150, -180, -210], [0, 3, 0, -3, 0], 'r')
    plt.xlabel('Phase [deg]')
    plt.xlabel('Magnitude [dB]')
    plt.title('Stability Margins' + title)
    plt.grid()
    plt.show()


def root_locus_alt(sys):
    tfs = tf(sys)
    r = roots(tfs.den[0][0])
    plt.plot(real(r), imag(r), marker="x", linestyle="None")
    plt.plot([-100, 100], [0, 0], linestyle='--', color='k')
    plt.plot([0, 0], [-100, 100], linestyle='--', color='k')
    plt.plot([0, -100], [0, -100*tan(arccos(0.4))], color='r')
    plt.plot([0, -100], [0, 100 * tan(arccos(0.4))], color='r')
    plt.grid()
    plt.show()

