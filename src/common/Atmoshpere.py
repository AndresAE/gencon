"""1976 standard atmosphere model of earth."""
from control import series
from scipy.signal import TransferFunction
from numpy import array, exp, pi, polymul, sqrt
gas_constant = 1716.5  # ft*lb/(slug*R)
gamma = 1.4  # []


class Atmosphere:
    def __init__(self, altitude):
        self.altitude = altitude

    def air_density(self):
        """return air density wrt altitude."""
        p = self.pressure()  # psf
        t = temperature(self.altitude)  # R
        rho = p / (gas_constant * t)  # lb/ft3
        return rho

    def pressure(self):
        """return ambient pressure wrt altitude."""
        t = temperature(self.altitude)  # R
        if self.altitude <= 36152:
            p = 2116 * (t / 518.6) ** 5.256  # psf
        else:
            p = 473.1 * exp(1.73 - 0.000048 * self.altitude)  # psf
        return p

    def speed_of_sound(self):
        """return speed of sound wrt altitude."""
        t = temperature(self.altitude)  # R
        a = sqrt(gamma * gas_constant * t)  # [ft/s]
        return a

    def viscosity(self):
        """return dynamic viscosity of air wrt altitude."""
        t_ref = temperature(0)  # R
        t = temperature(self.altitude)  # R
        s = 198.72  # R
        mu_ref = 3.737 * 10 ** (-7)  # [slug/(ft*s)]
        mu = mu_ref * ((t / t_ref) ** (3 / 2)) * (t_ref + s) / (t + s)  # [slug/(ft*s)]
        return mu


# Public Methods #######################################################################################################
def temperature(altitude):
    """return ambient temperature wrt altitude."""
    if altitude <= 36152:
        t = 59 - 0.00356 * altitude  # deg F
    else:
        t = -70  # deg F
    t = t + 459.7  # R
    return t


def von_karman_disturbance(v, alt, b):
    w_20 = 15
    sigma_u = von_karman_sigma_u(alt, w_20)
    sigma_v = von_karman_sigma_v(alt, w_20)
    sigma_w = von_karman_sigma_w(w_20)
    l_u = von_karman_l_u(alt)
    l_v = von_karman_l_v(alt)
    l_w = von_karman_l_w(alt)
    h_u = von_karman_u(v, sigma_u, l_u)
    h_v = von_karman_v(v, sigma_v, l_v)
    h_w = von_karman_w(v, sigma_w, l_w)
    h_p = von_karman_p(v, b, sigma_w, l_w)
    h_q = von_karman_q(v, b, h_w)
    h_r = von_karman_r(v, b, h_v)

    return h_u, h_v, h_w, h_p, h_q, h_r


def von_karman_p(v, b, sigma_w, l_w):
    h_p = TransferFunction([sigma_w*sqrt(0.8/v)*((pi/(4*b))**(1/6))], ((2*l_w)**(1/3))*array([4*b/(pi*v), 1]))
    return h_p


def von_karman_q(v, b, h_w):
    h_q = TransferFunction(polymul(array([1/v, 0]), h_w.num), polymul(array([4*b/(pi*v), 1]), h_w.den))
    return h_q


def von_karman_r(v, b, h_v):
    h_r = TransferFunction(polymul([1/v, 0], h_v.num), polymul([3*b/(pi*v), 1], h_v.den))
    return h_r


def von_karman_u(v, sigma_u, l_u):
    h_u = TransferFunction(sigma_u*sqrt(2*l_u/(pi*v))*array([0.25*l_u/v, 1]), [0.1987*(l_u/v)**2, 1.357*l_u/v, 1])
    return h_u


def von_karman_v(v, sigma_v, l_v):
    h_v = TransferFunction(sigma_v*sqrt(2*l_v/(pi*v))*array([0.3398*(2*l_v/v)**2, 2.7478*(2*l_v/v), 1]),
                           [0.1539*(2*l_v/v)**3, 1.9754*(2*l_v/v)**2, 2.9958*(2*l_v/v), 1])
    return h_v


def von_karman_w(v, sigma_w, l_w):
    h_w = TransferFunction(sigma_w * sqrt(2 * l_w / (pi * v)) *
                           array([0.3398 * (2 * l_w / v) ** 2, 2.7478 * (2 * l_w / v), 1]),
                           [0.1539 * (2 * l_w / v) ** 3, 1.9754 * (2 * l_w / v) ** 2, 2.9958 * (2 * l_w / v), 1])
    return h_w


def von_karman_sigma_u(h, w_20):
    sigma_w = von_karman_sigma_w(w_20)
    sigma_u = sigma_w / (0.177 + 0.000823 * h) ** 0.4
    return sigma_u


def von_karman_sigma_v(h, w_20):
    sigma_w = von_karman_sigma_w(w_20)
    sigma_v = sigma_w / (0.177 + 0.000823 * h) ** 0.4
    return sigma_v


def von_karman_sigma_w(w_20):
    sigma_w = 0.1 * w_20
    return sigma_w


def von_karman_l_u(h):
    l_u = h / (0.177 + 0.000823 * h) ** 1.2
    return l_u


def von_karman_l_v(h):
    l_v = h / (2 * (0.177 + 0.000823 * h) ** 1.2)
    return l_v


def von_karman_l_w(h):
    l_w = h / 2
    return l_w
