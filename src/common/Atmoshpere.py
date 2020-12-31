"""1976 standard atmosphere model of earth."""
from control import series, tf
from numpy import array, exp, pi, sqrt, interp
from scipy.interpolate import RectBivariateSpline
gas_constant = 1716.5  # ft*lb/(slug*R)
gamma = 1.4  # []
intensity = [1, 2, 3, 4, 5]
turb_altitude = [0, 5000, 25000, 35000, 45000, 65000, 80000]  # [ft]
turb_intensity = array([[5, 2.5, 0, 0, 0, 0, 0],
                        [6, 7, 3, 1, 0, 0, 0],
                        [8, 11, 7, 5, 4, 0, 0],
                        [12, 16, 10, 8, 8, 5, 2.5],
                        [15, 23, 20, 15, 14, 8, 5]])


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


def von_karman_disturbance(v, alt, b, level):
    w_20 = interp(level, [1, 3, 5], [15, 30, 45])

    f = RectBivariateSpline(intensity, turb_altitude, turb_intensity, kx=1)
    sigma_high = f(level, alt)
    sigma_high = float(sigma_high)

    sigma_u_low = von_karman_sigma_u(alt, w_20)
    sigma_v_low = von_karman_sigma_v(alt, w_20)
    sigma_w_low = von_karman_sigma_w(w_20)

    if alt < 1000:
        sigma_u = sigma_u_low
        sigma_v = sigma_v_low
        sigma_w = sigma_w_low
    elif (alt <= 2000) and (alt >= 1000):
        sigma_u = interp(alt, [1000, 2000], [sigma_u_low, sigma_high])
        sigma_v = interp(alt, [1000, 2000], [sigma_v_low, sigma_high])
        sigma_w = interp(alt, [1000, 2000], [sigma_w_low, sigma_high])
    else:
        sigma_u = sigma_high
        sigma_v = sigma_high
        sigma_w = sigma_high

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
    h_p = tf([sigma_w*sqrt(0.8/v)*((pi/(4*b))**(1/6))], ((2*l_w)**(1/3))*array([4*b/(pi*v), 1]))
    return h_p


def von_karman_q(v, b, h_w):
    h_q = series(tf([1/v, 0], [4*b/(pi*v), 1]), h_w)
    return h_q


def von_karman_r(v, b, h_v):
    h_r = series(tf([1/v, 0], [3*b/(pi*v), 1]), h_v)
    return h_r


def von_karman_u(v, sigma_u, l_u):
    h_u = tf(sigma_u*sqrt(2*l_u/(pi*v))*array([0.25*l_u/v, 1]), [0.1987*(l_u/v)**2, 1.357*l_u/v, 1])
    return h_u


def von_karman_v(v, sigma_v, l_v):
    h_v = tf(sigma_v*sqrt(2*l_v/(pi*v))*array([0.3398*(2*l_v/v)**2, 2.7478*(2*l_v/v), 1]),
             [0.1539*(2*l_v/v)**3, 1.9754*(2*l_v/v)**2, 2.9958*(2*l_v/v), 1])
    return h_v


def von_karman_w(v, sigma_w, l_w):
    h_w = tf(sigma_w * sqrt(2 * l_w / (pi * v)) *
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
    l_u_low = h / (0.177 + 0.000823 * h) ** 1.2
    l_u_high = 2500
    if h < 1000:
        l_u = l_u_low
    elif h < 2000:
        l_u = l_u_high
    else:
        l_u = interp(h, [1000, 2000], [l_u_low, l_u_high])

    return l_u


def von_karman_l_v(h):
    l_v_low = h / (2 * (0.177 + 0.000823 * h) ** 1.2)
    l_v_high = 2500 / 2
    if h < 1000:
        l_v = l_v_low
    elif h < 2000:
        l_v = l_v_high
    else:
        l_v = interp(h, [1000, 2000], [l_v_low, l_v_high])
    return l_v


def von_karman_l_w(h):
    l_w_low = h / 2
    l_w_high = 2500 / 2
    if h < 1000:
        l_w = l_w_low
    elif h < 2000:
        l_w = l_w_high
    else:
        l_w = interp(h, [1000, 2000], [l_w_low, l_w_high])
    return l_w
