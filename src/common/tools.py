from control import tf
from numpy import arctan, random, rad2deg, sqrt


def angle_of_attack(x):
    v = sqrt(x[0]**2+x[1]**2+x[2]**2)
    aoa = rad2deg(arctan(x[2] / v))
    return aoa


def angle_of_sideslip(x):
    aos = rad2deg(arctan(x[1] / x[0]))
    return aos


def dynamic_pressure(rho, x):
    v = speed(x)
    q = 0.5 * rho * v ** 2
    return q


def lead_lag(t_lead, t_lag):
    comp = tf([t_lead, 1], [t_lag, 1])
    return comp


def pade_1(tau):
    approx = tf([-tau/2, 1], [tau/2, 1])
    return approx


def speed(x):
    v = sqrt(x[0]**2+x[1]**2+x[2]**2)
    return v


def unit_noise():
    noise = (1 - 2 * random.rand())
    return noise
