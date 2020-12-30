from numpy import arctan, random, rad2deg, sqrt


def angle_of_attack(x):
    v = sqrt(x[0]**2+x[1]**2+x[2]**2)
    aoa = rad2deg(arctan(x[2] / v))
    return aoa


def angle_of_sideslip(x):
    aos = rad2deg(arctan(x[1] / x[0]))
    return aos


def unit_noise():
    noise = (1 - 2 * random.rand())
    return noise
