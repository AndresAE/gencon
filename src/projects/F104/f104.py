"""Lockheed F104 model."""
from numpy import array, cos, deg2rad, linalg
from common.Atmoshpere import Atmosphere
from common.equations_of_motion import nonlinear_eom, local_acceleration
from common.Gravity import Gravity
from common.rotations import body_to_wind, ned_to_body
from common.tools import angle_of_attack, angle_of_sideslip, dynamic_pressure, flight_path_angle, speed
pilot_station = array([6, 0, 2])
cg = array([25, 0, 2])


def model(x, u_in):
    """nonlinear f104 flight model."""
    m = 13600 / 32.174
    j = array([[3600, 0, 0], [0, 59000, 0], [0, 0, 60000]])
    s = 196
    b = 21.9
    c = 9.6

    c_d_0 = 0.1189
    c_d_a = 0.455
    c_d_q = 0
    c_d_de = 0
    c_l_0 = 0.24
    c_l_a = 3.44
    c_l_q = 2.3
    c_l_de = 0.684
    c_mp_0 = 0.03
    c_mp_a = -0.644
    c_mp_q = -5.84
    c_mp_de = -1.6
    c_y_b = -1.18
    c_y_p = 0
    c_y_r = 0
    c_y_da = 0
    c_y_dr = 0.329
    c_mr_b = -0.175
    c_mr_p = -0.285
    c_mr_r = 0.265
    c_mr_da = 0.0392
    c_mr_dr = 0.0448
    c_my_b = 0.507
    c_my_p = -0.144
    c_my_r = -0.753
    c_my_da = 0.0042
    c_my_dr = -0.1645
    f_x_t = 3500

    u = feedforward_control_law(u_in)

    g = Gravity(x[-1]).gravity()
    atm = Atmosphere(x[-1])
    v = speed(x)
    q = dynamic_pressure(atm.air_density(), x)

    alpha = deg2rad(angle_of_attack(x))
    beta = deg2rad(angle_of_sideslip(x))

    # states ###########################################################################################################
    w = array([0, 0, g * m])
    b_n2b = ned_to_body(x[3], x[4], x[5])
    w_b = b_n2b @ w
    f_x = - q * s * (c_d_0 + c_d_a * alpha + c_d_q * x[7]*c/(2*v) + c_d_de * u[1]) + f_x_t * u[3]
    f_y = q * s * (c_y_b * beta + c_y_p * x[6]*b/(2*v) + c_y_r * x[8]*b/(2*v) + c_y_da * u[0] + c_y_dr * u[2])
    f_z = - q * s * (c_l_0 + c_l_a * alpha + c_l_q * x[7]*c/(2*v) + c_l_de * u[1])
    m_x = q * s * b * (c_mr_b * beta + c_mr_p * x[6]*b/(2*v) + c_mr_r * x[8]*b/(2*v) + c_mr_da * u[0] + c_mr_dr * u[2])
    m_y = q * s * c * (c_mp_0 + c_mp_a * alpha + c_mp_q * x[7]*c/(2*v) + c_mp_de * u[1])
    m_z = q * s * b * (c_my_b * beta + c_my_p * x[6]*b/(2*v) + c_my_r * x[8]*b/(2*v) + c_my_da * u[0] + c_my_dr * u[2])
    fm = array([f_x, f_y, f_z, m_x, m_y, m_z])
    b_b2w = body_to_wind(alpha, beta)
    fm[0:3] = linalg.inv(b_b2w) @ fm[0:3]
    fm[3:6] = linalg.inv(b_b2w) @ fm[3:6]
    fm[0:3] = fm[0:3] + w_b

    # state derivative #################################################################################################
    x_dot = nonlinear_eom(x, m, j, fm)

    return x_dot


def feedforward_control_law(u, on=True):
    if on:
        u_out = u
    else:
        u_out = u
    return u_out


def feedback_control_law(x, y, u, w, on=True):
    """f104 control law."""
    k_p_da = 0
    k_psi_da = 0

    k_q_de = 0
    k_nz_de = 0  # 0
    k_alt_de = 0  # 0.001

    k_r_dr = 0
    if on:
        u[0] = u[0] + k_p_da * x[6] + k_psi_da * (x[5] - deg2rad(3))
        u[2] = u[2] + k_r_dr * x[8]
        n_z = 1 / cos(x[3])
        if x[-1] < 30:
            gamma_t = -3 / 50 * x[-1]
            gamma_i = flight_path_angle(x)
            cmd_gamma = -0 * (gamma_t - gamma_i)
        else:
            cmd_gamma = 0
        u[1] = u[1] + k_q_de * x[7] + k_nz_de * (y[1] - n_z) + k_alt_de * (x[-1] - 500) + cmd_gamma
    return u


def outputs(x_dot, x, u):
    g = Gravity(x[-1]).gravity()
    alpha = deg2rad(angle_of_attack(x))
    n_x = x_dot[0] / g
    n_y = x_dot[1] / g
    n_z = 1 + x_dot[2] / g
    a_p = local_acceleration(pilot_station, cg, x, x_dot)
    n_z_p = 1 + a_p[2] / g
    n_x_p = a_p[0] / g
    n_y_p = a_p[1] / g
    gamma = x[4] - alpha
    y = array([gamma, n_z, n_x, n_y, n_x_p, n_y_p, n_z_p])
    return y
