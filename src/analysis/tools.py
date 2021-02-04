from common.Atmoshpere import von_karman_disturbance
from common.tools import unit_noise
from control import forced_response
from matplotlib import pyplot as plt
from numpy import append, identity, size, shape, reshape, transpose, zeros


def nonlinear_eom_to_ss(model, outputs, x_0, u_0, y_0, dx=0.1, du=0.1):
    """aircraft system linearization routine."""
    """return jacobians a, b wrt to x_ss and output matrices c, and d wrt u_ss."""
    x = x_0
    u = u_0
    a = zeros((len(x_0), len(x_0)))
    b = zeros((len(x_0), len(u_0)))
    c = zeros((len(y_0), len(x_0)))
    d = zeros((len(y_0), len(u_0)))
    for ii in range(0, len(x_0)):
        x[ii] = x[ii] + dx
        dxdt_1 = model(x, u_0)
        x[ii] = x[ii] - dx

        x[ii] = x[ii] - dx
        dxdt_2 = model(x, u_0)
        ddx_dx = (dxdt_1 - dxdt_2)/(2*dx)
        a[:, ii] = transpose(ddx_dx)
        x[ii] = x[ii] + dx

    for ii in range(0, len(u_0)):
        u[ii] = u[ii] + du
        dxdt_1 = model(x_0, u)
        u[ii] = u[ii] - du

        u[ii] = u[ii] - du
        dxdt_2 = model(x_0, u)
        ddx_dx = (dxdt_1 - dxdt_2)/(2*du)
        b[:, ii] = transpose(ddx_dx)
        u[ii] = u[ii] + du

    for ii in range(0, len(x_0)):
        x[ii] = x[ii] + dx
        dxdt_1 = model(x, u_0)
        dydt_1 = outputs(dxdt_1, x, u_0)
        x[ii] = x[ii] - dx

        x[ii] = x[ii] - dx
        dxdt_2 = model(x, u_0)
        dydt_2 = outputs(dxdt_2, x, u_0)
        ddy_dx = (dydt_1 - dydt_2)/(2*dx)
        c[:, ii] = transpose(ddy_dx)
        x[ii] = x[ii] + dx

    for ii in range(0, len(u_0)):
        u[ii] = u[ii] + du
        dxdt_1 = model(x, u_0)
        dydt_1 = outputs(dxdt_1, x_0, u)
        u[ii] = u[ii] - du

        u[ii] = u[ii] - du
        dxdt_2 = model(x, u_0)
        dydt_2 = outputs(dxdt_2, x_0, u)
        ddy_dx = (dydt_1 - dydt_2)/(2*du)
        d[:, ii] = transpose(ddy_dx)
        u[ii] = u[ii] + du

    a_i = identity(len(x_0))
    b_i = zeros((len(x_0), len(u_0)))
    c = append(a_i, c, axis=0)
    d = append(b_i, d, axis=0)
    return a, b, c, d


def ss_selector(a, b, c, d, x, u, y):
    a_out = a[x, :]
    a_out = a_out[:, x]
    b_out = b[x, :]
    b_out = b_out[:, u]
    s = shape(a_out)
    if not y:
        c_out = c[x, :]
        c_out = c_out[:, x]
        d_out = d[x, :]
        d_out = d_out[:, u]
        b_out = reshape(b_out, (size(x), size(u)))
        d_out = reshape(d_out, (size(x), size(u)))
    else:
        c_out = c[append(x, s[0] + y), :]
        c_out = c_out[:, x]
        d_out = d[append(x, s[0] + y), :]
        d_out = d_out[:, u]
        b_out = reshape(b_out, (size(x), size(u)))
        d_out = reshape(d_out, (size(x) + size(y), size(u)))
    return a_out, b_out, c_out, d_out


def turbulence(v, altitude, b, intensity, t, iplot=0):
    noise = []
    for t_i in t:
        noise.append(unit_noise())

    h_u, h_v, h_w, h_p, h_q, h_r = von_karman_disturbance(v, altitude, b, intensity)
    t, u, x = forced_response(h_u, t, noise)
    t, v, x = forced_response(h_v, t, noise)
    t, w, x = forced_response(h_w, t, noise)
    t, p, x = forced_response(h_p, t, noise)
    t, q, x = forced_response(h_q, t, noise)
    t, r, x = forced_response(h_r, t, noise)

    if iplot:
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(t, u)
        plt.plot(t, v)
        plt.plot(t, w)
        plt.xlabel('linear rates [fps]')
        plt.grid()
        plt.subplot(2, 1, 2)
        plt.plot(t, p)
        plt.plot(t, q)
        plt.plot(t, r)
        plt.grid()
        plt.xlabel('angular rates [rps]')
        plt.xlabel('time [sec]')
    if intensity == 0:
        u = u * 0
        v = v * 0
        w = w * 0
        p = p * 0
        q = q * 0
        r = r * 0
    return u, v, w, p, q, r
