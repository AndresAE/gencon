from matplotlib import pyplot as plt
from numpy import cos, linspace, pi, random, size, zeros
from control import tf, c2d
from scipy.signal import lfilter

# Sample time
dt = 0.01
# Simulation time
t_end = 40
# Time vector
t = linspace(0, t_end, int(1/dt * t_end+1))
# Generate noise for tubulence
noise = 4 * (1-2*random.rand(size(t)))
# Filter turbulence as Dryden's turbulence model
tau_turb = 0.7
lpftf_turb = tf(1, [tau_turb, 1])
lpf_turb = c2d(lpftf_turb, dt, method='matched')
alpha_turb = lfilter(lpf_turb.num[0][0], lpf_turb.den[0][0], noise)
# Low pass filters definitions
tau = linspace(0.1, 1.1, int(1/0.2 * 1.1 + 1))
lpf = dict()
i = 0
for i_tau in tau:
    string = str(i)
    lpf[string] = c2d(tf(1, [i_tau, 1]), dt, 'matched')
    i = i + 1

# Construct the ansgle of attack measured by the Inertial Measurement Unit
alpha_imu = (1-cos(0.4*pi * t))
alpha_imu[t < 25] = 0
alpha_imu[(t > 27.5) & (t < 32.5)] = 2
alpha_imu[t > 35] = 0
alpha_imu = 4 + alpha_imu
# Generate the ADS' angle of attack with turbulence
alpha_turbulence = alpha_imu + alpha_turb
# Perform the SISO filtering for each one of the tau parameters
alpha_filtered = zeros([size(tau), size(t)])
turbulence_filtered = zeros([size(tau), size(t)])

i = 0
for i_tau in tau:
    string = str(i)
    alpha_filtered[i, :] = lfilter(lpf[string].num[0][0], lpf[string].den[0][0], alpha_turbulence)
    turbulence_filtered[i, :] = lfilter(lpf[string].num[0][0], lpf[string].den[0][0], alpha_turb)
    i = i + 1

# Plot the results
# -------------------------------------------------------------------------

plt.figure()
plt.plot(t-20, alpha_imu, 'k', label='alpha')
plt.plot(t-20, alpha_turbulence, 'k', label='alpha_turb')
plt.plot(t-20, alpha_filtered[5, :], 'r', label='alpha filtered')
plt.plot(t-20, alpha_imu + turbulence_filtered[5, :], 'b', label='alpha hybrid')
plt.xlabel('t [sec]')
plt.ylabel('Angle of Attack [deg]')
plt.legend()
plt.xlim([0, t[-1]-20])
plt.grid()
plt.show()
