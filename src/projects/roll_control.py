from control import bode, feedback, nichols, step_response, tf, tf2ss
from matplotlib import pyplot as plt
from common.Atmoshpere import von_karman_disturbance
from common.Pilot import Pilot

gain = 1
lead = 0
delay = 0.05
alt = 200
v = 200
b = 100
level = 3

phi_2_ail = tf([6.8], [1, 0.44, 0])
pilot = Pilot(gain, delay).delay
h_u, h_v, h_w, h_p, h_q, h_r = von_karman_disturbance(v, alt, b, level)

sys1 = tf2ss(h_p)
sys2 = feedback(phi_2_ail, pilot, -1)
t, y = step_response(sys2)
nichols(sys2)

plt.plot(t, y)
plt.show()
