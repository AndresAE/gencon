from control import tf, tf2ss
from src.common.Atmoshpere import von_karman_disturbance
from src.common.Pilot import Pilot
gain = 1
lead = 0.1
delay = 0.1
alt = 200
v = 200
b = 100
level = 3

phi_2_ail = tf([6.8], [1, 0.44, 0])
pilot = Pilot(gain, delay).delay_lead_lag(lead, 0)
h_u, h_v, h_w, h_p, h_q, h_r = von_karman_disturbance(v, alt, b, level)

sys1 = tf2ss(h_p)
a = 1