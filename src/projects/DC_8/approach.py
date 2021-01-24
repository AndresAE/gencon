from control import tf
from numpy import append, array, polymul
from src.analysis.criteria import nichols, bode_alt

den_log = polymul(array([1/(1.619**2), 2*0.522/1.619, 1]), array([1/(0.1635**2), 2*0.0606/0.1635, 1]))
num_u = polymul(348*(-0.641)*array([1/1.08, 1]), array([1/-35.3, 1]))
num_theta = polymul(-0.618*(-1.338)*array([1/0.0605, 1]), array([1/0.535, 1]))
num_q = append(num_theta, 0)

ele2theta = tf(num_theta, den_log)
bode_alt(ele2theta, 'Elevator to q')
