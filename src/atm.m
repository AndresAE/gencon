pkg load control
pkg load signal

a = [1,2,3,4;5,6,7,8;4,3,2,1;9,8,7,6];
b = [1;2;3;4];
c = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
d = [0;0;0;0];
sys = ss(a, b, c, d);
[num, den] = ss2tf(sys);