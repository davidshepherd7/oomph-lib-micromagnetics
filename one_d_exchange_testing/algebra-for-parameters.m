syms x t llg_precession_coeff llg_damping_coeff exchange_coeff;

%M = cos(2*pi*x) * [cos(t), 0, 0];
%M = cos(2*pi*x)* [cos(t), sin(t), 0]
%M = cos(t) *[cos(2*pi*x), cos(2*pi*x),0]
%M = [x* cos(t), 0, 0]
%M = [x*cos(t), (1-x)*cos(t), 0]
%M = [0, x*cos(t), 0 ] % doesn't work if we want phi bc non-zero '
%M = [sin(t), cos(t), 0] % again does not work if we want phi bc non-zero
M = sin(2*pi*x) * [cos(t), sin(t), 0]
%M = (x^2 - x) * [cos(t), sin(t), 0]
%M = cos(t) *[cos(2*pi*x), cos(2*pi*x),0]

  ddphi = 4*pi*diff(M(1),x);
dphi = int(ddphi,x);
phi = int(dphi,x)

%phi = x %% used to set phi when we want non-zero bc
%phi = 2*pi*cos(t)*x*x

dMdt = diff(M,t)

H_demag = [-diff(phi,x), 0, 0]

dM = diff(M(1),x);
H_ex = [exchange_coeff*diff(dM(1),x), 0, 0]

H = H_demag + H_ex
MxH = llg_precession_coeff*cross(M,H)
MxMxH = llg_damping_coeff*cross(M,cross(M,H))

  source = dMdt + MxH + MxMxH;

% check residual is zero:
dMdt + cross(M,H) + cross(M,cross(M,H)) - source == 0

ccode(source)
