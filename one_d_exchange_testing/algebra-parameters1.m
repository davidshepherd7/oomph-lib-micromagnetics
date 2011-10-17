syms M x t phi omega;

M = cos(2*pi*x) * [cos(t), 0, 0];

ddphi = 4*pi*diff(M(1),x)
dphi = int(ddphi,x)
phi = int(dphi,x)

dMdt = diff(M,t);

H_demag = [-diff(phi,x), 0, 0];

exchange_coeff = 1;
dM = diff(M(1),x);
H_ex = [exchange_coeff*diff(dM(1),x), 0, 0]

dMdt
H_demag
H = H_demag + H_ex
MxH = cross(M,H)
MxMxH = cross(M,cross(M,H))

C = dMdt + cross(M,H) + cross(M,cross(M,H))

% check residual is zero:
dMdt + cross(M,H) + cross(M,cross(M,H)) - C == 0
