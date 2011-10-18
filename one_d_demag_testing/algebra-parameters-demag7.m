syms M x t phi omega;

M = (x*x - x) * [ cos(omega*t), sin(omega*t), 0];

ddphi = 4*pi*diff(M(1),x)
dphi = int(ddphi,x)
phi = int(dphi,x)

dMdt = diff(M,t);

H_demag = [-diff(phi,x), 0, 0];

dMdt
H_demag
MxH = cross(M,H_demag)
MxMxH = cross(M,cross(M,H_demag))

C = dMdt + cross(M,H_demag) + cross(M,cross(M,H_demag))

% check residual is zero:
dMdt + cross(M,H_demag) + cross(M,cross(M,H_demag)) - C == 0
