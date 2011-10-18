syms M x t phi omega;

M = cos(omega*t) * [x, (1-x), 0];

ddphi = 4*pi*diff(M(1),x);
dphi = int(ddphi,x);
phi = int(dphi,x);

dMdt = diff(M,t);

H_demag = [-diff(phi,x), 0, 0];

dMdt
H_demag
MxH = cross(M,H_demag)
MxMxH = cross(M,cross(M,H_demag))

C = dMdt + cross(M,H_demag) + cross(M,cross(M,H_demag))

% check:
dMdt == -cross(M,H_demag) - cross(M,cross(M,H_demag)) + C
