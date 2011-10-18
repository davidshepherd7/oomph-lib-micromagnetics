syms M x t phi omega c d;

% Input the formula for M
M = [ sin(omega*t), cos(omega*t), 0];

% Calculate phi (works for the one D case ONLY)
ddphi = 4*pi*diff(M(1),x);
dphi = int(ddphi,x) + c;
phi = int(dphi,x) + d;

% Calculate other quantities (again for the one D case ONLY)
dMdt = diff(M,t);
H_demag = [-diff(phi,x), 0, 0];

% we have to manually apply boundary conditions - for phi = x we let c=1 d=0.
  subs(c,1);
subs(d,0);

% Output answers
dMdt
H_demag
MxH = cross(M,H_demag)
  MxMxH = cross(M,cross(M,H_demag))
  llg_source = dMdt + cross(M,H_demag) + cross(M,cross(M,H_demag))

% check:
dMdt == -cross(M,H_demag) - cross(M,cross(M,H_demag)) + llg_source
