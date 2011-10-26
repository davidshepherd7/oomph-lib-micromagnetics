syms x y t llg_precession_coeff llg_damping_coeff exchange_coeff;

%M = cos(2*pi*x) * [cos(t), 0, 0];
%M = cos(2*pi*x)* [cos(t), sin(t), 0]
%M = cos(t) *[cos(2*pi*x), cos(2*pi*x),0]
%M = [x* cos(t), 0, 0]
%M = [x*cos(t), (1-x)*cos(t), 0]
%M = [0, x*cos(t), 0 ] % doesn't work if we want phi bc non-zero '
%M = [sin(t), cos(t), 0] % again does not work if we want phi bc non-zero
%M = sin(2*pi*x) * [cos(t), sin(t), 0]
%M = (x^2 - x) * [cos(t), sin(t), 0]
%M = cos(t) *[cos(2*pi*x), cos(2*pi*x),0]
M = cos(t) *[cos(2*pi*x)*sin(2*pi*y), cos(2*pi*y)*sin(2*pi*x),0] % have to calculate phi by hand

  ddphi = 4*pi*(diff(M(1),x) + diff(M(2),y));
%dphi = int(ddphi,x);
%phi = int(dphi,x)

%phi = x %% used to set phi when we want non-zero bc
%phi = 2*pi*cos(t)*x*x

  phi = 2*cos(t)*sin(2*pi*x)*sin(2*pi*y) %% phi for M = cos(t) *[cos(2*pi*x)*sin(2*pi*y), cos(2*pi*y)*sin(2*pi*x),0]

dMdt = diff(M,t)

  H_demag = [-diff(phi,x), -diff(phi,y), 0]

  dMx = [diff(M(1),x), diff(M(2),x),0];
dMy = [diff(M(1),y), diff(M(2),y),0];
H_ex = exchange_coeff*[diff(dMx(1),x) + diff(dMy(1),y), diff(dMx(2),x) + diff(dMy(2),y) , 0]

H = H_demag + H_ex
MxH = llg_precession_coeff*cross(M,H)
MxMxH = llg_damping_coeff*cross(M,cross(M,H))

  source = dMdt + MxH + MxMxH;

% check residual is zero:
simple(dMdt + cross(M,H) + cross(M,cross(M,H)) - source) == 0

ccode(source)
