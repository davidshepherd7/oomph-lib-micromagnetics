syms x y z t llg_precession_coeff llg_damping_coeff exchange_coeff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M = cos(2*pi*x) * [cos(t), 0, 0];
%M = cos(2*pi*x)* [cos(t), sin(t), 0]
%M = cos(t) *[cos(2*pi*x), cos(2*pi*x),0]
%M = [x* cos(t), 0, 0]
%M = [x*cos(t), (1-x)*cos(t), 0]
M = [0, x*cos(t), 0 ]
%M = [sin(t), cos(t), 0] % again does not work if we want phi bc non-zero
%M = sin(2*pi*x) * [cos(t), sin(t), 0]
%M = (x^2 - x) * [cos(t), sin(t), 0]
%M = cos(t) *[cos(2*pi*x), cos(2*pi*x),0]
%M = cos(t) *[cos(2*pi*x)*sin(2*pi*y), cos(2*pi*y)*sin(2*pi*x),0] % have to calculate phi by hand
%M = cos(t) *[cos(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*z),sin(2*pi*x)*sin(2*pi*y)*cos(2*pi*z)] % have to calculate phi by hand
%M = t*t*(x*x - x)*[1,1,0]
%M = t*(x*x*x/3 - x*x/2)*[1,1,0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ddphi = 4*pi*(diff(M(1),x) + diff(M(2),y) + diff(M(3),z));
dphi = int(ddphi,x);
phi = int(dphi,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INPUT PHI (if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = x %% used to set phi when we want non-zero bc
%phi = 2*pi*cos(t)*x*x
%phi = 2*cos(t)*sin(2*pi*x)*sin(2*pi*y) %% phi for M = cos(t) *[cos(2*pi*x)*sin(2*pi*y), cos(2*pi*y)*sin(2*pi*x),0]
%phi = 2*cos(t)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
%% phi for M = cos(t) *[cos(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*z),sin(2*pi*x)*sin(2*pi*y)*cos(2*pi*z)]
%phi = 2*t*t*pi*((2/3)*x*x*x - x*x) % for M = t*t*(x*x - x)*[1,1,0]
%phi = (1/3)*pi*t*x*x*x*(x - 2) %% for M = t*(x*x*x/3 - x*x/2)*[1,1,0]
%phi = 2*cos(t)*sin(2*pi*x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dMdt = diff(M,t)

H_demag = [-diff(phi,x), -diff(phi,y), -diff(phi,z)]

dMx = [diff(M(1),x), diff(M(2),x), diff(M(3),x)];
dMy = [diff(M(1),y), diff(M(2),y), diff(M(3),y)];
dMz = [diff(M(1),z), diff(M(2),z), diff(M(3),z)];
H_ex = exchange_coeff*[diff(dMx(1),x) + diff(dMy(1),y) + diff(dMz(1),z), diff(dMx(2),x) + diff(dMy(2),y) + diff(dMz(2),z), diff(dMx(3),x) + diff(dMy(3),y) + diff(dMz(3),z)]

H = H_demag + H_ex
MxH = llg_precession_coeff*cross(M,H)
MxMxH = llg_damping_coeff*cross(M,cross(M,H))

  source = dMdt + MxH + MxMxH;

% check residual is zero:
simple(dMdt + MxH + MxMxH - source) == 0

ccode(M)

  ccode(simple(source))
