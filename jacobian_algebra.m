


% clear                                   %

% syms m1 m2 m3 hex1 hex2 hex3 alpha;

% m = [m1,m2,m3];
% hex = [hex1,hex2,hex3];

% % syms happ1 happ2 happ3 hca1 hca2 hca3
% % hms = [diff(phi,x),diff(phi,y),diff(phi,z)]
% % hca = [hca1,hca2,hca3];
% % happ = [happ1,happ2,happ3];
% % h = hms + happ + hex + hca;

% % The llg residual ignoring time dependence, assuming gyromagnetic
% % constant is normalised out and only field is exchange.
% llg_res = cross(m,hex) + alpha * cross(m,cross(m,hex));

% % Then the Jacobian is just:
% J_hex = [diff(llg_res,hex1)
%          diff(llg_res,hex2)
%          diff(llg_res,hex3)];

% % Similarly for the derivatives wrt m (just using total h for simplicity)
% syms h1 h2 h3;
% h = [h1,h2,h3];

% llg_res = cross(m,h) + alpha * cross(m,cross(m,h));

% J_m = [diff(llg_res,m1)
%        diff(llg_res,m2)
%        diff(llg_res,m3)];

% syms x y z
% %syms phi(x,y,z) % need matlab 2012a?
% syms dphidx dphidy dphidz
% phi = x*dphidx + y*dphidy + z*dphidz; % hack to do it with this version

% h_ms = [diff(phi,x),diff(phi,y),diff(phi,z)];

% llg_res = cross(m,h_ms) + alpha * cross(m,cross(m,h_ms));

% J_phi = diff(llg_res,phi)

% J = [J_phi
%      J_phi_1
%      J_m
%      J_hex];

% ccode(simple(J))

clear
syms m1 m2 m3 hex1 hex2 hex3 alpha phi phi1 happ hca h1 h2 h3;
m = [m1,m2,m3];
h = [h1,h2,h3];
hex = [hex1,hex2,hex3];

syms C F G w0; % C F and G are the corresponding entry in a variety of mass
             % matrix-style inner products. w0 is the weight of the present
             % time m in the time stepping scheme. See notebook 10/5/2012!

rllg = C*w0*[m1,m2,m3] + F*cross(m,h) + G*alpha*cross(m,cross(m,h));
% + other terms in earlier/later time steps which are irrelevant here

dllgdm = [diff(rllg(1),m1), diff(rllg(1),m2), diff(rllg(1),m3);
          diff(rllg(2),m1), diff(rllg(2),m2), diff(rllg(2),m3);
          diff(rllg(3),m1), diff(rllg(3),m2), diff(rllg(3),m3)]

% To get dependence on hex let h = hex then:
h = [hex1,hex2,hex3];
rllg = C*w0*[m1,m2,m3] + F*cross(m,h) + G*alpha*cross(m,cross(m,h));

dllgdhex = [diff(rllg(1),hex1), diff(rllg(1),hex2), diff(rllg(1),hex3);
            diff(rllg(2),hex1), diff(rllg(2),hex2), diff(rllg(2),hex3);
            diff(rllg(3),hex1), diff(rllg(3),hex2), diff(rllg(3),hex3)]

% Similarly for dependence on phi:
syms phi x y z;
h = - [diff(phi,x), diff(phi,y), diff(phi,z)]
rllg = C*w0*[m1,m2,m3] + F*cross(m,h) + G*alpha*cross(m,cross(m,h));

dllgdphi = [diff(rllg(1),phi);
            diff(rllg(2),phi);
            diff(rllg(3),phi)]
