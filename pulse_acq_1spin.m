%
% pulse acquire experiment, using product operator formalism
%

operators; % setup operators

% initial state is Iz
rho = Iz;

% pulse 90deg around the y-axis, should leave M in +x-direction
rho = arrow(rho, pi/2*Iy);

% chemical shift Hamiltonian, omega=offset in rad/s
omega = 2*pi*1;
Hcs = omega * Iz;

% do a "readout" - intermittant sampling of magnetization vectors
dt = 1/np; % sampling dwell-time
np = 128;  % number of points to record
% record M for ~1s
for ti=1:np
  m1(:,ti) = M(rho);
  % freely evolve - apply Hcs for dt
  rho = arrow(rho, Hcs*dt);
end

% should produce a sine, starting at M=(1,0,0)
t = [0:(np-1)]/np;
clf
plot(t, m1')
legend('Mx','My','Mz')
xlabel('s')
title('magnetization after 90y pulse')
