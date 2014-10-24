%
% pulse acquire experiment, using product operator formalism
%

operators; % setup operators

% initial state is I1z+I2z
rho = I1z + I2z;

% chemical shift of spins, omega=offset in rad/s
omega1 = 2*pi*1;
omega2 = 2*pi*-4;

% J-coupling between 1 and 2
J12 = .25; % in Hz, per convention

% Hamiltonians for the experiment
Hcs1 = omega1 * I1z;
Hcs2 = omega2 * I2z;
Hfree = Hcs1 + Hcs2 + 2*pi*J12*I1zI2z;

% pulse 90deg around the y-axis, should leave M in +x-direction
rho = arrow(rho, pi/2*I1y);

% do a "readout" - intermittant sampling of magnetization vectors
dt = 1/np; % sampling dwell-time
np = 128;  % number of points to record
% record M for ~1s
for ti=1:np
  m1(:,ti) = M1(rho);
  m2(:,ti) = M2(rho);
  mcoh(:,ti) = Mcoh2(rho);
  % freely evolve - apply Hfree for dt
  rho = arrow(rho, Hfree*dt);
end

% should produce a sine, starting at M=(1,0,0)
t = [0:(np-1)]/np;
clf
subplot(3,1,1)
plot(t, m1')
legend('Mx','My','Mz')
title('spin 1')
xlabel('s')
subplot(3,1,2)
plot(t, m2')
legend('Mx','My','Mz')
title('spin 2')
subplot(3,1,3)
plot(t, mcoh')
legend('I1xI2z', 'I1yI2z', 'I1zI2x', 'I1zI2y', ...
       'I1xI2x', 'I1xI2y', 'I1yI2x', 'I1yI2y', 'I1zI2z');
title('other coherences')
