%
% pulse acquire experiment, using product operator formalism
%

operators; % setup operators

clear m1 m2 s1 s2 mcoh

% initial state is I1z+I2z
rho = I1z + I2z;

% chemical shift of spins, omega=offset in rad/s
omega1 = 2*pi*20;
omega2 = 2*pi*-10;

% J-coupling between 1 and 2
J12 = .25; % in Hz, per convention
J12 = 4; % in Hz, per convention

% Hamiltonians for the coupled 2-spin experiment
Hcs1 = omega1 * I1z;
Hcs2 = omega2 * I2z;
Hfree = Hcs1 + Hcs2 + 2*pi*J12*I1zI2z;

% pulse 90deg around the y-axis, should leave M in +x-direction
rho = arrow(rho, pi/2*I1y);

% do a "readout" - intermittant sampling of magnetization vectors
np = 128;  % number of points to record
dt = 1/np; % sampling dwell-time
% record M for ~1s
for ti=1:np
  m1(ti,:) = M1(rho);
  m2(ti,:) = M2(rho);
  mcoh(:,ti) = Mcoh2(rho);
  % freely evolve - apply Hfree for dt
  rho = arrow(rho, Hfree*dt);
end

% should start at M1=(1,0,0)
t = [0:(np-1)]/np;
clf
subplot(5,1,1)
plot(t, m1')
legend('Mx','My','Mz')
title('spin 1')
xlabel('t (s)')

subplot(5,1,2)
plot(t, m2')
legend('Mx','My','Mz')
title('spin 2')
xlabel('t (s)')

% what's going on with the coherences?
subplot(5,1,3)
plot(t, real(mcoh)')
legend('I1xI2z', 'I1yI2z', 'I1zI2x', 'I1zI2y', ...
       'I1xI2x', 'I1xI2y', 'I1yI2x', 'I1yI2y', 'I1zI2z');
title('other coherences (real)')
xlabel('t (s)')

subplot(5,1,4)
plot(t, imag(mcoh)')
legend('I1xI2z', 'I1yI2z', 'I1zI2x', 'I1zI2y', ...
       'I1xI2x', 'I1xI2y', 'I1yI2x', 'I1yI2y', 'I1zI2z');
title('other coherences (imag)')
xlabel('t (s)')

% plot the frequency spectrum
subplot(5,1,5)
s1 = fftshift(fft(m1,[],1));
bw = 1/dt;
f = linspace(-bw/2,bw/2,np);
hold on
plot(f, real(s1(:,1)),'bx-')
plot(f, real(s1(:,2)),'g+-')
plot(f, real(s1(:,3)),'r*-')
plot(f, imag(s1(:,1)),'cx-')
plot(f, imag(s1(:,2)),'k+-')
plot(f, imag(s1(:,3)),'m*-')
legend('Mx_{re}','My_{re}','Mz_{re}','Mx_{im}','My_{im}','Mz_{im}')
title(['spin 1: spectrum after 90y pulse, omegaHz=' num2str(omega1/(2*pi)) ' J12=' num2str(J12)])
xlabel('f (Hz)')
