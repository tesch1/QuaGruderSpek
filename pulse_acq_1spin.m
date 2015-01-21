%
% pulse acquire experiment, using product operator formalism
%

operators; % setup product operators

% initial state is Iz
rho = Iz;

% pulse 90deg around the y-axis, should leave M in +x-direction
rho = arrow(rho, pi/2*Iy);

% chemical shift Hamiltonian, omega=offset in rad/s
omega = 2*pi*23;
Hcs = omega * Iz;

% do a "readout" - intermittant sampling of magnetization vectors
np = 128;  % number of points to record
dt = 1/np; % sampling dwell-time

% record M for ~1s
for ti=1:np
  m1(ti,:) = M(rho);
  % freely evolve - apply Hcs for dt
  rho = arrow(rho, Hcs*dt);
end

%
% plot the results
%

% should produce a sine, starting at M=(1,0,0)
clf
subplot(2,1,1)
t = [0:(np-1)]/np;
plot(t, m1)
legend('Mx','My','Mz')
title(['single spin: magnetization after 90y pulse, omegaHz=' num2str(omega/(2*pi))])
xlabel('t (s)')

% plot the frequency spectrum
subplot(2,1,2)
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
legend('Mx_re','My_re','Mz_re','Mx_im','My_im','Mz_im')
title(['single spin: spectrum after 90y pulse, omegaHz=' num2str(omega/(2*pi))])
xlabel('f (Hz)')
