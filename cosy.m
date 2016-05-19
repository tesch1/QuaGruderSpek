%
% pulse acquire experiment, using product operator formalism
%

operators; % setup operators

clear m1 m2 s1 s2 mcoh

% initial fully relaxed state is I1z+I2z
rho0 = I1z + I2z;

% chemical shift of spins, omega=offset in rad/s
omega1 = 2*pi*40;
omega2 = 2*pi*10;

% relaxation time constant
T2 = 0.5;

% J-coupling between 1 and 2
J12 = 3; % in Hz, per convention

% Hamiltonians for the coupled 2-spin experiment
Hcs1 = omega1 * I1z;
Hcs2 = omega2 * I2z;
Hfree = Hcs1 + Hcs2 + 2*pi*J12*I1zI2z;

% do a "readout" - intermittant sampling of magnetization vectors
np = 128;  % number of points to record
dt = 1/np; % sampling dwell-time

% record M for ~1s
U = makeU(Hfree, dt);
U90x = makeU(pi/2*(I1x+I2x), 1);

for ti1=1:np
  % rho relaxes back to resting state
  rho = rho0;
  % pulse 90deg around the x-axis
  rho = arrowU(rho, U90x);
  % evolve for t1
  t1 = dt * ti1;
  rho = arrow(rho, Hfree, t1);
  % second pulse
  rho = arrowU(rho, U90x);

  % readout
  for ti2=1:np
    % take measurement
    tcur = t1 + ti2*dt;
    m1(ti1,ti2) = meas2(rho) * exp((-tcur) / T2);
    % freely evolve - apply Hfree for dt
    %rho = arrow(rho, Hfree*dt); <- old way, calculate expm(-i*Hfree*dt) every time
    rho = arrowU(rho, U);
  end
end

% plot the frequency spectrum
% TODO: get rid of signal alias from COSY's mixed sin/cos modulation
clf
subplot(1,2,1)
contour(-np/2:np/2-1, -64:63,real(fftshift(fft2(m1))),30)
line([-np/2 np/2-1], [-np/2 np/2-1])
xlabel('Hz')
ylabel('Hz')
title('COSY Absorption Mode')

subplot(1,2,2)
contour(-np/2:np/2-1, -64:63,imag(fftshift(fft2(m1))),30)
line([-np/2 np/2-1], [-np/2 np/2-1])
xlabel('Hz')
ylabel('Hz')
%contour(imag(fftshift(fft2(m1))),20)
title('COSY Dispersion Mode')
colorbar

%xlabel('f1 (Hz)')
%ylabel('f2 (Hz)')

