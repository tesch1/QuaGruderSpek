% -*- matlab -*-
%
% Example of how to calculate quantum state propagation
% using Matlab.  Following the book "Understanding NMR
% Spectroscopy, second edition" by James Keeler
%
% Pulse-acq experiment (section 7.4.1, p. 151)
%
% 2016, Michael Tesch - michael.tesch@tum.de
%

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup single spin operators
%
% section 6.5.3 p. 111
Ix = 0.5*[0 1;
          1 0];
Iy = 0.5*[0 -1i;
          1i 0];
Iz = 0.5*[1 0;
          0 -1];
Ie  = eye(2);

% setup two spin operators
%
% section 7.4, p. 150

% z-magnetization
I1z = kron(Iz, Ie);
I2z = kron(Ie, Iz);

% in-phase
I1x = kron(Ix, Ie);
I1y = kron(Iy, Ie);
I2x = kron(Ie, Ix);
I2y = kron(Ie, Iy);

% anti-phase
I1xI2z = 2*kron(Ix, Iz);
I1yI2z = 2*kron(Iy, Iz);
I1zI2x = 2*kron(Iz, Ix);
I1zI2y = 2*kron(Iz, Iy);

% multiple-quantum coherence
I1xI2x = 2*kron(Ix, Ix);
I1xI2y = 2*kron(Ix, Iy);
I1yI2x = 2*kron(Iy, Ix);
I1yI2y = 2*kron(Iy, Iy);

% non-equilibrium population
I1zI2z = 2*kron(Iz, Iz);

% time propagation functions
%
% section 6.8.4, p. 129
arrow  = @(rho,H,t) expm(-1i * H * t) * rho * expm(1i * H * t);
arrowU = @(rho, U)  (U * rho * U');
makeU  = @(H,t)     expm(-1i * H * t);

% sample function - quadrature detection in Bruker is  < I- >
%
% section 13.6, p. 491
%
% take the real part to eliminate the tiny fp errors that end up in the imaginary part
%
%sample = @(rho)     (real(trace(rho*I1x)+trace(rho*I2x)) - i*real(trace(rho*I1y)+trace(rho*I2y)));
sample = @(rho)     (real(trace(rho*(I1x+I2x))) + i*real(trace(rho*(I1y+I2y))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup experiment
%

% fully relaxed system state is I1z + I2z
rho0 = I1z + I2z;

% chemical shift of spins, omega=offset in rad/s
omega1 = 2*pi * 40;
omega2 = 2*pi * 10;

% exponential decay rate
R = 2; % Hz
       %R = 0.5;

% J-coupling between 1 and 2
J12 = 3; % in Hz, per convention

% Hamiltonians for the coupled 2-spin experiment
Hcs1 = omega1 * I1z;
Hcs2 = omega2 * I2z;

% free evolution
Hfree = Hcs1 + Hcs2 + 2*pi*J12*I1zI2z;

% do a "readout" - sample the magnetization
%
% section 13.5.2, p. 490
np = 256;  % number of time points to record
dt = 1/np; % time between samples, total sample time = 1 s

% time propagators

% free evolution propagator
Udt = expm(-1i * Hfree * dt);

% ideal 90 degree pulse with y phase - "around y"
U90y = expm(-1i * pi/2 * (I1y+I2y) );

% pulse 90
rho = arrowU(rho0, U90y);

% readout
for ti=1:np
  % time elapsed since excitement - to roughly estimate relaxation
  tcur = (ti-1) * dt;
  % take measurement
  S(ti,1) = sample(rho) * exp(-tcur * R);
  % freely evolve to next sample time - apply Hfree for dt
  rho = arrowU(rho, Udt);
end

% processing

%Zfill = S;
%Zfill(:) = 0;
%S = [S ; Zfill];
spec = fft(S);
spec = fftshift(spec);
xticks = -size(spec,1)/2:size(spec,1)/2-1;
%xticks = xticks / 2;

%
% plot the spectrum
%
% p. 156, Fig 7.9

figure(4)
clf

subplot(3,1,1)
plot([real(S) imag(S)])
title('S')
legend('real(S)', 'imag(S)')

subplot(3,1,2)
hold on
plot(xticks, real(spec)./max(real(spec)))
spec2 = dct(real([S(2:end) ]));
np = length(spec2);
plot(linspace(0,np/2,np), spec2./max(abs(spec2)))
title('real(fft(S))')
xlabel('Hz')

subplot(3,1,3)
plot(xticks, imag(spec))
title('imag(fft(S))')
xlabel('Hz')
