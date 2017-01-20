% -*- matlab -*-
%
% Example of how to calculate quantum state propagation
% using Matlab.  Following the book "Understanding NMR
% Spectroscopy, second edition" by James Keeler
%
% COSY experiment (section 8.3, p. 190)
%
% 2016, Michael Tesch - michael.tesch@tum.de
%

% clear all variables from Matlab workspace
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time propagation functions
%
% section 6.8.4, p. 129
arrow  = @(rho,H,t) expm(-1i * H * t) * rho * expm(1i * H * t);
makeU  = @(H,t)     expm(-1i * H * t);
arrowU = @(rho, U)  (U * rho * U');

% sample function - quadrature detection in Bruker is of (I-)
%
% section 13.6, p. 491
%sample = @(rho)     (real(trace(rho*I1x)+trace(rho*I2x)) - i*real(trace(rho*I1y)+trace(rho*I2y)));
sample = @(rho)     (real(trace(rho*(I1x+I2x))) + i*real(trace(rho*(I1y+I2y))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup COSY experiment
%

% fully relaxed system state is I1z + I2z
rho0 = I1z + I2z;
%rho0 = I1z;

% chemical shift of spins, omega=offset in rad/s
omega1 = 2*pi * 40;
omega2 = 2*pi * 10;

% exponential decay rate
R = 4; % Hz

% J-coupling between 1 and 2
J12 = 5; % in Hz, per convention

% Hamiltonians for the coupled 2-spin experiment
Hcs1 = omega1 * I1z;
Hcs2 = omega2 * I2z;
% free evolution
Hfree = Hcs1 + Hcs2 + 2*pi*J12*I1zI2z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a "readout" - sample the magnetization over time
%
% section 13.5.2, p. 490
np = 128;  % number of time points to record
dt = 1.0/np; % time between samples, total sample time = 1.0 s

% time propagators

% free evolution between samples
Udt = expm(-1i * Hfree * dt);
% 90 degree pulse with x phase - "around x"
U90x = expm(-1i * pi/2 * (I1x+I2x) );

% 2-dimensional experiment, t1 and t2
for ti1=1:np
  % rho relaxed back to resting state
  rho = rho0;
  % pulse 90deg around the x-axis
  rho = arrowU(rho, U90x);
  % evolve for t1
  t1 = dt * (ti1-1);
  Ut1 = expm(-1i * Hfree * t1);
  rho = arrowU(rho, Ut1);
  % second pulse
  rho = arrowU(rho, U90x);

  % readout
  for ti2=1:np
    % time elapsed since excitement - to roughly estimate relaxation
    tcur = t1 + (ti2-1)*dt;
    % take measurement
    S(ti1,ti2) = sample(rho) * exp(-tcur * R);
    % freely evolve to next sample time - apply Hfree for dt
    rho = arrowU(rho, Udt);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing

figure(2)
clf

subplot(2,2,1)
% p. 192, "The cross-peak multiplet"
% second dimension FT
Stw = fft(S,[],2);
Stw = fftshift(Stw,2);
Stw = -1i * Stw;  % exp(-i * (pi/2)) : p. 193, eq 8.2 -> eq 8.3
% in matlab, the inverse dct has the right alignment
Sww = idct(real(Stw));

spec = Sww;
xticks = -size(spec,2)/2:size(spec,2)/2-1;
yticks = [0:size(spec,1)-1]/2;
proctext = 'COSY, cross peaks in "anti-phase square array"';
contour(xticks, yticks, spec, 30)
%imagesc(xticks, yticks, spec)
xlabel('\omega_2 (Hz)')
ylabel('\omega_1 (Hz)')
title(proctext)
colorbar

subplot(2,2,3)
surf(xticks,yticks,spec)
zoom(1.5);

subplot(2,2,2)
% p. 195, "The diagonal-peak multiplet"
% first dimension FT
spec = fft(S,[],2);
spec = fftshift(spec,2);
spec = dst(real(spec(2:end,:)));
yticks = [0:size(spec,1)-1]/2;
proctext = 'COSY, diagonal peaks in double absorption';
contour(xticks, yticks, real(spec), 30)
%imagesc(xticks, yticks, real(spec))
xlabel('\omega_2 (Hz)')
ylabel('\omega_1 (Hz)')
title(proctext)
colorbar

subplot(2,2,4)
surf(xticks,yticks,spec)
zoom(1.5);

