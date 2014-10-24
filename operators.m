%
% Quantenmechanische Grundlagen der NMR-Spektroskopie
% QuaGruderSpek
%

% setup operators

%useful constants
%RT2 = sqrt(2.0);
RT2 = sym('RT2');
RT2r = 1./RT2; % r for reciprocal

% symbolics
a=sym('a');
b=sym('b');
psi_ab = [a; b];

Iscale = 0.5;
IIscale = 2;
IIIscale = 2;
%Iscale = 1;
%IIscale = 1;

% cartesian operators for single 1/2-spin (2x2)
II  = eye(2);
Ix = Iscale*[0 1; 1 0];
Iy = Iscale*[0 -1i; 1i 0];
Iz = Iscale*[1 0; 0 -1];

% operators for two coupled 1/2-spin systems (4x4)
I1x = kron(Ix,II);
I1y = kron(Iy,II);
I1z = kron(Iz,II);
I2x = kron(II,Ix);
I2y = kron(II,Iy);
I2z = kron(II,Iz);
I1xI2x = IIscale*kron(Ix,Ix);
I1xI2y = IIscale*kron(Ix,Iy);
I1xI2z = IIscale*kron(Ix,Iz);
I1yI2x = IIscale*kron(Iy,Ix);
I1yI2y = IIscale*kron(Iy,Iy);
I1yI2z = IIscale*kron(Iy,Iz);
I1zI2x = IIscale*kron(Iz,Ix);
I1zI2y = IIscale*kron(Iz,Iy);
I1zI2z = IIscale*kron(Iz,Iz);

% operators for coupled 3 spin-1/2 sytem (8x8)


% the arrow operator
% ie: I1y --- 2*pi*I1x*t ---> I1y*cos(2*pi*t) + I1z * sin(2*pi*t)
%
% example: subject I1x to the Hamiltonian 2*pi*I1z for t=0.25):
% arrow(I1x, 2*pi*I1z*0.25)
%
arrow = @(op, H) (expm(-1i * H) * op * expm(1i * H));

% measure the expected 3d magnetization, given a density operator rho
% (due to floating-point error, imag part can be very small: force it to 0)
M = @(rho) real([trace(rho*Ix) ; trace(rho*Iy) ; trace(rho*Iz)]);

% M for spins 1 and 2 in a coupled 2-spin-1/2 system (single-quantum coherences)
M1 = @(rho) real([trace(rho*I1x) ; trace(rho*I1y) ; trace(rho*I1z)]);
M2 = @(rho) real([trace(rho*I2x) ; trace(rho*I2y) ; trace(rho*I2z)]);

% other magnetization coherences for 2-spin-1/2 system
Mcoh2 = @(rho) real([...
    trace(rho*I1xI2z) ; ... % anti-phase magnetization spin 1
    trace(rho*I1yI2z) ; ...
    trace(rho*I1zI2x) ; ... % anti-phase magnetization spin 2
    trace(rho*I1zI2y) ; ...
    trace(rho*I1xI2x) ; ... % multiple quantum coherences
    trace(rho*I1xI2y) ; ...
    trace(rho*I1yI2x) ; ...
    trace(rho*I1yI2y) ; ...
    trace(rho*I1zI2z) ]); % non-equlibrium population

% M for spins 1,2,3 in a coupled 3-spin-1/2
