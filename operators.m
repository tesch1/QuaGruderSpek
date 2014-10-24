%
% Quantenmechanische Grundlagen der NMR-Spektroskopie
% QuaGruderSpek
%

% operators

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
%Iscale = 1;
%IIscale = 1;

% cartesian operators for single 1/2-spin
II  = eye(2);
Ix = Iscale*[0 1; 1 0];
Iy = Iscale*[0 -1i; 1i 0];
Iz = Iscale*[1 0; 0 -1];

% operators for two coupled 1/2-spins
I1x = kron(Ix,II);
I1y = kron(Iy,II);
I1z = kron(Iz,II);
I2x = kron(II,Iy);
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

% the arrow operator, ie: I1y --- 2*pi*I1x*t ---> I1y*cos(2*pi*t) + I1z * sin(2*pi*t)
arrow = @(op, H) (expm(-1i * H) * op * expm(1i * H));

% example: subject I1x to the Hamiltonian 2*pi*I1z for t=0.25)
% arrow(I1x, 2*pi*I1z*0.25)

