%
% Quantenmechanische Grundlagen der NMR-Spektroskopie
%                  "QuaGruderSpek"
%
% 

% running this file creates operators for 1,2,3-spin systems
% in the MATLAB environment

% scaling factors
Iscale = 0.5;
IIscale = 2;
IIIscale = 4;

% cartesian operators for single 1/2-spin (2x2)
% E1 = unity operator, single spin
E1  = eye(2);
Ix = Iscale*[0 1; 1 0];
Iy = Iscale*[0 -1i; 1i 0];
Iz = Iscale*[1 0; 0 -1];

% operators for two coupled 1/2-spin systems (4x4)

% q=0

% E2 = unity operator, two spin
E2 = Iscale*kron(E1,E1);

% q=1

I1x = kron(Ix,E1);
I1y = kron(Iy,E1);
I1z = kron(Iz,E1);
I2x = kron(E1,Ix);
I2y = kron(E1,Iy);
I2z = kron(E1,Iz);

% q=2

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

% q=0

% E3 = unity operator, three spin
E3 = Iscale*kron(E1,kron(E1,E1));

% q=1

S1x = kron(E1,kron(Ix,E1));
S1y = kron(E1,kron(Iy,E1));
S1z = kron(E1,kron(Iz,E1));
S2x = kron(E1,kron(E1,Ix));
S2y = kron(E1,kron(E1,Ix));
S2z = kron(E1,kron(E1,Iz));
S3x = kron(kron(Ix,E1),E1);
S3y = kron(kron(Iy,E1),E1);
S3z = kron(kron(Iz,E1),E1);

% q=2

S1xS2xE = IIscale*kron(Ix,kron(Ix,E1));
S1xS2yE = IIscale*kron(Ix,kron(Iy,E1));
S1xS2zE = IIscale*kron(Ix,kron(Iz,E1));
S1xES3x = IIscale*kron(Ix,kron(E1,Ix));
S1xES3y = IIscale*kron(Ix,kron(E1,Iy));
S1xES3z = IIscale*kron(Ix,kron(E1,Iz));

S1yES3x = IIscale*kron(Iy,kron(E1,Ix));
S1yES3y = IIscale*kron(Iy,kron(E1,Iy));
S1yES3z = IIscale*kron(Iy,kron(E1,Iz));
S1yS2xE = IIscale*kron(Iy,kron(Ix,E1));
S1yS2yE = IIscale*kron(Iy,kron(Iz,E1));
S1yS2zE = IIscale*kron(Iy,kron(Iz,E1));

S1zES3x = IIscale*kron(Iz,kron(E1,Ix));
S1zES3y = IIscale*kron(Iz,kron(E1,Iy));
S1zES3z = IIscale*kron(Iz,kron(E1,Iz));
S1zS2xE = IIscale*kron(Iz,kron(Ix,E1));
S1zS2yE = IIscale*kron(Iz,kron(Iy,E1));
S1zS2zE = IIscale*kron(Iz,kron(Iz,E1));

ES2xS3x = IIscale*kron(E1,kron(Ix,Ix));
ES2xS3y = IIscale*kron(E1,kron(Ix,Iy));
ES2xS3z = IIscale*kron(E1,kron(Ix,Iz));
ES2yS3x = IIscale*kron(E1,kron(Iy,Ix));
ES2yS3z = IIscale*kron(E1,kron(Iy,Iy));
ES2yS3y = IIscale*kron(E1,kron(Iy,Iz));
ES2zS3x = IIscale*kron(E1,kron(Iz,Ix));
ES2zS3y = IIscale*kron(E1,kron(Iz,Iy));
ES2zS3z = IIscale*kron(E1,kron(Iz,Iz));

% q=3

S1xS2xS3x = IIIscale*kron(Ix,kron(Ix,Ix));
S1xS2yS3x = IIIscale*kron(Ix,kron(Iy,Ix));
S1xS2zS3x = IIIscale*kron(Ix,kron(Iz,Ix));
S1xS2xS3y = IIIscale*kron(Ix,kron(Ix,Iy));
S1xS2yS3y = IIIscale*kron(Ix,kron(Iy,Iy));
S1xS2zS3y = IIIscale*kron(Ix,kron(Iz,Iy));
S1xS2xS3z = IIIscale*kron(Ix,kron(Ix,Iy));
S1xS2yS3z = IIIscale*kron(Ix,kron(Iy,Iz));
S1xS2zS3z = IIIscale*kron(Ix,kron(Iz,Iz));

S1yS2xS3x = IIIscale*kron(Iy,kron(Ix,Ix));
S1yS2yS3x = IIIscale*kron(Iy,kron(Iy,Ix));
S1yS2zS3x = IIIscale*kron(Iy,kron(Iz,Ix));
S1yS2xS3y = IIIscale*kron(Iy,kron(Ix,Iy));
S1yS2yS3y = IIIscale*kron(Iy,kron(Iy,Iy));
S1yS2zS3y = IIIscale*kron(Iy,kron(Iz,Iy));
S1yS2xS3z = IIIscale*kron(Iy,kron(Ix,Iz));
S1yS2yS3z = IIIscale*kron(Iy,kron(Iy,Iz));
S1yS2zS3z = IIIscale*kron(Iy,kron(Iz,Iz));

S1zS2xS3x = IIIscale*kron(Iz,kron(Ix,Ix));
S1zS2yS3x = IIIscale*kron(Iz,kron(Iy,Ix));
S1zS2zS3x = IIIscale*kron(Iz,kron(Iz,Ix));
S1zS2xS3y = IIIscale*kron(Iz,kron(Ix,Iy));
S1zS2yS3y = IIIscale*kron(Iz,kron(Iy,Iy));
S1zS2zS3y = IIIscale*kron(Iz,kron(Iz,Iy));
S1zS2xS3z = IIIscale*kron(Iz,kron(Ix,Iz));
S1zS2yS3z = IIIscale*kron(Iz,kron(Iy,Iz));
S1zS2zS3z = IIIscale*kron(Iz,kron(Iz,Iz));

% the arrow operator
% ie: I1y --- 2*pi*I1x*t ---> I1y*cos(2*pi*t) + I1z * sin(2*pi*t)
%
% example: subject I1x to the Hamiltonian 2*pi*I1z for t=0.25):
% arrow(I1x, 2*pi*I1z*0.25)
%
arrow = @(op, H) (expm(-1i * H) * op * expm(1i * H));
makeU = @(H) expm(-1i * H);
arrowU = @(op, U) (U * op * U');

% measure the expected 3d magnetization, given a density operator rho
% (due to floating-point error, imag part can be very small, so
%  we use real() to force it to 0)
M = @(rho) real([trace(rho*Ix) trace(rho*Iy) trace(rho*Iz)]);

% M for spins 1 and 2 in a coupled 2-spin-1/2 system (single-quantum coherences)
M1 = @(rho) real([trace(rho*I1x) trace(rho*I1y) trace(rho*I1z)]);
M2 = @(rho) real([trace(rho*I2x) trace(rho*I2y) trace(rho*I2z)]);
N1x = @(rho) real([trace(rho*I1x)]);
N1y = @(rho) real([trace(rho*I1y)]);
N1z = @(rho) real([trace(rho*I1z)]);
N2x = @(rho) real([trace(rho*I2x)]);
N2y = @(rho) real([trace(rho*I2y)]); 
N2z = @(rho) real([trace(rho*I2z)]);

% other magnetization coherences for 2-spin-1/2 system
Mcoh2 = @(rho) real([...
    trace(rho*I1xI2z) ... % anti-phase magnetization spin 1
    trace(rho*I1yI2z) ...
    trace(rho*I1zI2x) ... % anti-phase magnetization spin 2
    trace(rho*I1zI2y) ...
    trace(rho*I1xI2x) ... % multiple quantum coherences
    trace(rho*I1xI2y) ...
    trace(rho*I1yI2x) ...
    trace(rho*I1yI2y) ...
    trace(rho*I1zI2z) ]); % non-equlibrium population

% constants to index into Mcoh2() generated measurements
idxI1xI2z = 1;
idxI1yI2z = 2;
idxI1zI2x = 3;
idxI1zI2y = 4;
idxI1xI2x = 5;
idxI1xI2y = 6;
idxI1yI2x = 7;
idxI1yI2y = 8;
idxI1zI2z = 9;

% M for spins 1,2,3 in a coupled 3-spin-1/2

MS1 = @(rho) real([trace(rho*S1x) trace(rho*S1y) trace(rho*S1z)])*IIscale;
MS2 = @(rho) real([trace(rho*S2x) trace(rho*S2y) trace(rho*S2z)])*IIscale;
MS3 = @(rho) real([trace(rho*S3x) trace(rho*S3y) trace(rho*S3z)])*IIscale;

% other magnetization coherences for 3-spin-1/2-system

Mcoh3 = @(rho) ([...
    trace(rho*S1xS2xE)...               % 2-spin-coherence
    trace(rho*S1xES3x)...
    trace(rho*S1yS2yE)...
    trace(rho*S1yES3y)...
    trace(rho*S1zS2zE)...
    trace(rho*S1zES2z)...
    trace(rho*ES2xS3x)...
    trace(rho*ES2yS3y)...
    trace(rho*ES2zS3z)...
    trace(rho*S1xS2yS3x)...
    trace(rho*S1xS2zS3x)...
    trace(rho*S1xS2xS3y)...
    trace(rho*S1xS2yS3y)...
    trace(rho*S1xS2xS3z)...
    trace(rho*S1xS2zS3z)...
    trace(rho*S1yS2xS3x)...
    trace(rho*S1yS2yS3x)...
    trace(rho*S1yS2xS3y)...
    trace(rho*S1yS2yS3z)...
    trace(rho*S1yS2zS3z)
    ]);
                                          % 3-spin-coherence
