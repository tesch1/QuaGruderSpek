% Ubungsaufgaben zur Vorlesung

% load the operators into matlab
operators

% 1)
psi = (1./sqrt(2.)) * [1 ; i];
% expectation of operator I{x,y,z}:
expectation_Ix = psi'*Ix*psi
expectation_Iy = psi'*Iy*psi
expectation_Iz = psi'*Iz*psi

% or...
M(psi*psi')

% 2)
