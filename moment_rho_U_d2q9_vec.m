function [Rho, U] = moment_rho_U__d2q9_vec(f)

% Define lattice velocity Xi (E)
E = [0, 1, 0, -1, 0, 1, -1, -1, 1;...
     0, 0, 1, 0, -1, 1, 1, -1, -1];

% Find Rho
Rho = sum(f,1);

% Find U
U = pagemtimes(E,f)./Rho;