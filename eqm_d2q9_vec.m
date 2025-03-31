function f_eq = eqm_d2q9_vec(Rho,U)
%% Function f_eq = eqm_d2q9_vec(Rho,U) computes the equalibrium PDF based on the given Density and Velocity at any location using D2Q9 lattice
%% Rho is the density and must be a 1 x M x N matrix
%% U is the velocity and must be a 2 x 9 x M x N matrix
%% f_eq is the equalibrium PDF and will be a 9 x M x N matrix

% Define latteral velocity xi (E), weight (w), speed of sound (c_s)
E = [0  0;
     1  0;
     0  1;
    -1  0;
     0 -1;
     1  1;
    -1  1;
    -1 -1;
     1 -1];
w = [4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36];
w = repmat(w, [1 size(Rho, [2 3])]);
% c_s = 1/sqrt(3);


f_eq = w.*Rho.*(1 + 3*(pagemtimes(E,U)) + 9*(pagemtimes(E,U)).^2/2 - 3*(dot(U,U)/2));