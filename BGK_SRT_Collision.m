clc; clear; close all

% define f vectors
f1 = [1.63, 0.61, 0.41, 0.27, 0.41, 0.15, 0.07, 0.07, 0.16];
f2 = [1.67, 0.42, 0.42, 0.42, 0.42, 0.1, 0.11, 0.1, 0.11];
f3 = [1.66, 0.5, 0.42, 0.35, 0.42, 0.12, 0.09, 0.08, 0.13];

% Define relaxation time constant
t = 0.6;    % must be > 0.5

% Compute
[Rho1, U1] = moment_rho_U_d2q9(f1);
[Rho2, U2] = moment_rho_U_d2q9(f2);
[Rho3, U3] = moment_rho_U_d2q9(f3);

f_eq1 = eqm_d2q9(Rho1, U1);
f_eq2 = eqm_d2q9(Rho2, U2);
f_eq3 = eqm_d2q9(Rho3, U3);

fstar1 = fstar(f1, f_eq1, t);
[Rhostar1, Ustar1] = moment_rho_U_d2q9(fstar1);

%%% PDF loop

% Define iterations
iter = 1e3;

% Initialize matricies
f = [f1; zeros(iter-1,9)];
Rho = zeros(iter,1);
U = zeros(2,iter);

% loop
for n = 1:iter-1
    [Rho(n), U(:,n)] = moment_rho_U_d2q9(f(n,:));
    f_eq = eqm_d2q9(Rho(n), U(:,n));
    f(n+1,:) = fstar(f(n,:), f_eq, t);
end
[Rho(iter), U(:,iter)] = moment_rho_U_d2q9(f(iter,:));

% Plot
iter_vec = 1:iter;

% f components
for n = 1:length(f1)
    figure
    plot(iter_vec,f(:,n))
    title(sprintf('f%d components over %d iterations', n, iter));
    xlabel('Iterations')
    ylabel(sprintf('f%d component values',n))
end

% Density
figure
plot(iter_vec,Rho);
title(sprintf('Density over %d iterations',iter));
xlabel('Iterations')
ylabel('Density')

% Velocity
xy = ['x', 'y'];
for n = 1:length(U1)
    figure
    plot(iter_vec,U(n,:))
    title(sprintf('U%s components over %d iterations', xy(n), iter));
    xlabel('Iterations')
    ylabel(sprintf('U%s values',xy(n)))
end