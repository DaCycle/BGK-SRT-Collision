clc; clear; close all;
% test

% Set M and N
M = 9;
N = 7;

% Define relaxation time constant
t = 0.6;    % must be > 0.5

% Make f
f = [1.63 1.67 1.66;
     0.61 0.42 0.50; 
     0.41 0.42 0.42;
     0.27 0.42 0.35;
     0.41 0.42 0.42;
     0.15 0.10 0.12;
     0.07 0.11 0.09;
     0.07 0.10 0.08;
     0.16 0.11 0.13];
f = repmat(f, [1 M/3 N]);

% Find denisity and velocity
[Rho, U] = moment_rho_U_d2q9_vec(f);

% Find f_eq
f_eq = eqm_d2q9_vec(Rho, U);

% Find f_star
f_star = fstar(f, f_eq, t);
[Rhostar, Ustar] = moment_rho_U_d2q9_vec(f_star);

%%% PDF loop

% Define iterations
iter = 1e3;

% Initialize matricies
f_loop = zeros(9, M, N, iter);
f_loop(:, :, :, 1) = f;
Rho_loop = zeros(1, M, N, iter);
U_loop = zeros(2, M, N, iter);

% loop
for n = 1:iter-1
    [Rho_loop(:,:,:,n), U_loop(:,:,:,n)] = moment_rho_U_d2q9_vec(f_loop(:,:,:,n));
    f_eq = eqm_d2q9_vec(Rho_loop(:,:,:,n), U_loop(:,:,:,n));
    f_loop(:,:,:,n+1) = fstar(f_loop(:,:,:,n), f_eq, t);
end
[Rho_loop(:,:,:,iter), U_loop(:,:,:,iter)] = moment_rho_U_d2q9_vec(f_loop(:,:,:,iter));

% % Plot
% iter_vec = 1:iter;
% 
% % f components
% for n = 1:length(f1)
%     figure
%     plot(iter_vec,f(:,n))
%     title(sprintf('f%d components over %d iterations', n, iter));
%     xlabel('Iterations')
%     ylabel(sprintf('f%d component values',n))
% end
% 
% % Density
% figure
% plot(iter_vec,Rho);
% title(sprintf('Density over %d iterations',iter));
% xlabel('Iterations')
% ylabel('Density')
% 
% % Velocity
% xy = ['x', 'y'];
% for n = 1:length(U1)
%     figure
%     plot(iter_vec,U(n,:))
%     title(sprintf('U%s components over %d iterations', xy(n), iter));
%     xlabel('Iterations')
%     ylabel(sprintf('U%s values',xy(n)))
% end
