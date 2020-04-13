%% APM 522 Group Project : Electrical Impedence Tomography
% Roberto Alvarez
% Bechir Amdouni
% Andrew Mullen
% Trevor Reckell
% Shan Zhu
clear;
clc;
N = 5;
%% Forward Problem
% calculate true sigma
[X,Y]= meshgrid(1:N+1,1:N+1);
sigmat = ((cos(X)) + sin(Y)).^2; % true conductivity
% setup boundary conditions
% [b1,b2,b3,b4] = forward_simulation(N,sigma);
b1 = abs(normrnd(1,0.05, N-2,1));
b2 = abs(normrnd(1,0.05, 1,N-2));
b3 = abs(normrnd(1,0.05, 1,N-2));
b4 = abs(normrnd(1,0.05, N-2,1));
% solve forward problem and find true solution
tol = 1e-14;
maxiter = 1e6;
ut = forward(N,sigmat,b1,b2,b3,b4,tol,maxiter); % true solution
%% Inverse Problem
% Gauss Newton Method
% 1 - Guess initial value for conductivity (sigma0)
% 2 - Find voltage (u) using forward model
% 3 - Calculate residual: (r_i = u_i - u_true)
% 4 - Update sigma: sigma_(i+1) = sigma_i + (J'J)^-1 * J' * r_i
% 5 - Repeat 2-4 until we are near the true value for sigma
sigma0 = 2*ones(N+1); % initial guess
u = forward(N,sigma0,b1,b2,b3,b4,tol,maxiter); % solve with initial guess
r = abs(u - ut); % calculate residual
% J = zeros((N-1)*(N-1)); % allocate space for jacobian
% for i = 1:(N-1)*(N-1)
%     J(:,i)=u_sol(i)*sigma0*A(:,i); % calculate jacobian
% end
% sigma = sigma0 + (J'*J)^-1*(J'*r); % update sigma
% sigma_tol = 1e-20; % how far we want reconstruction to be from true value
% it = 0; % iteration number
% while r > sigma_tol
%     u = forward(N,sigma,b1,b2,b3,b4,tol,maxiter); % solve forward problem
%     r = abs(u-ut); % calculate residual
%     for i = 1:(N-1)*(N-1)
%         J(:,i) = u_sol(i)*sigma*A(:,i); % update jacobian
%     end
%     sigma = sigma + (J'*J)^-1*(J'*r); % update sigma
%     it = it+1; % update iteration number
%     if it == 100000 % prevent infinite loop
%         break
%     end
% end
% % plot reconstructed sigma and true sigma
figure
surf(X,Y,sigmat)
title('True Conductivity')
% figure
% surf(X,Y,sigma)
% title('Reconstructed Conductivity')