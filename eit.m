%% APM 520 EIT
% EIT group

%% Inverse problem
% Gauss Newton Method
N = 3;
u = zeros(N+1, N+1);
% step one: boudary conditions
u(N+1,:) = normrnd(1,0.09,N+1,1);
u(1,:)   = normrnd(1,0.09,N+1,1);
u(:,1)   = normrnd(1,0.09,N+1,1);
u(:,N+1) = normrnd(1,0.09,N+1,1);
% step two: guess initial conductivity (sigma0)
sigma = zeros(N*N, N*N);
<<<<<<< HEAD
% step one: initial guess sigma 0
for i= 2:N*N
    for j = 2:N*N
=======
for i= 1:N*N
    for j = 1:N*N
>>>>>>> 3f8fd622f1fd385c292eeea7e1f1433d056d8f2d
        sigma(i,j) = 1;
    end
end
% step three: find voltage (u0) by solving forward problem
[u_sol, residual] = Forward(N,sigma,u);
% step four: compute residual (r_i = f(sigma_i, u_i) - y_i)
r0 = norm(u_sol - u);
% step five: update sigma (sigma_i+1 = sigma_i + (J'J)^(-1)*J'*r_i)
J = [u_sol(1), u_sol(2), 0, u_sol(4),0,0,0,0,0;
    u_sol(1), u_sol(2), u_sol(3), 0, u_sol(5), 0,0,0,0;
    0, u_sol(2), u_sol(3),0,0, u_sol(6), 0,0,0;
    u_sol(1),0,0, u_sol(4), u_sol(5),0,  u_sol(7), 0,0;
    0,u_sol(2), 0, u_sol(4), u_sol(5), u_sol(6),0,u_sol(8),0;
    0,0,u_sol(3), 0, u_sol(5), u_sol(6), 0,0,u_sol(9);
    0,0,0,u_sol(4), 0, 0, u_sol(7), u_sol(8), 0;
    0,0,0,0,u_sol(5),0,  u_sol(7), u_sol(8), u_sol(9);
    0, 0, 0, 0, 0, u_sol(6), 0, u_sol(8), u_sol(9)];
sigma1 = sigma + (J'*J)^-1 *J'*r0;

% plot
k = (1:N*N)';
x = (k-1)/N; y = (k-1)/N;
figure
%clf
surf(x,y,sigma1)

%%
function [u_s, residual] = Forward(N, sigma,u)
%laplace1(N) solves the 2D Laplace equation 
%        u_xx + u_yy = 0 
%    on the unit square with N dx by N dy, with dx = dy = h. 
%    Gridpoints are labeled i = 1,...,N+1 and j = 1,...,N+1.
%    In u(i,j), i labels x and j labels y.
%    BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 1.
%    Uses Gauss-Seidel iteration.

h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;

sum = 0;
for i = 2:N
    for j = 2:N
        residual = sigma(i,j)*u(i,j)+sigma(i+1,j)*u(i+1,j)+sigma(i-1,j)*u(i-1,j)+sigma(i,j+1)*u(i,j+1)+sigma(i,j-1)*u(i,j-1); 
        sum = sum + abs(residual);
    end
end
normresidual = 0.5; % could divide by (N-1)^2 here and below in line 39
EPSILON = 10^-5*h^2*normresidual;
iter = 0;
while normresidual > EPSILON
    iter = iter + 1;
    sum = 0;
    for i = 2:N
        for j = 2:N
            residual = -4*sigma(i,j)*u(i,j)+sigma(i+1,j)*u(i+1,j)+sigma(i-1,j)*u(i-1,j)+sigma(i,j+1)*u(i,j+1)+sigma(i,j-1)*u(i,j-1);
            sum = sum + abs(residual);
            u(i,j) = u(i,j) + residual/4;  
        end
    end
    normresidual = sum;
end
iterations = iter;

u_s = u;
end