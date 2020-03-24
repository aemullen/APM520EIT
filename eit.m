%% Electrical Impedence Tomography
% Roberto Alvarez, Bechir Amdouni, Andrew Mullen, Trevor Reckell, Shan Zhu

% Outline
% 1 - Add noise to boundary conditions
% 2 - Guess sigma zero (initial value for conductivity)
% 3 - Find u0 (voltage) using forward model
% 4 - Calculate residual r_0 (r_i = f(sigma i, ui) - yi )
% 5 - Using u0, update sigma:
%         sigma_i = sigma_0 + ((J'J)^-1) * J'  * r_i
% 6 - Find u1 using forward model
% 7 - Find r1 and so on

%% Inverse problem
% N = 50;

function [u, sigma] = eit(N, tol, maxiter)
% Gauss Newton Method
u = zeros(N+1, N+1);
sigma = zeros(N+1, N+1);
residual = 0;
iter = 0;
% step one: boudary conditions
u(N+1, :) = normrnd(1, 0.09, N+1, 1);
u(1, :) = normrnd(1, 0.09, N+1, 1);
u(:, 1) = normrnd(1, 0.09, N+1, 1);
u(:, N+1) = normrnd(1, 0.09, N+1, 1);
% step two: initial guess for sigma
for i = 1:N+1
    for j = 1:N+1
        sigma(i,j) = 1;
    end
end

while residual > tol
% step three: solve for the voltage using forward model
u_sol = forward(N, sigma, u);
% step four: compute residual
residual = abs(u_sol - u);
% step five: update sigma
J = jac(residual);
sigma = sigma + (J'*J)^(-1)*J'*residual;
if iter == maxiter
    break
end
iter = iter + 1;
end
end
% figure
% surf(x,y,u);
% shading flat;
% xlabel('x','FontSize',24); 
% ylabel('y','FontSize',24); 
% zlabel('Potential','FontSize',24);
% 
% figure
% hh = pcolor(x,y,u);
% set(hh,'edgecolor','none','facecolor','interp');
% axis equal;
% axis off;
% set(gca,'fontsize',24);
% title('Potential');
% colormap jet;
% colorbar;

%% Solve the forward problem
function u = forward(N, sigma, u)
% Solve the 2D forward model,
% 
%   div(sigma grad(u)) = 0,
%
% on the unit square with N dx by N dy, with dx = dy = h. 
% Gridpoints are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% Uses Gauss-Seidel iteration.

h = 1/N;
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;

sum = 0;
for i = 2:N
    for j = 2:N
        A = sigma(i,j)*u(i,j)+sigma(i+1,j)*u(i+1,j)+...
            sigma(i-1,j)*u(i-1,j)+sigma(i,j+1)*u(i,j+1)+...
            sigma(i,j-1)*u(i,j-1); 
        sum = sum + abs(A);
    end
end
normresidual = 0.5;
EPSILON = 10^-5*h^2*normresidual;
while normresidual > EPSILON
    sum = 0;
    for i = 2:N
        for j = 2:N
            A = -4*sigma(i,j)*u(i,j)+sigma(i+1,j)*u(i+1,j)+...
                sigma(i-1,j)*u(i-1,j)+sigma(i,j+1)*u(i,j+1)+...
                sigma(i,j-1)*u(i,j-1);
            sum = sum + abs(A);
            u(i,j) = u(i,j) + A/4;  
        end
    end
    normresidual = sum;
end
end

function J = jac(residual)
% Calculate Jacobian with respect to spatial derivatives.
J = zeros(size(residual));
end