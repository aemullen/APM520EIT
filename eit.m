% 1- Add noise to boundary conditions
% 2- Guess sigma zero (initial value for conductivity)
% 3- Find u0 (voltage) using forward model
% 4- Calculate r0 (residual) (ri = f(sigma i, ui) - yi )
% 5- Using u0, update sigma:
%         sigma 1 = sigma 0 + ((J'J)^-1) * J'  * r0
% 6 - Find u1 using forward model
% 7- Find r1 and so on

%% Inverse problem
N = 50;
u = zeros(N+1, N+1);
% Boudary conditions
k = (1:N+1)';
n = normrnd(1,0.09,N+1,1);
u(N+1,k) = normrnd(1,0.09,N+1,1); % nonzero BC 
u(1,k) = normrnd(1,0.09,N+1,1);
u(k,1) = normrnd(1,0.09,N+1,1); % nonzero BC 
u(k,N+1) = normrnd(1,0.09,N+1,1);


% Gauss Newton Method


sigma = zeros(N+1, N+1);
% step one: initial guess sigma 0
for i= 1:N+1
    for j = 1:N+1
        sigma(i,j) = 1;
        %u = Forward(N,sigma);
        
    end
     
end
u_sol = Forward(N,sigma,u);

% Step 4

r1 = abs(u - u_sol);
% sigma1 = sigma + ((J'J)^-1) * J'  * r0     % Need to find the Jacobian

%u1 = Forward(N, sigma1, u_sol);

%%
function u = Forward(N, sigma,u)
%laplace1(N) solves the 2D Laplace equation 
%        u_xx + u_yy = 0 
%    on the unit square with N dx by N dy, with dx = dy = h. 
%    Gridpoints are labeled i = 1,...,N+1 and j = 1,...,N+1.
%    In u(i,j), i labels x and j labels y.
%    BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 1.
%    Uses Gauss-Seidel iteration.

%tic
h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
% u = zeros(N+1,N+1);
% u(N+1,k) = 1; % nonzero BC 
% u(1,k) = 1;
% u(k,1) = 1; % nonzero BC 
% u(k,N+1) = 1;

sum = 0;
%sigma = 0.39;
% for i = 2:N
%     for j = 2:N
%         residual = sigma(i,j)*u(i,j)+sigma(i+1,j)*u(i+1,j)+sigma(i-1,j)*u(i-1,j)+sigma(i,j+1)*u(i,j+1)+sigma(i,j-1)*u(i,j-1); 
%         sum = sum + abs(residual);
%     end
% end
normresidual = 0.5; % could divide by (N-1)^2 here and below in line 39
EPSILON = 10^-5*h^2*normresidual;
iter = 0;
while normresidual > EPSILON
    iter = iter + 1;
    sum = 0;
    for i = 2:N
        for j = 2:N
            A = -(1/2)*(4*sigma(i,j) +sigma(i+1,j)+ sigma(i-1,j)+...
            sigma(i,j+1)+sigma(i,j-1)) *u(i,j)+sigma(i+1,j)*u(i+1,j)+...
            sigma(i-1,j)*u(i-1,j)+sigma(i,j+1)*u(i,j+1)+sigma(i,j-1)*u(i,j-1);
            sum = sum + abs(A);
            u(i,j) = u(i,j) + A/((1/2)*(4*sigma(i,j) +sigma(i+1,j)+...
            sigma(i-1,j)+ sigma(i,j+1)+sigma(i,j-1)));  
        end
    end
    normresidual = sum;
end
iterations = iter
%toc

%u = u'; % transpose to plot

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

end

%function J = Jacob(A,u)
