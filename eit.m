% 1- Guess sigma zero (initial value for conductivity)
% 2- Find u0 (voltage) using forward model
% 3- Add noise to b
% 4- Calculate r0 (residual) (ri = f(sigma i, ui) - yi )
% 5- Using u0, update sigma:
%         sigma 1 = sigma 0 + ((J'J)^-1) * J'  * r0
% 6 - Find u1 using forward model
% 7- Find r1 and so on




function Forward(N)
%laplace1(N) solves the 2D Laplace equation 
%        u_xx + u_yy = 0 
%    on the unit square with N dx by N dy, with dx = dy = h. 
%    Gridpoints are labeled i = 1,...,N+1 and j = 1,...,N+1.
%    In u(i,j), i labels x and j labels y.
%    BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 1.
%    Uses Gauss-Seidel iteration.

tic
h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
u = zeros(N+1,N+1);
u(N+1,k) = 1; % nonzero BC u(x=1,y) = 1
% u(k,1) = 0.5; % nonzero BC u(x,y=0) = 0.5

sum = 0;
sigma = 0.39;
for i = 2:N
    for j = 2:N
        residual = -0.5*(8*sigma)*u(i,j)+sigma*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)); 
        sum = sum + abs(residual);
    end
end
normresidual = sum; % could divide by (N-1)^2 here and below in line 39
EPSILON = 10^-5*h^2*normresidual;
iter = 0;
while normresidual > EPSILON
    iter = iter + 1;
    sum = 0;
    for i = 2:N
        for j = 2:N
            residual = -4*sigma*u(i,j)+sigma*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
            sum = sum + abs(residual);
            u(i,j) = u(i,j) + residual/4;  
        end
    end
    normresidual = sum;
end
iterations = iter
toc

u = u'; % transpose to plot

figure
surf(x,y,u);
shading flat;
xlabel('x','FontSize',24); 
ylabel('y','FontSize',24); 
zlabel('Potential','FontSize',24);

figure
hh = pcolor(x,y,u);
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis off;
set(gca,'fontsize',24);
title('Potential');
colormap jet;
colorbar;

end