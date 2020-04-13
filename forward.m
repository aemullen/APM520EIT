function u = forward(N,sigma,b1,b2,b3,b4,tol,maxiter)
% forward solves the forward EIT model:
% 
%        div(sigma(grad(u))) = 0
% 
%  iteratively on the (N-1) x (N-1) unit square. 
%  Gridpoints are labeled i = 1,...,N+1 and j = 1,...,N+1.
%  In u(i,j), i labels x and j labels y.
% 
% Parameters: N - number of grid points
%             sigma - conductivity
%             b1 - left boundary
%             b2 - top boundary
%             b3 - bottom boundary
%             b4 - right boundary
%             tol - tolerance for solver
%             maxiter - maximum number of iterations
%
% setup domain
u = zeros(N+1);
% assign boundary conditions
u(2:N-1,1) = b1;
u(1,2:N-1) = b2;
u(N+1,2:N-1) = b3;
u(2:N-1,N+1) = b4;
normresidual = 1;
iter = 0;
while tol > normresidual
    iter = iter + 1;
    sum = 0;
    for i = 2:N
        for j = 2:N
            residual = -(sigma(i+1,j)+sigma(i-1,j)+sigma(i,j+1)+sigma(i,j-1)+4*sigma(i,j))*u(i,j)...
                +(sigma(i+1,j)+sigma(i,j))*u(i+1,j)...
                +(sigma(i-1,j)+sigma(i,j))*u(i-1,j)+...
                +(sigma(i,j+1)+sigma(i,j))*u(i,j+1)...
                +(sigma(i,j-1)+sigma(i,j))*u(i,j-1); 
            sum = sum + abs(residual);
            u(i,j) = u(i,j) + residual/(sigma(i+1,j) + sigma(i-1,j) + sigma(i,j+1)+ sigma(i,j-1) + 4*sigma(i,j));   
        end
    end
    normresidual = sum;
    if iter == maxiter
        break
    end
end
% return a vector
u = reshape(u(2:N,2:N),(N-1)*(N-1),1);
end