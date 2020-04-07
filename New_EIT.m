%% Forward problem using sigmas = 1
clear all; close all; clc;
N=9;
M=ones(N,N);
d=ones(1,N)*-4;
d2=ones(1,N-1);
d3=ones(1,N^2-3);
a=diag(d)+diag(d2,1)+diag(d2,-1);
ACell = repmat({a}, 1, N);
A = blkdiag(ACell{:})+diag(d3,3)+diag(d3,-3);
u=zeros(N*N,1);
b = normrnd(1,0.5, N*N,1);
u_sol1 = A\b;

% Backward problem using Gauss Newton Method
sigma0 = .8*ones(N*N,1);
u_sol = (0.8*A)\b;
r = abs(u_sol1 - u_sol);
J = zeros(N*N,N*N);
for i = 1:N*N
    J(:,i)=u_sol(i)*0.8*A(:,i);
end
%sigma = sigma0 + (J'*J)\ (J'*r);
sigma = sigma0 + J\r;
tol = 1e-12;
it = 0;
while r >tol
    
    %sigmat = sigma;
    ut = u_sol;
    u_sol = (sigma(N)*A)\b;
    r = abs(u_sol-ut);
    for i = 1:N*N
        J(:,i)=u_sol(i)*sigma(1)*A(:,i);
    end
    %sigma = sigmat + (J'*J)\ (J'*r);
    sigma = sigma + J\r;
    it = it+1;
end      
it       
[x,y] = meshgrid(1:N,1:N);
figure;
surf(x,y,reshape(sigma,N,N))
         
        
        
