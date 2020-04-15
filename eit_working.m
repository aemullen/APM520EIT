% function eit()

N = 9;

sigmat = ((cos(1:N)) + sin(1:N)).^2;

b = normrnd(1.0, 0.005, N, 1)*1e-3;
A = diffMat(sigmat);
u0 = A\b;

sigma = 0.8*ones(N,1);
A = diffMat(sigma);
u_sol = A\b;
r = abs(u_sol - u0);

J = [-4*u_sol(1)+u_sol(2)+u_sol(4),-u_sol(1)+u_sol(2),0,-u_sol(1)+u_sol(4),0,0,0,0,0;
    u_sol(1),-4*u_sol(2)+u_sol(1)+u_sol(3) + u_sol(5),-u_sol(2)+u_sol(3),0,-u_sol(2)+u_sol(4),0,0,0,0;
    0,u_sol(2)-u_sol(3),-4*u_sol(3)+u_sol(2)+u_sol(6),0,0,-u_sol(3)+u_sol(6),0,0,0;
    u_sol(1)-u_sol(4),0,0,-4*u_sol(4)+u_sol(1)+u_sol(5)+u_sol(7),-u_sol(4)+u_sol(5),0,-u_sol(4)+u_sol(7),0,0;
    0,u_sol(2)-u_sol(5),0,u_sol(4)-u_sol(5),-4*u_sol(5)+u_sol(2)+u_sol(4)+u_sol(6)+u_sol(8),-u_sol(5)+u_sol(6),0,-u_sol(5)+u_sol(8),0;
    0,0,u_sol(3)-u_sol(6),0,u_sol(5)-u_sol(6),-4*u_sol(6)+u_sol(3)+u_sol(5),0,0,-u_sol(6)+u_sol(9);
    0,0,0,u_sol(4)-u_sol(7),0,0,-4*u_sol(7)+u_sol(4)+u_sol(8),-u_sol(7)+u_sol(8),0;
    0,0,0,0,u_sol(5)-u_sol(8),0,u_sol(7)-u_sol(8),-4*u_sol(8)+u_sol(5)+u_sol(7)+u_sol(9),-u_sol(8)+u_sol(9);
    0,0,0,0,0,u_sol(6)-u_sol(9),0,u_sol(8)-u_sol(9),-4*u_sol(9)+u_sol(6)+u_sol(8)];

sigma = sigma + (J'*J)\(J'*r);
% end

function A = diffMat(sigma)
A = [-(sigma(2)+sigma(4)+4*sigma(1)), (sigma(1)+sigma(2)), 0, (sigma(4)+sigma(1)), 0, 0, 0, 0, 0;
     (sigma(1)+sigma(2)), -(sigma(3)+sigma(5)+4*sigma(2)), (sigma(3)+sigma(2)), 0, (sigma(5)+sigma(2)), 0, 0, 0, 0;
     0, (sigma(2)+sigma(3)), -(sigma(2)+sigma(6)+4*sigma(3)), 0, 0, (sigma(3)+sigma(6)), 0, 0, 0;
     (sigma(1)+sigma(4)), 0, 0, -(sigma(1)+sigma(5)+sigma(7)+4*sigma(4)), (sigma(4)+sigma(5)), 0, (sigma(4)+sigma(7)), 0, 0;
     0, (sigma(2)+sigma(5)), 0, (sigma(4)+sigma(5)), -(sigma(2)+sigma(4)+sigma(6)+sigma(8)+4*sigma(5)), (sigma(5)+sigma(6)), 0, (sigma(5)+sigma(8)), 0;
     0, 0, (sigma(3)+sigma(6)), 0, (sigma(5)+sigma(6)), -(sigma(3)+sigma(5)+sigma(9)+4*sigma(6)), 0, 0, (sigma(6)+sigma(9));
     0, 0, 0, (sigma(4)+sigma(7)), 0, 0, -(sigma(4)+sigma(8)+4*sigma(7)), (sigma(7)+sigma(8)), 0;
     0, 0, 0, 0, (sigma(5)+sigma(8)), 0, (sigma(7)+sigma(8)), -(sigma(5)+sigma(7)+sigma(9)+4*sigma(8)), (sigma(8)+sigma(9));
     0, 0, 0, 0, 0, (sigma(6)+sigma(9)), 0, (sigma(8)+sigma(9)), -(sigma(6)+sigma(8)+4*sigma(9))];
end