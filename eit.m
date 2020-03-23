% 1- Guess sigma zero (initial value for conductivity)
% 2- Find u0 (voltage) using forward model
% 3- Add noise to b
% 4- Calculate r0 (residual) (ri = f(sigma i, ui) - yi )
% 5- Using u0, update sigma:
%         sigma 1 = sigma 0 + ((J'J)^-1) * J'  * r0
% 6 - Find u1 using forward model
% 7- Find r1 and so on