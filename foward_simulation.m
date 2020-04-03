% This experiment has 16 electrodes around the domain that pump and
% measure currents. The domain is a 40x40 square of water with a 5x5
% copper square in the center(ish). Each electrode is a cylinder of
% radius 1 centimeter and length of 5 millimeter.
N = 16; % number of electrodes
% conductivity of water from
% https://www.lenntech.com/applications/ultrapure/conductivity/water-conductivity.htm
n = 40; % mesh size
domain = ones(n)*(5.5E-6);
% conductivity of copper from
% https://www.thoughtco.com/electrical-conductivity-in-metals-2340117
domain(20:25,20:25) = 5.98e7;
% resistivity of gold from
% https://www.thoughtco.com/electrical-conductivity-in-metals-2340117
r = 2.44e-8;
% we assume circular electrodes
radius = 1e-2; % 1 centimeter
length = 5e-4; % 5 millimeters
area = 2*pi*radius;
% resistance of electrodes
% https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity
R = (radius*length)/area;
% current coming in from one electrodes
current = normrnd(1, 0.0005)*10^-3; % 1 mA
% voltage on that electrode
V = current/R;
% convert current to current density I/A
J = current/area;
% calculate electric field from Ohm's Law
% div(u) = J/sigma
u = zeros(n+1);
u(5) = V;
h = 1/n;
sum = 0;
for i = 2:n
    for j = 2:n
        residual = J/domain(i,j)-4*u(i,j)+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1);
        sum = sum + abs(residual);
    end
end
normresidual = sum; % could divide by (N-1)^2 here and below in line 39
EPSILON = 10^-5*h^2*normresidual;
iter = 0;
while normresidual > EPSILON
    iter = iter + 1;
    sum = 0;
    for i = 2:n
        for j = 2:n
            residual = J/domain(i,j)-4*u(i,j)+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1);
            sum = sum + abs(residual);
            u(i,j) = u(i,j) + residual/4;   
        end
    end
    normresidual = sum;
end
% plot
[x, y] = meshgrid((1:(n+1))/n);
figure
surf(x, y, u')