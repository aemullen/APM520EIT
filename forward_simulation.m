function forward_simulation(N, sigma)
% forward_simulation will simulate experimental conditions for given sigma
% and return boundary condtions.
% 
% This simulation has 12 electrodes around the domain that measure currents.
% Each electrode is a cylinder of radius 1 centimeter and length of 5 millimeter.
% 1 mA of current is sent into the domain from one electrode.
% The rest of the elctrodes mearsure the resulting current at their end
% and convert to voltage using Ohm's law.
% 
% resistivity of gold from
% https://www.thoughtco.com/electrical-conductivity-in-metals-2340117
r = 2.44e-8;
% we assume circular electrodes
radius = 1e-2; % 1 centimeter
length = 5e-4; % 5 millimeters
area = 2*pi*radius;
% resistance of electrodes
% https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity
R = (r*length)/area;
% current coming in from one electrode
current = normrnd(1, 0.00005)*10^-3; % 1 mA + noise
% convert current to voltage using Ohm's Law
u = zeros(N+1);
u(2) = current*R;
% calculate electric field from div(u) = 0
iter = 0;
normresidual = 1;
while normresidual > 1e-12
    iter = iter + 1;
    sum = 0;
    for i = 2:n
        for j = 2:n
            residual = -4*u(i,j)+u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1);
            sum = sum + abs(residual);
            u(i,j) = u(i,j) + residual/4;   
        end
    end
    normresidual = sum;
    if iter == 100000
        break
    end
end
end