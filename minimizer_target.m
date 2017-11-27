function [total_energy] = minimizer_target(constants, positions, phos_state)

S = constants(1);  % spring const
B = constants(2); % bending const
k = constants(3); % resting spring length
theta = constants(4); % dephosphor angle

total_energy = 0;
for i = 2 : (size(positions, 2)-1)
    a = 3
    spring_energy = 0.5*S*(positions(2,i) - k);
    
    alpha = atan( abs(positions(2,i)-positions(2,i-1))/abs(positions(1,i)-positions(1,i-1)) );
    beta = atan( abs(positions(2,i+1)-positions(2,i))/abs(positions(1,i+1)-positions(1,i)) );
    dimer_theta = pi - (alpha+beta);
    adjusted_angle = dimer_theta - (1-phos_state(i))*theta;
    bending_energy = -B*cos(adjusted_angle);
    
    dimer_energy = spring_energy + bending_energy;
    total_energy = total_energy + dimer_energy;
end

end