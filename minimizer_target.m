function [total_energy] = minimizer_target(positions, phos_state)

S = obj.e_params.S; % spring const
B = obj.e_params.B; % bending const
k = obj.e_params.k; % resting spring length
theta = obj.e_params.theta; % dephosphor angle

total_energy = 0;
for i = 2 : (size(positions, 2)-1)
    spring_energy = 0.5*S*(positions(i, 2) - k);
    
    alpha = atan( abs(positions(i,2)-positions(i-1,2))/abs(positions(i,1)-postiions(i-1,1)) );
    beta = atan( abs(positions(i+1,2)-positions(i,2))/abs(positions(i+1,1)-postiions(i,1)) );
    dimer_theta = pi - (alpha+beta);
    adjusted_angle = dimer_theta - (1-phos_state(i))*theta;
    bending_energy = -B*cos(adjusted_angle);
    
    dimer_energy = spring_energy + bending_energy;
    total_energy = total_energy + dimer_energy;
end

end