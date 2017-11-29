function [full_minimization_function] = minimizer_target(constants, num_dimers, phos_state)

% pull out constants
S = constants.S;  % spring const
B = constants.B; % bending const
k = constants.k; % resting spring length
theta = constants.theta; % dephosphor angle

% define x and y in one large p vector
%p(1:num_dimers) -> x values; p(num_dimers+1, 2*num_dimers) -> y values

% create minimization function
spring_total = ''; bend_total = '';
for i = 1 : num_dimers
    spring_part = strcat('0.5*(',num2str(S),')*(x(',num2str(num_dimers+i),')-',num2str(k),')');
    spring_total = strcat(spring_total,spring_part,'+');
end
for i = 2 : num_dimers-1    
    alpha = strcat('atan(abs(x(',num2str(num_dimers+i),')-x(',num2str(num_dimers+i-1),'))/abs(x(',num2str(i),')-x(',num2str(i-1),')))');
    beta = strcat('atan(abs(x(',num2str(num_dimers+i+1),')-x(',num2str(num_dimers+i),'))/abs(x(',num2str(i+1),')-x(',num2str(i),')))');
    dimer_theta = strcat('pi-(',alpha,'+',beta,')');
    if phos_state(i) == 0
        adjusted_angle = strcat(dimer_theta,'-',num2str(theta));
    else
        adjusted_angle = dimer_theta;
    end
    bend_part = strcat('-(',num2str(B),')*cos(',adjusted_angle,')');
    bend_total = strcat(bend_total,bend_part,'+');
end
spring_total = spring_total(1:end-1);
bend_total = bend_total(1:end-1);
energy_total = strcat('(',spring_total,') + (',bend_total, ')');
 
%full_minimization_string = strcat(full_var_string,energy_total);
full_minimization_string = strcat('@(x)',energy_total);
full_minimization_function = str2func(full_minimization_string);
end









%                     total_e= 0;
%                     for i = 2 : (size(positions, 2)-1)
%                         spring_energy = 0.5*S*(positions(2,i) - k);
% 
%                         alpha = atan( abs(temp_positions(2,i)-temp_positions(2,i-1))/abs(temp_positions(1,i)-temp_positions(1,i-1)) );
%                         beta = atan( abs(temp_positions(2,i+1)-temp_positions(2,i))/abs(temp_positions(1,i+1)-temp_positions(1,i)) );
%                         dimer_theta = pi - (alpha+beta);
%                         adjusted_angle = dimer_theta - (1-phos_state(i))*theta;
%                         bending_energy = -B*cos(adjusted_angle);
% 
%                         dimer_energy = spring_energy + bending_energy;
%                         total_e = total_e + dimer_energy;
%                     end
