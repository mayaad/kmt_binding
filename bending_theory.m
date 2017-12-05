function [y_pos] = bending_theory(phos_state, e_params, x_pos, dimer_length, smooth_window)
B = e_params.B;
S = e_params.S;
a = dimer_length;

num_dimers = length(x_pos);
theta = zeros(1,num_dimers);
theta(phos_state == 0) = 23*pi/180;

kappa = theta/a;
q = (S/B)^(1/4)/sqrt(2);

y_pos = kappa./(2*q^2*(1+a*S/(4*B*q^3))).*exp(-q*x_pos)- kappa./(2*q^2).*exp(-q*x_pos).*sin(q*x_pos);
y_pos = smooth(y_pos, smooth_window);
end