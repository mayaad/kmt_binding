function [y_pos] = bending_theory(phos_state, e_params, x_pos, dimer_length, smooth_window)
%{

returns y positions of a bending microtubule predicted by theory.

Parameters
----------
phos_state: vector
    phosphorylation state at a given time for all dimers of microtubule
e_params: structure containing parameters for microtubule energy 
    minimization:
    1) S = spring constant
    2) B = bending rigidity 
    3) k = resting spring length
    4) theta = preferred angle for dephosphorylated tubulin
x_pos: vector
    x-coordinates of postitions of dimers at a given time
dimer_length: float
    length between centers of masses of dimers
smooth-window: int (odd)
    length if window accross which to do the smoothing

Returns
-------
y_pos: vector
    y-coordinates of dimers in bending microtubule at a given time

%}


B = e_params.B;
S = e_params.S;
a = dimer_length;
x_pos = x_pos*a; % convert x_position from dimer number to units of dimer length

num_dimers = length(x_pos);
theta = zeros(1,num_dimers);
theta(phos_state == 0) = 23*pi/180;

kappa = theta/a;
q = (S/B)^(1/4)/sqrt(2);
y_pos = kappa./(2*q^2*(1+a*S/(4*B*q^3))).*exp(-q*x_pos).*cos(q*x_pos)- kappa./(2*q^2).*exp(-q*x_pos).*sin(q*x_pos);
y_pos = smooth(y_pos, smooth_window);
end