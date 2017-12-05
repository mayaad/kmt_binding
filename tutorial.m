%{
Tutorial for running a kinetochore microtubule binding simulation
%}
clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% choose parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_time_steps = 50;             % number of time steps to run simulation
num_hec1 = 10;                   % number of hec1 proteins in kinetochore
tether_length = 10e-9;           % length of tether attaching hec1 to kinetochore
prob_bind= 0.4;                  % probability related to k_bind
prob_unbind=1e-4;                % probability related to k_unbind
binding_distance = 5e-9;         % distance at which hec1 binds to microtubule
num_dimers = 20;                 % number of dimers in a microtubule
dimer_length = 1;%6e             % 6 nm per tubulin dimer (6e-9 m)
hec1_step = 0.5e-9;              % length of each step taken by hec1 in random walk
mt_phosphor_params = [0.01, 0.1];  % microtubule phosphorylation probabilities
                                 % [prob of phos, prob of dephos]
e_params.S = 0.1;%1000;% 0.38e9; % spring constant for microtubule (0.38 - 2 GPa)
e_params.B = 1; %1e-20;%7e-25;   % bending rigidity for microtubule (7e-25 - 7e-23 Nm^2)
e_params.k = 0; %10 %1e-9;       % resting spring length between substrate and dimer
e_params.theta = 23*(pi/180);    % preferred angle for gdp (dephos) tubulin

max_pos_delta = 0.5;             % maximum change in dimer position at each time step
monitor_minimization = 0;        % binary determining whether to plotprogress of minimization


%%%%%%%%%%%%%%%% initialize the kinetochore and microtubule %%%%%%%%%%%%%%%
[kinetochore, microtubule] = initialize_kmt(num_time_steps, num_hec1,...
    tether_length, num_dimers, dimer_length, mt_phosphor_params, e_params);


%%%%%%%%%%%%%%%%%%%%%%%%% run the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let the microtubule change phosphorylation state
microtubule.phosphorylate()

% option 1: calculate the microtubule dimer positions based on 
% energy minimization dyanamics 
%microtubule.curve_min_energy(monitor_minimization)

% option 2: calculate the microtubule dimer positions based on 
% curvature predicted by theory
smooth_window = floor(num_dimers/5);
microtubule.curve_theory(smooth_window)

% let the kinetochore diffuse and bind and unbind from the microtubule
kinetochore.diffuse_bind_unbind(microtubule,prob_bind, prob_unbind,...
                                binding_distance, hec1_step)

% calculate the fraction bound for each time step
fraction_bound = kinetochore.calc_fraction_bound();

%%%%%%%%%%%%%%%%%%%%%%%%%% plot  results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

microtubule.plot_dimer_positions()
microtubule.plot_phos_state()
kinetochore.plot_hec1_trajectories()
kinetochore.plot_fraction_bound()


