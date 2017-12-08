%{
Tutorial for running a kinetochore microtubule binding simulation
%}
clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% choose parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kinetochore parameters
num_hec1 = 10;                % number of hec1 proteins in kinetochore
tether_length = 100e-9;        % length of tether attaching hec1 to kinetochore (m)
binding_distance = 5e-9;      % distance at which hec1 binds to microtubule (m)
hec1_step = 6e-9;              % length of each step taken by hec1 in random walk (m)
prob_on_dephos = 3.1968e-4; 
prob_on_phos = 0.0024;
prob_off_dephos = 4.5e-5;
prob_off_phos = 1.8e-5;
prob_bind=[ones(num_hec1,1)*prob_on_dephos,ones(num_hec1,1)*prob_on_phos].*100;
prob_unbind=[ones(num_hec1,1)*prob_off_dephos,ones(num_hec1,1)*prob_off_phos].*100;
    % above parameters govern probabilities of binding to microtubule given
    % phosphorylation state of hec1 molecule
prob_phos= 0.5; 
prob_dephos= 1;
    % n.b. there is no good data on phosphorylation probabilities of hec1
    % by auroraB

% microtubule parameters
num_dimers = 20;              % number of dimers in a microtubule
dimer_length = 6e-9;          % length per tubulin dimer (6 nm)
e_params.S = 0.38e9;          % spring constant for microtubule (0.38 - 2 GPa)
e_params.B = 7e-23;           % bending rigidity for microtubule (7e-25 - 7e-23 Nm^2)
e_params.k = 1e-9;            % resting spring length between substrate and dimer
e_params.theta = 23*(pi/180); % preferred angle for gdp (dephos) tubulin
mt_phosphor_params = [0.1125, 0.01125]; % microtubule prob [phos, dephos]


% simulation parameters
num_time_steps = 1000;          % number of time steps to run simulation
time_step = 1.8e-4; %sec
time = (1 : 1 : num_time_steps).*time_step;
num_runs = 50;                % number of runs over which to avg sim
% max_pos_delta = 0.5;          % maximum change in dimer position at each time step
% monitor_minimization = 0;     % binary determining whether to plotprogress of minimization
save_rep = 0; % set to 1 if representative 4-panel fig should be saved
plot_avg = 0; % set to 1 if figure with average and fit should be saved

            
%%%%%%%%%%%%%%%% initialize the kinetochore and microtubule %%%%%%%%%%%%%%%
[kinetochore, microtubule] = initialize_kmt(num_time_steps, num_hec1,...
    tether_length, num_dimers, dimer_length, mt_phosphor_params, e_params);


%%%%%%%%%%%%%%%%%%%%%%%%% run the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameter to vary
S = linspace(0.38e7, 0.38e10, 50);

tau_store = zeros(length(S), 1);
coeff_store = zeros(length(S), 1);

% run a number of simulations to average and fit average to see parameters
% change over time
for a = 1 : 1: length(S);
    e_params.S = S(a);
    
    for i = 1 : num_runs
        
        % let the microtubule change phosphorylation state
        microtubule.phosphorylate()
        
        % option 1: calculate the microtubule dimer positions based on
        % energy minimization dyanamics
        %microtubule.curve_min_energy(monitor_minimization)
        
        % option 2: calculate the microtubule dimer positions based on
        % curvature predicted by theory
        smooth_window = floor(num_dimers/5); % set to 1 if you want no smoothing
        microtubule.curve_theory(smooth_window)
        
        % let the kinetochore diffuse and bind and unbind from the microtubule
        kinetochore.diffuse_bind_unbind(microtubule,prob_bind, prob_unbind,...
            prob_phos, prob_dephos, binding_distance,...
            hec1_step, dimer_length)
        
        % calculate the fraction bound for each time step
        fraction_bound = kinetochore.calc_fraction_bound();
        
        %calculate the phosphorylated fraction for each timestep
        fraction_phos = kinetochore.calc_fraction_phos();
        
        % plot  results 
        if i == 1 && save_rep == 1
            plot_all(microtubule, kinetochore, time)
            fig = gcf;
            saveas(fig, 'representative', 'fig')
        end
        
        total_simulation(i).microtubule = microtubule;
        total_simulation(i).kinetochore = kinetochore;
        total_simulation(i).fraction_bound = fraction_bound;
        total_simulation(i).fraction_phos = fraction_phos;
        
        disp(strcat('completed cycle ', num2str(i)));
    end
    % average runs and store fit parameters
    fit_params = avg_frac_bound(total_simulation, time, plot_avg);
    % fit is (A(1)./(1+exp(-A(2).*(x-A(3)))))
    tau_store(a) = 1/fit_params(2);
    coeff_store(a) = fit_params(1);    
end

figure
plot(S, tau_store, '.-', 'MarkerSize', 20)
xlabel('S'); ylabel('tau'); title('fit time scale varying with S')
figure
plot(S, coeff_store, '.-', 'MarkerSize', 20)
xlabel('S'); ylabel('A'); title('fit coefficient varying with S')


