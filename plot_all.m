function plot_all(microtubule, kinetochore, time)

figure

subplot(2,2,1)
% plots the dimer positions over time

% parameters pulled from position matrix
num_time_steps = size(microtubule.dimer_positions,3);

% plot
hold on
for time_step=1:num_time_steps
    plot(microtubule.dimer_positions(1,:,time_step)*microtubule.dimer_length*(10^9), microtubule.dimer_positions(2,:,time_step)*(10^9),'.', 'MarkerSize',20)
end
xlabel('x-position (nm)')
ylabel('y-position (nm)')
title('microtubule shape')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
% plot microtubule phosphorylation state over time

% parameters pulled from position matrix
num_time_steps = size(microtubule.dimer_positions,3);

hold on
for time_step = 1 : num_time_steps
    plot(microtubule.dimer_positions(1,:,time_step)*microtubule.dimer_length*(10^9), microtubule.phos_state(:, :, time_step))
end
xlabel('x-position (nm)')
ylabel('phosphorylation state (0-GDP, 1-GTP)')
title('microtubule phosphorylation state')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
% plot the trajectories of the hec1 proteins

% get some parameters from position matries
num_time_steps = size(kinetochore.hec1_positions,3);
num_hec1 = size(kinetochore.hec1_positions,2);

% plot
hold on
for hec1=1:num_hec1
    x = reshape(kinetochore.hec1_positions(1,hec1,:),[1,num_time_steps]);
    y = reshape(kinetochore.hec1_positions(2,hec1,:),[1,num_time_steps]);
    z = reshape(kinetochore.hec1_positions(3,hec1,:),[1,num_time_steps]);
    plot3(x*(10^9), y*(10^9), z*(10^9))
    hold on
end
title('hec1 trajectories')
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view(45,20)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)
% plot the fraction bound over time

% calculate fraction bound
fraction_bound = kinetochore.calc_fraction_bound();
% plot
plot(time, fraction_bound)
xlabel('time (sec)')
ylabel('bound fraction')
title('bound fraction over time')


end