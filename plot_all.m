function plot_all(microtubule, kinetochore)

subplot(2,2,1)
% plots the dimer positions over time

% parameters pulled from position matrix
num_time_steps = size(microtubule.dimer_positions,3);

% plot
hold on
for time_step=1:num_time_steps
    plot(microtubule.dimer_positions(1,:,time_step), microtubule.dimer_positions(2,:,time_step),'.', 'MarkerSize',20)
    %plot(microtubule.dimer_positions(1,:,time_step), microtubule.dimer_positions(2,:,time_step))
end
xlabel('x-position')
ylabel('y-position')
title('dimer positions')
hold off





subplot(2,2,2)
% plot microtubule phosphorylation state over time

% parameters pulled from position matrix
num_time_steps = size(microtubule.dimer_positions,3);

hold on
for time_step = 1 : num_time_steps
    plot(microtubule.phos_state(:, :, time_step))
end
xlabel('x-position')
ylabel('phosphorylation state (0-GDP, 1-GTP)')
title('microtubule phosphorylation state')
hold off



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
    plot3(x, y, z)
    hold on
end
title('hec1 trajectories')
xlabel('x')
ylabel('y')
zlabel('z')
view(45,20)
hold off





subplot(2,2,4)
% plot the fraction bound over time

% calculate fraction bound
fraction_bound = kinetochore.calc_fraction_bound();

% plot
plot(fraction_bound)
xlabel('time')
ylabel('bound fraction')
title('bound fraction over time')


end