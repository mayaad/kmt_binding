%{
Tutorial for running a kinetochore microtubule binding simulation
%}

% choose initial parameters
num_time_steps = 5;
num_hec1 = 4;
tether_length = 10;
prob_bind=0.4; % related to k_bind
prob_unbind=1e-4; % related to k_unbind
binding_distance = 2.5;
num_dimers = 10;
dimer_length = 1;

% initialize the kinetochore and microtubule
[kinetochore, microtubule] = initialize_kmt(num_time_steps, num_hec1,...
    tether_length, num_dimers, dimer_length);


% let the microtubule curve
microtubule.curve()

% plot the microtubule positions over time
figure
hold on
for time_step=1:num_time_steps
    plot(microtubule.dimer_positions(1,:,time_step), microtubule.dimer_positions(2,:,time_step))
end
xlabel('x-position')
ylabel('y-position')
hold off

% let the kinetochore diffuse
kinetochore.diffuse()

% plot the trajectories of the hec1 proteins
figure
for hec1=1:num_hec1
    x = reshape(kinetochore.hec1_positions(1,hec1,:),[1,5]);
    y = reshape(kinetochore.hec1_positions(2,hec1,:),[1,5]);
    z = reshape(kinetochore.hec1_positions(3,hec1,:),[1,5]);
    plot3(x, y, z)
    hold on
end

% TODO: let the kinetochore bind and unbind

% calculate the fraction bound for each time step
fraction_bound = kinetochore.calc_fraction_bound();

% plot the fraction bound over time
figure
plot(num_time_steps,fraction_bound)
xlabel('time')
ylabel('bound fraction')

