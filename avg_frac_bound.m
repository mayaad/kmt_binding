function [A_fit] = avg_frac_bound(results_struct, time, plot_avg, plot_fit)

% deal with averaging all simulation results
% average fraction bound over time
average_fraction_bound = zeros(length(results_struct(1).fraction_bound), 1);
for j = 1 : length(average_fraction_bound)
    total = 0;
    for i = 1 : length(results_struct)
        total = total+results_struct(i).fraction_bound(j);
    end
    average_fraction_bound(j) = total / length(results_struct);
end

if plot_avg == 1   
    figure
    for i = 1 : length(results_struct)
        hold on
        plot(results_struct(i).fraction_bound, 'b')
    end
    plot(average_fraction_bound, 'r', 'LineWidth', 6)
    xlabel('time step'); ylabel('fraction bound')
    title('average hec1 fraction bound to microtubule')
    fig = gcf;
    saveas(fig, 'avg_frac_bound', 'fig')
end

% % fit average to a saturating exponential to pull out a time scale
%logistic = @(A, x)((A(1).*exp(A(2)*(x-A(3))./((1+exp(-A(2).*(x-A(3)))));
logistic = @(A,x)(A(1).*(1-exp(-A(2).*x)));
A0 = [0.2, 1000]; % initial guess for parameters
A_fit = nlinfit(time', average_fraction_bound, logistic, A0);
fit_fxn = feval(logistic, A_fit, time);

if plot_fit == 1
    figure
    for i = 1 : length(results_struct)
        hold on
        plot(time, results_struct(i).fraction_bound, 'b')
    end
    plot(time, average_fraction_bound, 'r', 'LineWidth', 4)
    plot(time, fit_fxn, 'g', 'LineWidth', 4)
    xlabel('time step'); ylabel('fraction bound')
    title('blue: individual traces; red: average; green: fit')
end
end

