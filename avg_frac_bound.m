function [k, A] = avg_frac_bound(results_struct, plot_var)

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

if plot_var == 1   
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
    close all
end

% fit average to a saturating exponential to pull out a time scale
time = 1 : 1 : length(average_fraction_bound);
fit_result = fit(time', average_fraction_bound, 'exp1');

k = fit_result.b; A = fit_result.a;

end

