function [positions] = microtubule(params)

% input: parameters determining microtubule evolution in time 
% output: position list of all dimers in the microtubule
% function: evolve microtubule 

    t_steps = params.t_steps;
    mt_length = params.n_dimers;
    
    first_position = zeros(mt_length,2);
    for i = 1 : mt_length
        first_position(i,:) = [i, 0]; %x,y position
    end
    positions = zeros(mt_length,2,t_steps);
    positions(:, :, 1) = first_position;
    
%     % first option. a straight microtubule that stays straight
%     figure
%     for i = 2 : t_steps
%         positions(:, :, i) = positions(:, :, 1);
%         plot(positions(:, 1, i), positions(:, 2, i))
%         hold on
%     end
%     title('option 1. straight microtubule that stays straight over time.')
    
    % second option. a straight microtubule that bends over time.
    max_bend = mt_length / 5;
    figure
    for i = 2 : t_steps
        
        for j = 1 : mt_length
            positions(j, 2, i) = positions(j, 2, i-1) + max_bend/j;
            positions(j, 1, i) = positions(j, 1, i-1);
        end
        
        plot(positions(:, 1, i), positions(:, 2, i))
        hold on
    end
    title('option 2. straight microtubule that bends over time.')
    
    % third option. get more creative and add some stochasticity!
    
   
end