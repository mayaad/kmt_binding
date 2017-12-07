
classdef Microtubule < handle
    %{
    
    stores data associated with a microtubule
    
    Parameters
    ----------
    dimer_positions: 3d matrix
        positions of the dimers of the microtubule. Has shape
        (2, number of dimers, number of time steps). The 2 comes from the
        fact that the microtubule position can be defined by only two
        dimensions, because microtubule bending takes place in one plane
    dimer_length: float
        length of the the dimers makting up the microtubule
    phosphor: structure containing
        1) phos_state: vector
            phosphorylation state (0 GDP, 1 GTP) of each dimer
        2) params: vector
            probabilities: [p(phos|dephos), p(dephos|phos)]
    energy_params: structure containing parameters for microtubule energy
        minimization:
        1) S = spring constant
        2) B = bending rigidity
        3) k = resting spring length
        4) theta = preferred angle for dephosphorylated tubulin
    %}
    
    properties
        dimer_positions
        dimer_length
        phos_state
        phos_params
        e_params
    end
    methods
        function obj = Microtubule(dimer_positions, dimer_length, phosphor, energy_params)
            % constructor function
            obj.dimer_positions = dimer_positions;
            obj.dimer_length = dimer_length;
            obj.phos_state = phosphor.phos_state;
            obj.phos_params = phosphor.params;
            obj.e_params = energy_params;
        end
        
        function curve(obj)
            % curves the microtubule over time
            
            % get some parameters from the positions matrices
            num_dimers = size(obj.dimer_positions,2);
            num_time_steps = size(obj.dimer_positions,3);
            
            max_bend = num_dimers/5;
            
            % dimer y-positions are shifted up at each timestep according
            % to the curvature defined by max_bend
            obj.dimer_positions(2,:,2:end) = obj.dimer_positions(2,:,1:end-1)...
                + cumsum(reshape(repmat(max_bend./(1:obj.dimer_length:num_dimers),...
                num_time_steps-1,1)',[1,num_dimers,num_time_steps-1]),3);
            
            % dimer x-positions are shifted back one timestep
            obj.dimer_positions(1,:,2:end) = obj.dimer_positions(1,:,1:end-1);
            
        end
        
        function phosphorylate(obj)
            % updates phosphorylation state of microtubule dimers over time
            % 1- GTPtubulin (phosphorylated), 0- GDPtubulin (dephosphorylated)
            
            % phosphorylation probabilities
            prob_phos = obj.phos_params(1); % p(phos|dephos)
            prob_dephos = obj.phos_params(2); % p(dephos|phos)
            
            % parameters pulled from position matrix
            num_dimers = size(obj.dimer_positions,2);
            num_time_steps = size(obj.dimer_positions,3);
            
            for i = 2 : num_time_steps
                for j = 1 : num_dimers
                    change_variable = rand();
                    
                    if obj.phos_state(1, j, i-1) == 1 % dimer was GTP
                        if change_variable > (1-prob_dephos)
                            obj.phos_state(1, j, i) = 0; % dimer becomes GDP
                        else
                            obj.phos_state(1, j, i) = 1; % dimer stays GTP
                        end
                    elseif obj.phos_state(1, j, i-1) == 0 % dimer was GDP
                        if change_variable > (1-prob_phos)
                            obj.phos_state(1, j, i) = 1; % dimer becomes GTP
                        else
                            obj.phos_state(1, j, i) = 0; % dimer stays GDP
                        end
                    end
                end
            end
        end
        
        function curve_min_energy(microtubule, monitor_minimization)
            
            % get some parameters from the positions matrices
            num_dimers = size(microtubule.dimer_positions,2);
            num_time_steps = size(microtubule.dimer_positions,3);
            
            for tstep = 2 : num_time_steps

                phos_state = microtubule.phos_state(:, :, tstep);
                energy_function = minimizer_target(microtubule.e_params, num_dimers, phos_state);

                % create input vector of initial guess for positions (use previous
                % positions) -> need to reshape initial guess vector to be 1 row with
                % 2*num_dimers columns
                % guess(1:num_dimers) -> x values; guess(num_dimers+1, 2*num_dimers) -> y values
                init_guess = [microtubule.dimer_positions(1,:, tstep-1), microtubule.dimer_positions(2,:, tstep-1)];

                if monitor_minimization == 1
                   options = optimset('PlotFcns',@optimplotfval); 
                   [pos, energy] = fminsearch(energy_function, init_guess, options);
                else
                   [pos, energy] = fminsearch(energy_function, init_guess);
                end

            %     % update dimer positions and enforce a maximum change in position
            %     for j = 1 : num_dimers   
            %         if abs(pos(j)-init_guess(j)) > max_pos_delta
            %             pos(j) = init_guess(j) + max_pos_delta * (pos(j)/abs(pos(j)));
            %         end
            %         if abs(pos(j+num_dimers)-init_guess(j+num_dimers)) > max_pos_delta
            %             pos(j+num_dimers) = init_guess(j+num_dimers) + max_pos_delta * (pos(j+num_dimers)/abs(pos(j+num_dimers)));
            %         end
            %     end
            %     
            %     microtubule.dimer_positions(:,:,tstep) = [pos(1:num_dimers); pos(num_dimers+1:end)];
                microtubule.dimer_positions(:,:,tstep) = [microtubule.dimer_positions(1,:,tstep); pos(num_dimers+1:end)];
            end
        end
        
        function curve_theory(microtubule, smooth_window)
            % get some parameters from the positions matrices
            num_dimers = size(microtubule.dimer_positions,2);
            num_time_steps = size(microtubule.dimer_positions,3);
            
            for tstep = 2 : num_time_steps
                % update y-positions based on theory
                microtubule.dimer_positions(2,:,tstep) = ...
                    bending_theory(microtubule.phos_state(:,:,tstep),...
                                   microtubule.e_params,...
                                   microtubule.dimer_positions(1,:,tstep),...
                                   microtubule.dimer_length, smooth_window);
            end
            
        end
        
        function plot_dimer_positions(obj)
            % plots the dimer positions over time 
            
            % parameters pulled from position matrix
            num_time_steps = size(obj.dimer_positions,3);
            
            % plot
            figure
            hold on
            for time_step=1:num_time_steps
                plot(obj.dimer_positions(1,:,time_step)*obj.dimer_length, obj.dimer_positions(2,:,time_step),'.', 'MarkerSize',20)
                %plot(obj.dimer_positions(1,:,time_step), obj.dimer_positions(2,:,time_step))
            end
            xlabel('x-position (m)')
            ylabel('y-position (m)')
            title('dimer positions')
            hold off
        end
        
        function plot_phos_state(obj) 
            % plot microtubule phosphorylation state over time
            
            % parameters pulled from position matrix
            num_time_steps = size(obj.dimer_positions,3);

            figure
            hold on
            for time_step = 1 : num_time_steps
                plot(obj.phos_state(:, :, time_step))
            end
            xlabel('x-position')
            ylabel('phosphorylation state (0-GDP, 1-GTP)')
            title('microtubule phosphorylation state')
            hold off
            
        end
        
    end
end


