 
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
    %}
    
    properties
        dimer_positions
        dimer_length
        phos_state
        phos_params
    end
    methods
        function obj = Microtubule(dimer_positions, dimer_length, phosphor)
            % constructor function
            obj.dimer_positions = dimer_positions;
            obj.dimer_length = dimer_length;
            obj.phos_state = phosphor.phos_state;
            obj.phos_params = phosphor.params;
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
                        if change_variable > prob_dephos
                            obj.phos_state(1, j, i) = 0; % dimer becomes GDP
                        end
                    elseif obj.phos_state(1, j, i-1) == 0 % dimer was GDP
                        if change_variable > prob_phos
                            obj.phos_state(1, j, i) = 0; % dimer becomes GTP
                        end
                    end
                end
            end
        end
        
    end
end
    
