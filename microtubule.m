 
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
    end
    methods
        function obj = Microtubule(dimer_positions, dimer_length)
            % constructor function
            obj.dimer_positions = dimer_positions;
            obj.dimer_length = dimer_length;
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
    end
end
    