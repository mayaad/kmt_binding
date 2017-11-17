
classdef Kinetochore < handle
    %{
    
    Stores data associated with the a kinetochore and contains functions
    associated with kinetochore dynamics
    
    Parameters
    ----------
    hec1_positions: 3d matrix
        positions of the hec1 proteins of the kinetochore. Has shape 
        (3, number of hec1 proteins, number of time steps)
    hec1_bound: 2d matrix 
        keeps track of whether hec1 proteins are bound or unbound to the 
        microtubule. Has shape (number of hec1 proteins, number of time steps)
    tether_length: double
        length of tether attaching hec1 proteins to the kinetochore
    
    %}

    properties
        hec1_positions
        hec1_bound
        tether_length
    end
    methods
        function obj = Kinetochore(hec1_positions, hec1_bound, tether_length)
            % constructor function
            obj.hec1_positions = hec1_positions;
            obj.hec1_bound = hec1_bound;
            obj.tether_length = tether_length;
        end
        
        function diffuse(obj)
            % generates positions of tethered random walk
            
            % get some parameters from position matries
            num_time_steps = size(obj.hec1_positions,3);
            num_hec1 = size(obj.hec1_positions,2);
            
            % generate random numbers for random walk
            rand_steps = rand(3, num_hec1, num_time_steps);
            
            % TODO: get rid of this loop
            for time=2:num_time_steps
                % generate the random walk
                % TODO: replace +/- 1 with a step size
                rand_step_time = rand_steps(:,:,time);
                rand_step_time(rand_step_time<0.5) = -1;
                rand_step_time(rand_step_time>=0.5) = +1;              
                obj.hec1_positions(:,:,time) = obj.hec1_positions(:,:,time-1)+...
                    rand_step_time;
                
                % If the hec1 is outside the tether, bring it back to where it
                % was the timestep before
                distance = obj.hec1_positions(1,:,time).^2 + ...
                    obj.hec1_positions(2,:,time).^2 +...
                    obj.hec1_positions(3,:,time).^2;
                ind = find(distance > obj.tether_length);
                obj.hec1_positions(:,ind,time) = obj.hec1_positions(:,ind,time-1);             
            end            
        end
        
        function bind_unbind(obj, microtubule)
            % TODO: write this function, may have to combine with
            % diffuse() or call diffuse() in this function
            
            
        end
        function fraction_bound = calc_fraction_bound(obj)
                %{
    
                Calculates a bound fraction for the kinetochore object

                Parameters
                ----------
                kinetochore: Kinetochore object
                    object of the Kinetochore class 

                Returns
                -------
                fraction_bound: double
                    fraction of hec1 proteins that are bound to a microtubule

                %}
                fraction_bound = sum(obj.hec1_bound)/length(obj.hec1_bound);
            
        end
        
    end
end



