
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
        hec1_phos
    end
    methods
        function obj = Kinetochore(hec1_positions, hec1_bound, tether_length, hec1_phos)
            % constructor function
            obj.hec1_positions = hec1_positions;
            obj.hec1_bound = hec1_bound;
            obj.tether_length = tether_length;
            obj.hec1_phos=hec1_phos;
        end
        
        function diffuse_bind_unbind(obj, microtubule, prob_bind, prob_unbind, prob_phos, prob_dephos, binding_distance, hec1_step, dimer_length)
            %{
            
            Calculates diffusion and binding/unbinding of hec1 proteins in
            the kinetochore. Does not return a value, but updates
            kinetochore.hec1_positions and kinetochore.hec1_bound

            Parameters
            ----------
            kinetochore: Kinetochore object
                kinetochore whose binding we are simulating
            microtubule: Microtubule object
                microtbule that the kinetochore is binding to
            prob_bind: float
                prbability of binding, related to k_bind
            prob_unbind: float
                probability of unbinding, related to k_unbind
            prob_phos: float
                prbability of phosphorylating
            prob_dephos: float
                probability of dephosphorylating
            binding_distance: float
                distance at which kinetochore can bind to microtubule
            hec1_step: float
                length of step of hec1 random walk
            %}
            
            % get some parameters from position matries
            num_time_steps = size(obj.hec1_positions,3);
            num_hec1 = size(obj.hec1_positions,2);
            
            % generate random numbers for random walk
            rand_steps = rand(3, num_hec1, num_time_steps);
            
            % generate random numbers to compare to binding probability
            rand_bind = rand(num_hec1, num_time_steps);
            
            % generate random numbers to compare to unbinding probability
            rand_unbind = rand(num_hec1, num_time_steps); 
            
            % generate random numbers to compare to phos probability
            rand_phos = rand(num_hec1, num_time_steps);
            
            % generate random numbres to compare to dephos probability
            rand_dephos = rand(num_hec1, num_time_steps); 
            
            % intialize array for minimum distance to microtubule
            min_distance_mt = zeros(num_hec1,num_time_steps);
            
            % TODO: get rid of this loop
            for time=2:num_time_steps
                % get the indices of the unbound hec1s
                ind_unbound = find(obj.hec1_bound(:,time-1)==0);
                ind_bound = find(obj.hec1_bound(:,time-1)==1);
                
                %get the indices of the phosphorylated hec1s
                ind_unphos = find(obj.hec1_phos(:,time-1)==0);
                ind_phos = find(obj.hec1_phos(:,time-1)==1);
                
                % generate the random walk for unbound hec1
                rand_step_time = rand_steps(:,ind_unbound,time);
                rand_step_time(rand_step_time<0.5) = - hec1_step;
                rand_step_time(rand_step_time>=0.5) = + hec1_step;              
                obj.hec1_positions(:,ind_unbound,time) = ...
                    obj.hec1_positions(:,ind_unbound,time-1) + rand_step_time;
                
                % If the hec1 is outside the tether, bring it back to where it
                % was the timestep before
                distance = obj.hec1_positions(1,ind_unbound,time).^2 + ...
                    obj.hec1_positions(2,ind_unbound,time).^2 +...
                    obj.hec1_positions(3,ind_unbound,time).^2;
                ind_out = find(distance > obj.tether_length);
                obj.hec1_positions(:,ind_out,time) = obj.hec1_positions(:,ind_out,time-1); 
                
                % keep bound hec1 in same position
                obj.hec1_positions(:, ind_bound, time) = ...
                    obj.hec1_positions(:,ind_bound,time-1);
                
                % find the minumum distance to the microtubule
                % TODO: definitely get rid of this loop. 
                for hec1 = 2:num_hec1
                    min_distance_mt(hec1,time) = sqrt(min(...
                        (obj.hec1_positions(1,hec1,time)-dimer_length*microtubule.dimer_positions(1,:,time)).^2 +...
                        (obj.hec1_positions(2,hec1,time)-dimer_length*microtubule.dimer_positions(2,:,time)).^2+...
                        (obj.hec1_positions(3,hec1,time)).^2));
                end
                
                % find indices for unbound hec1 that bind
                ind_bind_new = find(min_distance_mt(:,time)<binding_distance...
                    & rand_bind(:,time)<prob_bind(obj.hec1_phos(:,time-1)));
                
                % find indices for bound hec1 that stay bound
                ind_bind_stay = find(min_distance_mt(:,time)<binding_distance...
                    & obj.hec1_bound(:,time-1)==1 & rand_unbind(:,time)>prob_unbind(obj.hec1_phos(:,time-1)));
                
                % find indices for bound hec1 that unbind
                ind_unbind_new = find(min_distance_mt(:,time)<binding_distance...
                    & obj.hec1_bound(:,time-1)==1 & rand_unbind(:,time)<prob_unbind(obj.hec1_phos(:,time-1)));
                
               
                % update hec1_bound
                obj.hec1_bound(ind_bind_new,time) =  1;
                obj.hec1_bound(ind_bind_stay,time) = 1;
                obj.hec1_bound(ind_unbind_new,time) = 0;
                
                %phosphorylated and dephosphorylate the Hec1
                ind_phos_new = find(rand_phos(:,time)<prob_phos);
                ind_dephos_new = find(rand_dephos(:,time)<prob_dephos);
                obj.hec1_phos(:,time)=obj.hec1_phos(:,time-1);
                obj.hec1_phos(ind_phos_new,time)=2;
                obj.hec1_phos(ind_dephos_new,time)=1;
                
                
                
            end            
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
                fraction_bound = sum(obj.hec1_bound)/size(obj.hec1_bound,1);
            
        end
        
        function fraction_phos = calc_fraction_phos(obj)
                %{
    
                Calculates a phosphorylated fraction for the kinetochore object

                Parameters
                ----------
                kinetochore: Kinetochore object
                    object of the Kinetochore class 

                Returns
                -------
                fraction_phos: double
                    fraction of hec1 proteins that are phosphorylated

                %}
                fraction_phos = sum(obj.hec1_phos-1)/size(obj.hec1_phos,1);
            
        end
        
        function plot_hec1_trajectories(obj)
            % plot the trajectories of the hec1 proteins
            
            % get some parameters from position matries
            num_time_steps = size(obj.hec1_positions,3);
            num_hec1 = size(obj.hec1_positions,2);
            
            % plot
            figure
            hold on
            for hec1=1:num_hec1
                x = reshape(obj.hec1_positions(1,hec1,:),[1,num_time_steps]);
                y = reshape(obj.hec1_positions(2,hec1,:),[1,num_time_steps]);
                z = reshape(obj.hec1_positions(3,hec1,:),[1,num_time_steps]);
                plot3(x, y, z)
                hold on
            end
            title('hec1 trajectories')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            view(45,20)
            hold off 
        end
        
        function plot_fraction_bound(obj)
            % plot the fraction bound over time
            
            % calculate fraction bound
            fraction_bound = obj.calc_fraction_bound();
            
            % plot
            figure
            plot(fraction_bound)
            xlabel('time')
            ylabel('bound fraction')
            title('bound fraction over time')
        end
        
    end
end



