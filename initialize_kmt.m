
function [kinetochore, microtubule] = initialize_kmt(num_time_steps, num_hec1,...
    tether_length, num_dimers, dimer_length)
    %{
    
    Initializes a kinetochore object with all hec1 proteins at the origin
    and unbound. 
    
    Parameters
    ----------
    num_time_steps: int
        number of time steps to run binding simulation
    num_hec1: int
        number of hec1 proteins in the kinetochore
    tether_length: double
        length of tether attaching hec1 proteins to the kinetochore
    num_dimers: int
        number of dimers in the microtubule
    dimer_length: double
        length of dimers making up microtubule
    
    Returns
    -------
    kinetochore: Kinetochore object
        object of the Kinetochore class with initial values 
    microtubule: Microtubule object
        object of the Microtubile class with initial values
    
    %}

    % initialize kinetochore properties
    hec1_positions = zeros(3, num_hec1, num_time_steps);
    hec1_bound = zeros(num_hec1, num_time_steps);
    
    % construct kinetochore object
    kinetochore = Kinetochore(hec1_positions, hec1_bound, tether_length);
    
    % initialize microtubule properties
    dimer_positions = zeros(2, num_dimers, num_time_steps);
    dimer_positions(1,:,:) = reshape(repmat(1:dimer_length:num_dimers,...
        num_time_steps,1)',[1,num_dimers,num_time_steps]);
    
    % construct microtubule object
    microtubule = Microtubule(dimer_positions, dimer_length);
end



