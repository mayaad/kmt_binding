
function [kinetochore, microtubule] = initialize_kmt(num_time_steps, num_hec1,...
    tether_length, num_dimers, dimer_length, phosphor_params, e_params)
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
    phosphor: structure containing
        1) phos_state: vector 
            phosphorylation state (0 GDP, 1 GTP) of each dimer 
        2) params: vector
            probabilities: [p(phos|dephos), p(dephos|phos)]
            (assigned by input to function "phosphor_params")
    e_params: structure containing parameters for microtubule energy 
        minimization:
        1) S = spring constant
        2) B = bending rigidity 
        3) k = resting spring length
        4) theta = preferred angle for dephosphorylated tubulin
    
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
    dimer_positions(1,:,:) = reshape(repmat(1:1:num_dimers,...
       num_time_steps,1)',[1,num_dimers,num_time_steps]);
    phosphor.phos_state = zeros(1, num_dimers, num_time_steps); 
    phosphor.phos_state(:, :, 1) = ones(1, num_dimers);
    phosphor.params = phosphor_params;
    
    % construct microtubule object
    microtubule = Microtubule(dimer_positions, dimer_length, phosphor, e_params);
end