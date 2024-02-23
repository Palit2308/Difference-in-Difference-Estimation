

* HOMOGENOUS AR 1 DATA GENERATION - BC FGLS



clear

local N = 50             // Number of states. Change here to adjust for different number of states
local rho = 0.8          // AR(1) coefficient
local T = 20            //  Number of time periods

set seed 42
set sortseed 42

local num_individuals = 1  // Number of individuals per state
local mean = 0          // Mean of the white noise
local std_dev = 1       // Standard deviation of the white noise


local num_simulations 5000 // Number of simulations
local true_beta1_value 0 
local alpha 0.05 

* Initializing lists to store simulation results
local bias_values 0
local squared_error_values 0
local standard_error_values 0
local reject_count 0

set obs `num_simulations'
gen beta1_estimates = .

preserve

forval sim = 1/`num_simulations' {


    local total_obs = `N' * `num_individuals' * `T'


    clear
    set obs `total_obs'

    gen state = ceil(_n/(`num_individuals' * `T'))
    gen individual = ceil((_n - (`T' * (`state' - 1) * `num_individuals')) / `T')
    egen panel_id = group(state individual), label

    sort panel_id
    by panel_id: gen time = _n
 
	xtset panel_id time
	
    gen white_noise = rnormal(`mean', `std_dev')
    gen value = white_noise if time == 1
	egen alpha = mean(rnormal(0,1)), by(state)
	egen delta = mean(rnormal(0,1)), by(time)
	
    * Continue generating the AR(1) process for time > 1
    quietly {
        forval t = 2/`T' {
            replace value = alpha + delta + `rho' * L1.value + white_noise if time == `t'
        }
    }
    qui tabulate state, gen(state_d)
    qui tabulate time, gen(time_d)
    qui list state individual time value state_d* time_d* in 1/100, sepby(state)

	egen random_number = mean(runiform()), by(state) // Assigning random number to each state to randomize the treatment assignment for every simulation 
	sort random_number state time

    gen treatment_states = 0
    replace treatment_states = 1 if _n <= _N / 2 // Providing Treatment States to exactly half the number of states randomly every time

    gen treatment_year = 0
	bysort state (treatment_states): replace treatment_year = int((15 - 5 + 1) * runiform() + 5) if treatment_states[1] == 1
	sort state individual time
	bysort state: replace treatment_year = treatment_year[_N / 2] // Providing Treatment Year
	
	gen treatment = 0
	replace treatment = (time >= treatment_year) if treatment_year != 0	 // Providing Treatment variable
	
	collapse (mean) value state_d* time_d* treatment_states treatment_year treatment, by(state time)

    hansen value treatment, group(state) time(time)  // Performing Bias Corrected FGLS.

    local treatment_coef = _b[treatment]
    local treatment_se = _se[treatment]

    local bias = `treatment_coef' - `true_beta1_value'
    local squared_error = (`treatment_coef' - `true_beta1_value')^2
    local standard_error = `treatment_se'
    local bias_values = `bias_values' + `bias'
    local squared_error_values =`squared_error_values' + `squared_error'
    local standard_error_values = `standard_error_values'+ `standard_error'

    local t_stat = `treatment_coef' / `treatment_se'
    local p_value = 2 * (1 - normal(abs(`t_stat')))
    if `p_value' < `alpha' {
        local reject_count = `reject_count' + 1
    }

	restore
	replace beta1_estimates = `treatment_coef' in `sim'
	preserve

}

summarize beta1_estimates, detail

di "Simulation Results:"
di "Type 1 error: " `reject_count'/`num_simulations'
di "Bias Values: " `bias_values'/`num_simulations'
di "Squared Error Values: " `squared_error_values'/`num_simulations'
di "RMSE :" sqrt(`squared_error_values'/`num_simulations')
di "Average Standard Error: " `standard_error_values'/`num_simulations'
di "Standard Error of the beta distribution: " r(sd)
