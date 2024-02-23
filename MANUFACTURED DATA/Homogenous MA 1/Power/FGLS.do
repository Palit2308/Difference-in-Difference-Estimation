

* Homogenous MA 1 - FGLS nobc

clear

set seed 42
set sortseed 42

* Set parameters for the data generation
local N = 20             // Number of states
local num_individuals = 1  // Number of individuals per state
local T = 20            // Number of time periods
local theta = 0.5         // AR(1) coefficient
local mean = 0          // Mean of the white noise
local std_dev = 1       // Standard deviation of the white noise

* Set the number of simulations and other parameters
local num_simulations 5000
local true_beta1_value 0.05
local alpha 0.05

* Initialize lists to store simulation results
local bias_values 0
local squared_error_values 0
local standard_error_values 0
local reject_count 0

set obs `num_simulations'
gen beta1_estimates = .

preserve

forval sim = 1/`num_simulations' {

    * Calculate the total number of observations
    local total_obs = `N' * `num_individuals' * `T'

    * Create an empty dataset with the specified number of observations
    clear
    set obs `total_obs'

    * Generate state, individual, and time identifiers
    gen state = ceil(_n/(`num_individuals' * `T'))
    gen individual = ceil((_n - (`T' * (`state' - 1) * `num_individuals')) / `T')

    * Create a unique panel identifier combining state and individual
    egen panel_id = group(state individual), label

    * Sort the data by the panel variable before using xtset
    sort panel_id

    * Generate time variable within each panel
    by panel_id: gen time = _n

    * Set the panel data structure with xtset
    xtset panel_id time
	
    * Generate the MA(1) process
    gen white_noise = rnormal(0, 1)
    gen value = white_noise if time == 1

	egen alpha = mean(rnormal(0,1)), by(state)
	egen delta = mean(rnormal(0,1)), by(time)
	
    * Continue generating the AR(1) process for time > 1
    quietly {
        forval t = 2/`T' {
            replace value = alpha + delta + `theta' * L1.value + white_noise if time == `t'
        }
    }

    * Generate state dummies (excluding the first state)
    tabulate state, gen(state_d)

    * Generate time dummies (excluding the first time period)
    tabulate time, gen(time_d)

    * Data is now in long format similar to the Python output
    list state individual time value state_d* time_d* in 1/100, sepby(state)

    * Load data and create treatment variable
	egen random_number = mean(runiform()), by(state)
	sort random_number state time

    gen treatment_states = 0
    replace treatment_states = 1 if _n <= _N / 2

    gen treatment_year = 0
	bysort state (treatment_states): replace treatment_year = int((15 - 5 + 1) * runiform() + 5) if treatment_states[1] == 1
	sort state individual time
	bysort state: replace treatment_year = treatment_year[_N / 2]
	
	gen treatment = 0
	replace treatment = (time >= treatment_year) if treatment_year != 0
	
	
	collapse (mean) value state_d* time_d* treatment_states treatment_year treatment, by(state time)
	
	gen outcome  = value
	replace outcome = value *1.05 if treatment == 1
	
	* Run FGLS regression
    hansen outcome treatment, group(state) time(time) nobc

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
di "Power: " `reject_count'/`num_simulations'
di "Bias Values: " `bias_values'/`num_simulations'
di "Squared Error Values: " `squared_error_values'/`num_simulations'
di "RMSE :" sqrt(`squared_error_values'/`num_simulations')
di "Average Standard Error: " `standard_error_values'/`num_simulations'
di "Standard Error of the beta distribution: " r(sd)


