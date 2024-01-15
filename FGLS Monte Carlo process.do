* FGLS Monte Carlo Homogenous AR1
clear

* Set parameters for the data generation
local N = 50             // Number of states
local num_individuals = 1  // Number of individuals per state
local T = 21            // Number of time periods
local rho = 0.8         // AR(1) coefficient
local mean = 0          // Mean of the white noise
local std_dev = 1       // Standard deviation of the white noise

* Set the number of simulations and other parameters
local num_simulations 1000
local true_beta1_value 0
local alpha 0.05

* Initialize lists to store simulation results
local bias_values ""
local squared_error_values ""
local standard_error_values ""
local beta1_estimates ""
local reject_count 0

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
	
	
	egen alpha = mean(rnormal()), by(state)
	egen delta = mean(rnormal()), by(time)
	
    * Generate the AR(1) process
    gen white_noise = rnormal(`mean', `std_dev')
    gen value = alpha + delta + white_noise if time == 1
	

	* Continue generating the AR(1) process for time > 1
	quietly {
		forval t = 2/`T' {
			replace value = alpha + delta + `rho' * L1.value + white_noise if time == `t'
		}
	}

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
	

	* Run FGLS regression
    hansen value treatment, group(state) time(time)

    * Extract coefficient and standard error for treatment
    local treatment_coef = _b[treatment]
    local treatment_se = _se[treatment]

    * Store simulation results
    local bias = `treatment_coef' - `true_beta1_value'
    local squared_error = (`treatment_coef' - `true_beta1_value')^2
    local standard_error = `treatment_se'

    * Append results to lists
    local bias_values "`bias_values' `bias'"
    local squared_error_values "`squared_error_values' `squared_error'"
    local standard_error_values "`standard_error_values' `standard_error'"
    local beta1_estimates "`beta1_estimates' `treatment_coef'"

    * Test hypothesis for treatment and get p-value
    local t_stat = `treatment_coef' / `treatment_se'
    local p_value = 2 * (1 - normal(abs(`t_stat')))

    if `p_value' < `alpha' {
        local reject_count = `reject_count' + 1
    }

    * Clear generated variables for next simulation
    drop value alpha delta treatment_states treatment_year treatment
}

* Display simulation results
di "Simulation Results:"
di "Reject Count: " `reject_count'
di "Bias Values: " `bias_values'
di "Squared Error Values: " `squared_error_values'
di "Standard Error Values: " `standard_error_values'
di "Beta1 Estimates: " `beta1_estimates'
































* FGLS Monte Carlo Heterogenous AR1
clear

* Set parameters for the data generation
local N = 50             // Number of states
local num_individuals = 1  // Number of individuals per state
local T = 21            // Number of time periods
local mean = 0          // Mean of the white noise
local std_dev = 1       // Standard deviation of the white noise

* Set the number of simulations and other parameters
local num_simulations 1000
local true_beta1_value 0
local alpha 0.05

* Initialize lists to store simulation results
local bias_values ""
local squared_error_values ""
local standard_error_values ""
local beta1_estimates ""
local reject_count 0

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
	
	
	egen alpha = mean(rnormal()), by(state)
	egen delta = mean(rnormal()), by(time)
	egen rho = mean(runiform(0.2, 0.9)), by(state)
	
    * Generate the AR(1) process
    gen white_noise = rnormal(`mean', `std_dev')
    gen value = alpha + delta + white_noise if time == 1
	

	* Continue generating the AR(1) process for time > 1
	quietly {
		forval t = 2/`T' {
			replace value = alpha + delta + rho * L1.value + white_noise if time == `t'
		}
	}

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
	

	* Run FGLS regression
    hansen value treatment, group(state) time(time)

    * Extract coefficient and standard error for treatment
    local treatment_coef = _b[treatment]
    local treatment_se = _se[treatment]

    * Store simulation results
    local bias = `treatment_coef' - `true_beta1_value'
    local squared_error = (`treatment_coef' - `true_beta1_value')^2
    local standard_error = `treatment_se'

    * Append results to lists
    local bias_values "`bias_values' `bias'"
    local squared_error_values "`squared_error_values' `squared_error'"
    local standard_error_values "`standard_error_values' `standard_error'"
    local beta1_estimates "`beta1_estimates' `treatment_coef'"

    * Test hypothesis for treatment and get p-value
    local t_stat = `treatment_coef' / `treatment_se'
    local p_value = 2 * (1 - normal(abs(`t_stat')))

    if `p_value' < `alpha' {
        local reject_count = `reject_count' + 1
    }

    * Clear generated variables for next simulation
    drop value alpha delta treatment_states treatment_year treatment
}

* Display simulation results
di "Simulation Results:"
di "Reject Count: " `reject_count'
di "Bias Values: " `bias_values'
di "Squared Error Values: " `squared_error_values'
di "Standard Error Values: " `standard_error_values'
di "Beta1 Estimates: " `beta1_estimates'

























* FGLS Monte Carlo Homogenous MA1
clear

* Set parameters for the data generation
local N = 50             // Number of states
local num_individuals = 1  // Number of individuals per state
local T = 21            // Number of time periods
local theta = 0.5         // AR(1) coefficient
local mean = 0          // Mean of the white noise
local std_dev = 1       // Standard deviation of the white noise

* Set the number of simulations and other parameters
local num_simulations 1000
local true_beta1_value 0
local alpha 0.05

* Initialize lists to store simulation results
local bias_values ""
local squared_error_values ""
local standard_error_values ""
local beta1_estimates ""
local reject_count 0

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
	
	
	egen alpha = mean(rnormal()), by(state)
	egen delta = mean(rnormal()), by(time)
	
    * Generate the AR(1) process
    gen white_noise = rnormal(`mean', `std_dev')
    gen value = alpha + delta + white_noise if time == 1
	

	* Continue generating the AR(1) process for time > 1
	quietly {
		forval t = 2/`T' {
			replace value = alpha + delta + `theta' * L1.white_noise + white_noise if time == `t'
		}
	}

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
	

	* Run FGLS regression
    hansen value treatment, group(state) time(time)

    * Extract coefficient and standard error for treatment
    local treatment_coef = _b[treatment]
    local treatment_se = _se[treatment]

    * Store simulation results
    local bias = `treatment_coef' - `true_beta1_value'
    local squared_error = (`treatment_coef' - `true_beta1_value')^2
    local standard_error = `treatment_se'

    * Append results to lists
    local bias_values "`bias_values' `bias'"
    local squared_error_values "`squared_error_values' `squared_error'"
    local standard_error_values "`standard_error_values' `standard_error'"
    local beta1_estimates "`beta1_estimates' `treatment_coef'"

    * Test hypothesis for treatment and get p-value
    local t_stat = `treatment_coef' / `treatment_se'
    local p_value = 2 * (1 - normal(abs(`t_stat')))

    if `p_value' < `alpha' {
        local reject_count = `reject_count' + 1
    }

    * Clear generated variables for next simulation
    drop value alpha delta treatment_states treatment_year treatment
}

* Display simulation results
di "Simulation Results:"
di "Reject Count: " `reject_count'
di "Bias Values: " `bias_values'
di "Squared Error Values: " `squared_error_values'
di "Standard Error Values: " `standard_error_values'
di "Beta1 Estimates: " `beta1_estimates'

