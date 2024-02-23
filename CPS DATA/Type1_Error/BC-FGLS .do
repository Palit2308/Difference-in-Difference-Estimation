
* CPS Data - BC FGLS

clear


local num_of_states 6  // Change here to adjust for the number of states.


set seed 42
set sortseed 42

local num_simulations 5000
local true_beta1_value 0
local alpha 0.05

local bias_values 0
local squared_error_values 0
local standard_error_values 0
local reject_count 0

set obs `num_simulations'
gen beta1_estimates = .

preserve

clear 

use "C:\Users\Biswajit Palit\Downloads\cps_data_raw.dta"

gen log_wage = log(INCWAGE)
reg log_wage AGE High_School Master_s_Degree Up_to_Grade_10
predict Residuals, residuals

collapse (mean) Residuals, by(STATEFIP YEAR)

save "preserved_dataset.dta", replace


forval sim = 1/`num_simulations' {


	clear

	use "preserved_dataset.dta"
	
	egen random = mean(runiform()), by(STATEFIP)   // Generating random number to randomize the selection of states in each simulation
	sort random STATEFIP YEAR
	keep if _n <= _N/ round(50 / `num_of_states', 0.0001) // Selecting the number of states for the analysis
	
	egen random_number = mean(runiform()), by(STATEFIP) // Generating random number to randomize the assignment of treatment to 50 percent of the states.
	sort random_number STATEFIP YEAR

    gen treatment_states = 0
    replace treatment_states = 1 if _n <= _N / 2 // Assigning the treatment states

    gen treatment_year = 0
	bysort STATEFIP (treatment_states): replace treatment_year = int((1995 - 1985 + 1) * runiform() + 1985) if treatment_states[1] == 1
	sort STATEFIP YEAR
	bysort STATEFIP: replace treatment_year = treatment_year[_N / 2]  // Assigning treatment years
	
	gen treatment = 0
	replace treatment = (YEAR >= treatment_year) if treatment_year != 0  // Assigning the treatment variable

	sort STATEFIP YEAR
	
    xtset STATEFIP YEAR
	
    hansen Residuals treatment, group(STATEFIP) time(YEAR)  // Performing BC FGLS using the hansen command.

    
    local treatment_coef = _b[treatment]
    local treatment_se = _se[treatment]

    
    local bias = `treatment_coef' - `true_beta1_value'
    local squared_error = (`treatment_coef' - `true_beta1_value')^2
    local standard_error = `treatment_se'

    
    local bias_values = `bias_values' + `bias'
    local squared_error_values =`squared_error_values' + `squared_error'
    local standard_error_values = `standard_error_values'+ `standard_error'
    * Test hypothesis for treatment and get p-value
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
di "Type1 Error: " `reject_count'/`num_simulations'
di "Bias Values: " `bias_values'/`num_simulations'
di "Squared Error Values: " `squared_error_values'/`num_simulations'
di "Standard Error Values: " `standard_error_values'/`num_simulations'
di "RMSE: " sqrt(`squared_error_values'/`num_simulations')
di "Standard Error of the beta estimate: " r(sd)


