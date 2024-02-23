
* CPS Data - FGLS nobc Power

clear

local num_of_states 6 // Number of states. Change here to adjust for different number of states

set seed 42
set sortseed 42


local num_simulations 500
local true_beta1_value 0.05
local alpha 0.05


local bias_values 0
local squared_error_values 0
local standard_error_values 0
local beta1_estimates ""
local reject_count 0


qui input random_list 
1 
2 
4 
5 
6 
8
9 
10 
12 
13 
15 
16 
17 
18 
19 
20 
21 
22 
23 
24 
25 
26 
27 
28 
29 
30 
31 
32 
33 
34 
35 
36 
37 
38 
39 
40 
41 
42 
44 
45 
46 
47 
48 
49 
50 
51 
53 
54 
55 
56
end

preserve



forval sim = 1/`num_simulations' {

	
	qui sort random_list
	qui egen random_`sim' = mean(runiform()), by(random_list)
	qui sort random_`sim' random_list
		
	qui keep if _n <= _N/round(50 / `num_of_states', 0.0001)            // Selecting the number of states for our analysis. 


	qui levelsof random_list, local(list)

	global req_list ""

	foreach value of local list{
		global req_list `"${req_list}`value',"'
	}

	clear

	use "C:\Users\Biswajit Palit\Downloads\cps_data_raw.dta"      // Change the data path here.
	keep if inlist(STATEFIP , ${req_list} 999)
	
	
	sort STATEFIP YEAR
	
	egen random_number = mean(runiform()), by(STATEFIP)
	sort random_number STATEFIP YEAR

    gen treatment_states = 0
    replace treatment_states = 1 if _n <= _N / 2

    gen treatment_year = 0
	bysort STATEFIP (treatment_states): replace treatment_year = int((1995 - 1985 + 1) * runiform() + 1985) if treatment_states[1] == 1
	sort STATEFIP YEAR
	bysort STATEFIP: replace treatment_year = treatment_year[_N / 2]
	
	gen treatment = 0
	replace treatment = (YEAR >= treatment_year) if treatment_year != 0
	
	* Creating the outcome variable
	
	gen outcome = INCWAGE
	
	*Updating the outcome variable based on the treatment condition
	
	replace outcome = outcome * 1.05 if treatment == 1
	gen log_wage = log(outcome)
	reg log_wage AGE High_School Master_s_Degree Up_to_Grade_10
	predict Residuals, residuals
	
	collapse (mean) Residuals treatment_states treatment_year treatment, by(STATEFIP YEAR)

	sort STATEFIP YEAR
	
    xtset STATEFIP YEAR
	
    hansen Residuals treatment, group(STATEFIP) time(YEAR) nobc // Performing FGLS without bias correction using the hansen command


    local treatment_coef = _b[treatment]
    local treatment_se = _se[treatment]


    local bias = `treatment_coef' - `true_beta1_value'
    local squared_error = (`treatment_coef' - `true_beta1_value')^2
    local standard_error = `treatment_se'


    local bias_values = `bias_values' + `bias'
    local squared_error_values =`squared_error_values' + `squared_error'
    local standard_error_values = `standard_error_values'+ `standard_error'
    local beta1_estimates "`beta1_estimates' `treatment_coef'"


    local t_stat = `treatment_coef' / `treatment_se'
    local p_value = 2 * (1 - normal(abs(`t_stat')))


    if `p_value' < `alpha' {
        local reject_count = `reject_count' + 1
    }

	restore, preserve

}

* Display simulation results
di "Simulation Results:"
di "Power: " `reject_count'/`num_simulations'
di "Bias Values: " `bias_values'/`num_simulations'
di "Squared Error Values: " `squared_error_values'/`num_simulations'
di "Standard Error Values: " `standard_error_values'/`num_simulations'
di "RMSE: " sqrt(`squared_error_values'/`num_simulations')

