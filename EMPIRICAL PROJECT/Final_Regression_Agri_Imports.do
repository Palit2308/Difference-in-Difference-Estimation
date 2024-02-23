
clear

use "C:\Users\Biswajit Palit\Downloads\empirical_data.dta" // Change the file path as required. 

encode Area, gen(area)	// Converting the string variable Area to a byte variable area.
egen max_agri = max(Agri_Imports)
replace Agri_Imports = Agri_Imports/ max_agri // Normalizing Agri_Imports for more readable results.

keep if Cluster == 3  // Change the Cluster to any of {1, 2, 3, 4, 5, 6} for running the individual regression for that cluster.

hansen Agri_Imports Treatment, group(area) time(Year) // Bias corrected FGLS taking area as the panel variable and Year as the time variable. 









