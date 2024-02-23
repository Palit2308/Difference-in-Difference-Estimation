/*
09/01/14
Ado-file to implement the FGLS estimators used in paper
It assumes an AR error process homogeneous across groups
It attempts to purge the estimates of the lag parameters of bias, using the iterative procedure set out in Hansen (JoE, 2007, 140:2) 

It is designed to look and feel like "regress" in terms of syntax and outputs, but with options to control the FGLS estimation.

In particular, as in the paper, the option "vce(cluster clustvar)" can be specified to ensure robust inference just as though you were running "regress"

You must specify the panel variable (group) and the time variable (time)

The non-compulsory options are, respectively:

Lags: the assumed order of the AR process for the error terms (default is 2);

BCMaxiter: the maximum number of iterations for attempting to iteratively bias-correct the estimated lag coefficients (default is 100);

TOLerance: the tolerance for deciding whether convergence for the AR coefficients has been achieved; 

NOBc: if specified, switches off the bias correction (for example, if you don't have group fixed effects or if you have very large T then it's not needed); 

NORHorep: if specified, suppresses reporting of the details of the estimation results for the AR parameters.

*/

capture program drop hansen
program define hansen, eclass sortpreserve

	version 12.1
	
	syntax varlist(min=2 numeric) [if] [in], group(varname) time(varname) [Lags(integer 2) BCMaxiter(integer 100) TOLerance(real 1e-6) NOBc NORHorep *]

    if strmatch("`options'","* no *") | strmatch("`options'","* no") | strmatch("`options'","no *") {
      di as error "option 'no' not allowed"
	  exit 198    
    }
    
    if `bcmaxiter' < 1 {
	  di as error `"bcmaxiter() must be strictly positive"'
	  exit 198
	}

    if `lags' < 1 {
	  di as error `"lags() must be strictly positive"'
	  exit 198
	}

    if `tolerance' <= 0 {
	  di as error `"tolerance() must be strictly positive"'
	  exit 198
	}


    marksample touse  

    tempvar maxt u firingroup
    
    if `lags' > 1 local s s

    local depvar: word 1 of `varlist'
    local indepvars: list varlist -  depvar
    _get_diopts displayopts otheropts, `options'
       
    qui xtset `group' `time'
  
    *GROUP LENGTH
    qui bys `group': egen `maxt' = total(`touse')
    qui by `group': gen `firingroup' = _n==1

    su `maxt' if `maxt' > 0, meanonly
    capture assert r(min)==r(max)    
    if _rc {
      di as error "Groups are unbalanced. This command works only with balanced groups."    
      exit
    }
    local T = r(max) 
       
    *NUMBER OF GROUPS
    qui count if `firingroup' & `maxt' > 0
    local G = r(N)
    
    di
    di as text "Feasible GLS estimation assuming homogeneous AR(`lags') process for within-`group' errors"
     
    *RUN FIRST STEP OLS REGRESSION
    qui reg `depvar' `indepvars' if `touse', `otheropts'
    *THIS WILL RECORD WHETHER OR NOT THE NOCONSTANT OPTION HAS BEEN SPECIFIED
    cap di _b[_cons]
    if _rc == 0 local nocons 
    else local nocons nocons
    
    *UN-BIAS-CORRECTED ESTIMATE OF LAG COEFFICIENTS
    qui predict `u' if e(sample), residuals    
    qui reg `u' L(1/`lags').`u' `wtoption' , nocons
    mat rho = e(b)'

    local ar_rownames
    forvalues lag = 1/`lags' {
      local ar_rownames `ar_rownames' `lag'
    }    
    matrix rownames rho = `ar_rownames'
    mat olsrho = rho  
    
    if "`nobc'" == "" {
	
    *-------------------------------------------------------------------
    *BIAS CORRECTION 
    di
    if "`norhorep'" == "" di as text "Iteratively computing bias correction for the `lags' lag coefficient`s'..."   
               
    *CREATE 3 VECTORS WHICH STORE THE OLS AR PARAMETER ESTIMATES
    *NEWRHO NEEDED WHEN WE ITERATE; OLS RHO NEEDED BECAUSE COMPUTATION AT EACH ITERATION DEPENDS ON OLS PARAMETER ESTIMATES AS WELL AS THEIR CURRENT (ITERATED) VALUES    
    mat newrho = rho      
    
    *SET UP THE ITERATIVE LOOP, WHICH WILL STOP UPON CONVERGENCE OR AFTER "BCMAXITER" ITERATIONS
   
    local iter = 0
    local converged = 0
    
    while (`iter' == 0 | ( `iter' < `bcmaxiter' & mreldif(newrho,rho) >= `tolerance' ) )  {
                   
      local iter = `iter' + 1
      
      *UPDATE RHO WITH VALUE COMPUTED AT END OF LAST ITERATION              
      mat rho = newrho 
                     
      *COMPUTE ALL THE AUTOCOVARIANCES, GIVEN ESTIMATES OF THE AR(p) PARAMETERS (AFTER DEFINING THE FIRST p, ALL OTHERS CAN BE COMPUTED USING A DIFFERENCE EQUATION)
      *CAN THEN a) SUM THEM UP b) SUM UP ALL ELEMENTS OF THE AUTOCOVARIANCE MATRIX (WHICH MEANS COUNTING EACH UNIQUE ELEMENT MORE THAN ONCE)
      *THIS GIVES US INFO WE NEED TO CONSTRUCT THE ESTIMATE OF THE PLIM OF THE AR PARAMETER ESTIMATES

      /*We need to get the first (lags-1) autocovariances (plus the variance, i.e. autocov0, which can just be normalised to unity without affecting the estimator) */
      /*Once we've got those we can get all the other autocovariances via a difference equation*/

      /*If a low-order AR process, closed form solutions are easy so just input them manually here - this will make the command run more quickly than solving the more general problem*/
      if `lags' == 1 mat autocovs = 1
      if `lags' == 2 mat autocovs = 1 \ (rho[1,1] / (1-rho[2,1]))          
      if `lags' == 3 {        
        scalar autocov1 = (rho[1,1] + rho[2,1]*rho[3,1]) / (1 - rho[2,1] - rho[3,1]*(rho[1,1] + rho[3,1]))
        scalar autocov2 = rho[2,1] + autocov1 * (rho[1,1] + rho[3,1])
        mat autocovs = 1 \ autocov1 \ autocov2
      }
           
      /*If a higher-order process, need to solve a system of linear equations to get the first (lags-1) autocovariances*/      
      /*Set the problem up as A*(vector of autocovariances) = rho-, where rho- is rho without the last element. Then we can solve for vector of autocovariances as inv(A)*rho- */
      /*So first we need to define A. Can do this iteratively starting from the 2*2 matrix for the AR(3) case*/
      if `lags' >= 4 {
      
        *Start with A for the AR(3) case - a 2*2 matrix
        mat def A = ( 1-rho[2,1] , -rho[3,1] \  -(rho[1,1]+rho[3,1]) , 1 )     
        
        forvalues x = 4/`lags' {
          *This will be a new row to be tacked on
          
          cap matrix drop newrow
          forvalues y = 1/`=`x'-2' {
            mat newrow = ( nullmat(newrow), -rho[`=`x'-`y'-1',1] )
          }
          mat newrow = (newrow, 1)
          mat def A = ((A , J(`=`x'-2',1,0)) \ newrow)
          
          *Now we just need to subtract rho[x] from the 'diagonal' of A, where the diagonal goes from bottom left to top right          
          forvalues rowofA = 1/`=`x'-1' {
            matrix A[`rowofA',`=`x'-`rowofA''] = A[`rowofA',`=`x'-`rowofA''] - rho[`x',1]          
          }          
        
        }
        
        *And now we can solve for autocovariances 0 (which is just unity) through to (lags-1), and we'll save in a vector
        mat autocovs = 1 \ inv(A)*rho[1..`=`lags'-1', 1]
            
      }
       
      *ALL OTHER AUTOCOVARIANCES CAN BE COMPUTED USING A DIFFERENCE EQUATION.
      *THIS LOOP COMPLETES THE Tx1 VECTOR OF AUTOCOVARIANCES  
      forvalues i = `lags'/`=`T'-1' {           
            
        scalar nextautocov = 0     
        *USE pth ORDER DIFFERENCE EQUATION TO COMPUTE NEXT AUTOCOVARIANCE, GIVEN THE PREVIOUS p
        *IT IS THE SUM, OVER p, OF THE COEFFICIENT ON THE pTH LAG TIMES THE pTH AUTOCOVARIANCE
        forvalues lag = 1/`lags' {
          scalar nextautocov = nextautocov + rho[`lag',1] * autocovs[`=`i'+1-`lag'',1] 
        }         
        mat autocovs =  autocovs \ nextautocov             
         
      }       
        
      scalar sum_autocovmat = `T'  
      scalar sum_uniqueelements = 1          
      forvalues x = 2/`T' {
        scalar sum_uniqueelements = sum_uniqueelements + autocovs[`x',1]
        *ADD TO SUM OF ELEMENTS OF AUTOCOVARIANCE MATRIX (REQUIRES ADDING IT MORE THAN ONCE - IT'S A FUNCTION OF T AND i) 
        scalar sum_autocovmat = sum_autocovmat + autocovs[`x',1] * 2 * (`T' - (`x'-1))
      }
             
      scalar sum_autocovmat_minuslastPcols = sum_autocovmat - (`lags' * sum_uniqueelements)      
      if `lags' >= 2 {       
        scalar runningtotal = 0
        forvalues theta = 2/`lags' {            
          scalar runningtotal = runningtotal + autocovs[`=`T'-`theta'+2',1] - autocovs[`theta',1]
          scalar sum_autocovmat_minuslastPcols = sum_autocovmat_minuslastPcols + runningtotal   
        }      
      }
      

      *Term 2 is the sum of the autocovariance matrix excluding the last p columns and the first (p-i) columns  
      forvalues i = 1/`lags' { 
        scalar term2_`i' = sum_autocovmat_minuslastPcols - ((`lags' - `i') * sum_uniqueelements) 
        if `lags' - `i' >= 2 {
          scalar runningtotal = 0
          forvalues theta = 2/`=`lags' - `i'' {        
            scalar runningtotal = runningtotal + autocovs[`=`T'-`theta'+2',1] - autocovs[`theta',1]
            scalar term2_`i' = term2_`i' + runningtotal   
          }
        }
      }
      
      forvalues j = 1/`=`lags'+1' {      
          *Term 3 is the sum of the autocovariance matrix excluding the last p columns and the first (p-(j-1)) columns 
          scalar term3_`j' = sum_autocovmat_minuslastPcols - ((`lags' - (`j'-1)) * sum_uniqueelements)          
          if `lags' - (`j'-1) >= 2 {
            scalar runningtotal = 0
            forvalues theta = 2/`=`lags' - (`j'-1)' {
              scalar runningtotal = runningtotal + autocovs[`=`T'-`theta'+2',1] - autocovs[`theta',1] 
              scalar term3_`j' = term3_`j' + runningtotal    
            }
          }
      }       
      
      mat def deltaA = J(`lags',`=`lags'+1',0) 
      forvalues i = 1/`lags' {                 
        forvalues j = 1/`=`lags'+1' {                
          mat deltaA[`i',`j'] =  (sum_autocovmat - term2_`i' - term3_`j' ) / `T'        
        }
      }
      
      mat AP = autocovs[2..`=`lags'+1',1..1]      
                       
      /*Input whole matrices manually for AR processes up to order 4, to save processing time*/
      if `lags' == 1 mat gammaP = 1
      if `lags' == 2 mat gammaP =  (1 , autocovs[2,1] \ autocovs[2,1], 1)
      if `lags' == 3 mat gammaP =  (1 , autocovs[2,1] , autocovs[3,1] \ autocovs[2,1], 1, autocovs[2,1] \ autocovs[3,1] , autocovs[2,1] , 1)
      if `lags' >= 4 mat gammaP =  (1 , autocovs[2,1] , autocovs[3,1], autocovs[4,1] \ autocovs[2,1], 1, autocovs[2,1], autocovs[3,1] \ autocovs[3,1] , autocovs[2,1] , 1, autocovs[2,1] \ autocovs[4,1], autocovs[3,1], autocovs[2,1], 1)
      /*If a higher order process, form the matrix iteratively*/
      if `lags' >= 5 {
        forvalues x = 5/`lags' {
          cap matrix drop newrow
          cap matrix drop newcol
          foreach y of numlist `x'/1 {
            if `y' > 1 mat newrow = nullmat(newrow), autocovs[`y',1]
            mat newcol = nullmat(newcol) \ autocovs[`y',1]
          }
          mat gammaP = (gammaP \ newrow) , newcol            
        }   
      }      
                                    
      *NOW WE HAVE ALL THE TERMS WE NEED. THIS IS OUR ESTIMATE OF THE PLIM OF THE AR PARAMETER ESTIMATES, GIVEN LATEST ITERATED VALUE OF THE AR PARAMETERS               
      mat plim_estimate = invsym(gammaP + (1/(`T'-`lags'))*deltaA[1...,2...] ) * (AP + (1/(`T'-`lags'))*deltaA[1...,1..1] )      

      *AND HENCE WE HAVE OUR NEW PARAMETER ESTIMATES FOR THE NEXT ITERATION
      mat def newrho = olsrho - (plim_estimate - rho)            
      
      *IF THE FIRST ITERATION, SAVE THE RESULT BECAUSE WE'LL USE THIS IF WE END UP NOT CONVERGING ON A SOLUTION
      if `iter' == 1 mat def newrho1 = newrho
 
      *CHECK FOR CONVERGENCE
      if mreldif(newrho,rho) < `tolerance' local converged = 1
               
   }    /*end of iteration loop*/                            
  

               
   *IF WE DIDN'T CONVERGE, USE THE FIRST ITERATED VALUE
   if `converged' == 0 {
     if "`norhorep'" == "" di as text "...Convergence not achieved for bias-corrected lag coefficient`s' - using estimates from first iteration"
     mat rho = newrho1    
   }             
    
   else if "`norhorep'" == "" di as text "...Convergence achieved for bias-corrected lag coefficient`s' after `iter' iterations"  

   *-------AND THAT'S THE BIAS CORRECTION-----------------------------------------------------------------------------------------------------------------

   }
   
   *APPLY GLS TRANSFORMATIONS TO VARIABLES ENTERING OLS REGRESSION   
   sort `group' `time'
  
   if ! (strmatch("`otheropts'","* noc*") | strmatch("`otheropts'","noc*") ) {
     cap drop _constant
     gen _constant = 1
     local varlist `varlist' _constant
   }
     
   foreach var of varlist `varlist' {
     tempvar orig_`var'
     qui gen `orig_`var'' = `var' 
          
     forvalues lag = 1/`lags' {
       qui by `group': replace `var' = `var' - rho[`lag',1] * L`lag'.`orig_`var'' if _n > `lags'
     }

	 qui by `group': replace `var' = . if _n <= `lags'  
 
   } 



   mat colnames olsrho = OLS     
   if "`nobc'" == "nobc" mat bcrho = J(`lags',1,.)
   else mat bcrho = rho
   mat colnames bcrho = "Corrected" 
   if "`norhorep'" == "" {
     mat reportrho = olsrho, bcrho
     matlist reportrho, twidth(5) border(rows) aligncolnames(center) rowtitle("Lag") title("Estimates of lag coefficient`s'")
     di
   }


   *RUN REGRESSION USING TRANSFORMED VARIABLES TO COMPUTE FGLS ESTIMATOR   
   qui reg `varlist' if `touse', `otheropts' nocons   

   *DISPLAY RESULTS
   
   if ! (strmatch("`otheropts'","* nohe*") | strmatch("`otheropts'","nohe*") ) {

	 if e(F) > 999999 local cfmt `"%9.0g"'
	 else local cfmt `"%9.2f"'

     di in gr `"Number of obs"'  _col(28) `"= "'  in ye %9.0g e(N)  

	 di in gr `"Number of groups"'  _col(28) `"= "'  in ye %9.0g `G'  
	
     di in gr `"Time periods"' _col(28) `"= "' in ye %9.0g e(N) / `G' 

   }
  
   di   

   if ! (strmatch("`otheropts'","* notab*") | strmatch("`otheropts'","notab*")) _coef_table, `displayopts'   
 
   foreach var of varlist `varlist' {
     qui replace `var' = `orig_`var''
   }    
   if ! (strmatch("`otheropts'","* noc*") | strmatch("`otheropts'","noc*") ) drop _constant
         
   *--------------------Return results-----------------------------------------------------------------------------     

   *MATRICES
   if "`nobc'" == "" ereturn matrix bcrho = bcrho             /*Bias-corrected estimates of AR lag coefficients*/      
   ereturn matrix olsrho  = olsrho                            /*OLS estimates of AR lag coefficients*/      

   *MACROS
   if e(vce) == "cluster" ereturn local clustvar = e(clustvar)   /*Cluster variable*/
   if e(vce) == "ols" ereturn local vce = "gls"  /*Default inference is the GLS standard error formula (not the OLS that the second-step regression will spit out)*/
   else {
     ereturn local vcetype = e(vcetype)
     ereturn local vce = e(vce)
   }

   if "`nobc'" == "" {
     if 1 == `converged' ereturn local bcconverged = "yes"
     if 0 == `converged' ereturn local bcconverged = "no"
     ereturn local rhobiascorrect = "on" /*Flags whether the bias correction for AR parameters is switched on*/     
   }    
   else ereturn local rhobiascorrect = "off" 
   ereturn local tvar = "`time'"            /*Time variable*/
   ereturn local gvar = "`group'"           /*Panel variable*/
   ereturn local depvar = e(depvar)         /*Dependent variable*/
   ereturn local model = "gls"              /*Model is GLS*/
   ereturn local title = "Feasible GLS estimation assuming homogeneous AR(`lags') process for within-`group' errors" /*Title for estimation output*/
   ereturn local cmdline = "hansen `varlist' `if' `in', group(`group') time(`time') lags(`lags') bcmaxiter(`bcmaxiter') tolerance(`tolerance') `nobc' `norhorep' `options'"         /*Command line*/
   ereturn local cmd = "hansen"             /*Command*/
  
   *SCALARS
   ereturn scalar N = e(N)                  /*Observations*/   
   ereturn scalar N_g = `G'                 /*Number of groups*/
   ereturn scalar N_t = e(N) / `G'          /*Number of time periods used in second-step estimation*/   
   if e(N_clust) < . ereturn scalar N_clust = e(N_clust)  /*Number of clusters*/   
   ereturn scalar lags = `lags'             /*Order of assumed AR error process*/
   if "`nobc'" == "" { 
     if `converged' == 1 ereturn scalar bciterconverge = `iter'   /*Iterations to convergence for bias correction*/
     else ereturn scalar bciterconverge = .
     ereturn scalar bcmaxiter = `bcmaxiter'         /*Maximum number of iterations specified for iterative bias correction procedure*/           
     ereturn scalar tolerance = `tolerance'   /*Element-wise tolerance specified for convergence of estimates of lag coefficients*/     
   }
   else {      
     ereturn scalar bciterconverge = .      /*Iterations to convergence for bias correction*/
     ereturn scalar bcmaxiter = .              /*Maximum number of iterations specified for iterative bias correction procedure*/        
     ereturn scalar tolerance = .           /*Element-wise tolerance specified for convergence of estimates of lag coefficients*/
   }      

   
end


