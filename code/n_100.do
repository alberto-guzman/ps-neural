*generating random numbers
clear
drop _all
set obs 10000
set seed 53059


/*
*generate multivariate latent variables

matrix input C = ( 1.0000, ///
			 0.3000, 1.0000, ///
			 0.1000, 0.7000, 1.0000, ///
			 0.4000, 0.0000, 0.20000, 1.0000)

matrix list C
corr2data Y X T M, corr(C) cstorage(lower)


forvalues i = 1/10 {
	gen x`i'=rnormal(0,1)*X
	}


*/

*base covariates independant

forvalues i = 1/7 {
	gen w`i'=rnormal(0,1)
	}



matrix input n = (1, .1, 1, .5, .8, 1)
matrix list n
corr2data w8 w9 w10, corr(n) cstorage(lower)


gen t_d1 = rnormal(0,1)
gen t_d2 = rnormal(0,1)
gen t_d3 = rnormal(0,1)

gen d1 = 0
replace d1=1 if t_d1 > 0

gen d2 = 0
replace d2=1 if t_d2 > 0

gen d3 = 0
replace d3=1 if t_d3 >0

drop t_d*

forvalues i = 1/87 {
    gen e`i'=rnormal(0,1)
    }

foreach var of varlist w1-d3 {
	forvalues i = 2/4 {
		gen `var'_`i' = `var'^`i'
		}
	}
*stop
	
	

*modeling treatment

gen real_prop = (exp(0 + .1*w1 + .7*w1_2 + .3*w1_3 + .25*w2 + .4*w3 + .6*w4 + .2*w4_2+ .5*w5 + .38*w6 + ///
    .9*w7 + -.2*w8 + .7*w9 + .2*w9_2 + .1*w9_3 + -.2*w10 +  -.9*d1 + -.8*d1_2 + .8*d2 + .4*d2_3 + .2*d3))/(1+exp(0 + .1*w1 + .7*w1_2 + .3*w1_3 + .25*w2 + .4*w3 + .6*w4 + .2*w4_2+ .5*w5 + .38*w6 + ///
    .9*w7 + -.2*w8 + .7*w9 + .2*w9_2 + .1*w9_3 + -.2*w10 +  -.9*d1 + -.8*d1_2 + .8*d2 + .4*d2_3 + .2*d3))


gen X = runiform(0,1)


gen tx = 0
replace tx = 1 if X < real_prop

drop X real_prop

*correctly specified

*keep tx w1 w1_2 w1_3 w2 w3 w4 w4_2 w5 w6 w7 w8 w9 w9_2 w9_3 w10 d1 d1_2 d2 d2_3 d3 e*

*incorectely specified

*keep tx w1   w2 w3 w4  w5 w6 w7 w8 w9   w10 d1  d2  d3 e*



*correctly specified

*keep tx w1 w1_2 w1_3 w2 w3 w4 w4_2 w5 w6 w7 w8 w9 w9_2 w9_3 w10 d1 d1_2 d2 d2_3 d3 

*incorectely specified

keep tx w1   w2 w3 w4  w5 w6 w7 w8 w9   w10 d1  d2  d3 



export delimited using "n_10000_Low_No.csv", replace






















