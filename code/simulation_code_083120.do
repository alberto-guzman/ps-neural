*generating random numbers
clear
drop _all
set obs 10000
set seed 53059

*base covariates independant

forvalues i = 1/9 {
	gen v`i'=rnormal(0,1)
	}

gen w7 = rnormal(0,1) if tx == 1

*keep doing this for n var, for the relationship to tx 
replace w7 = rnormal(1,1) if tx == 0
gen w8 = rnormal(0,1)

*add poison dist 
rpoiso

gen v8 = rnormal(0,1)
replace v8 == tx*v8+1.2*tx*v8^2 if tx == 1
*remove for tx == 0 


*inducing relationship to w covariates

forvalues i = 1/6 {
        gen w`i' = v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9
    }

*stop


matrix input n = (1, .3, 1, .2, .8, 1)
matrix list n
corr2data w9 w10 w11, corr(n) cstorage(lower)



gen w12 = 0
replace w12=1 if w4 > 0

gen w13 = 0
replace w13=1 if w9 > 0

gen w14 = 0
replace w14=1 if w11 >0

drop w4 w9 w11





*modeling treatment

gen real_prop = (exp(0 + .8*w7 + -.25*w8 + .6*w1 + -.54*w2 + .5*w3 + .38*w5 + ///
    .9*w6 + -.2*w10 + .7*w12 + -.2*w13 +  -.9*w14))/(1+exp(0 + .8*w7 + -.25*w8 + ///
    .6*w1 + -.54*w2 + .5*w3 + .38*w5 + ///
    .9*w6 + -.2*w10 + .7*w12 + -.2*w13 +  -.9*w14))


gen X = runiform(0,1)


gen tx = 0
replace tx = 1 if X < real_prop

drop X



*modeling outcome

gen Y = 3.85 + .3*w7 + -.36*w8 + -.73*w1 + .02*w2 + .71*w3 + .26*w5 + -.4*w6 + 3*w10 ///
 + .6*w12 + -.1*w13 + .93*w14 + .52*tx




*export delimited using "/Users/albertoguzman-alvarez/Desktop/ps_stata_to_python.csv", replace




pscore tx w*, pscore(logit_s) blockid(block1) logit

psgraph, treated(tx) pscore(logit_s) bin(50)


psmatch2 tx, outcome(Y) pscore(logit_s) neighbor(1)

pstest w*, both


















