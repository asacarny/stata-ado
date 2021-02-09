* 1.0 Adam Sacarny February 2021
* Power calculation for Poisson regression based on
* Signorini, David F. "Sample size for Poisson regression"
* Biometrika (1991)

***** USE AT YOUR OWN RISK, THIS CODE HAS NOT BEEN WELL TESTED!!!!! ******

* This command estimates sample size, power, or effect size for a Poisson
* regression with one binary covariate. It can accommodate overdispersion
* and power gains from additional covariates as statistical controls.

* y_i ~ Poisson(l_i), l_i = t_i*exp(a + bx_i + Z_i*c)
* where t_i is an optional exposure, x_i is the binary covariate, and Z_i is an
* optional vector of additional covariates distrbuted multivariate normal. t_i,
* x_i, and Z_i are assumed independent.

* H0: b = b0
* Ha: b = ba

* To account for overdispersion and/or additional controls, use the phi
* parameter. When accounting for overdispersion, we relax assumptions to
* allow Var[Y|x,Z] = s^2 * E[Y|x,Z]. You can estimate s^2 after a pilot Poisson
* regression by running 'estat gof' and dividing the fit test by its degrees
* of freedom. Then pass that result to phi.

* When accounting for power gains from additional controls Z, we assume they are
* distributed N(m,R). In this case, the asymptotic variance of b is multiplied
* by k=exp(-c'm - (1/2)c'Rc) to account for power gains from the controls. You
* can estimate k by estimating the mean and variance-covariance matrix of Z,
* then by running a pilot Poisson regression to extract the projected
* coefficients. Then calculate k and pass the result to phi.

* Overdispersion and additional controls can be accounted by calculating s^2 and
* k as above and then passing (s^2 * k) to phi.

* Usage:

* estimate sample size
* power poisson a b0 ba [, power(numlist) options ]

* estimate power
* power poisson a b0 ba, n(numlist) [ options ]

* estimate effect size
* power poisson a b0, n(numlist) power(numlist) [ options ]

* a is the constant term
* b0 is the coefficient under the null hypothesis
* ba is the coefficient under the alternative hypothesis

* Options

* This is a user-written power command, so it supports many standard stata
* power options. Type 'help power' for more. The below options are particular
* this program, but others may be supported if they are implemented by stata's
* wrappers around this command (e.g. 'graph')

* alpha(numlist) - significance level. default is alpha(0.05)
* note that all tests are two-sided. to get the one-sided test, double the
* significance level.

* power(numlist) - power. default is power(0.8). required to compute effect size

* beta(numlist) - type II error rate (1-power). can be specified instead of
* power(). default is beta(0.2)

* n(numlist) - total sample size. required to compute power or effect size

* mu_t(real) - mean of exposure parameter. default is mu_t(1)

* pi(real) - mean of binary covariate. default is pi(0.5)

* phi(real) - poisson overdispersion parameter (see sigma^2 in power notes)
* and/or multiplier to account for power gains from statistical controls (see
* kappa in power notes).  default is phi(1) (no overdispersion/no gains from
* statistical controls)

* direction(upper|lower) - when graphing effect sizes, graph the minimum
* detectable effect ba > b0 (upper) or ba < b0 (lower). default is upper.

* options for the optimizer when estimating effect size

* opt_lower(real) - when estimating effect size, the inner range for the
* optimizer will be b0 + opt_lower for minimum detectable effect > b0 and
* b0 - opt_lower for minimum detectable effect < b0. default opt_lower(0.001)

* opt_upper(real) - when estimating effect size, the outer range for the
* optimizer will be b0 + opt_upper for minimum detectable effect > b0 and
* b0 - opt_upper for minimum detectable effect < b0. default opt_upper(10)

* trace - when estimating effect size, enable tracing of the optimizer

// initializer
capture program drop power_cmd_poisson_init
program power_cmd_poisson_init, sclass

	syntax [, ///
		Power(string) Beta(string) ///
		n(string) direction(string) ///
		* ///
	]
	
	
	sreturn clear
	//sreturn local pss_numopts "STDDiff"
	
	if ( (!missing("`power'") | !missing("`beta'")) & !missing("`n'") ) {
		* effect size mode
		sreturn local pss_colnames "b0 ba_lower ba_upper"
		
		* set target parameter - this is what will be graphed with the
		* 'graph' option
		if (missing("`direction'") | "`direction'"=="upper") {
			sreturn local pss_target "ba_upper"
		}
		else {
			sreturn local pss_target "ba_lower"
		}
	}
	else {
		* sample size or power mode
		sreturn local pss_colnames "b0 ba"
	}
	
end

capture program drop power_cmd_poisson
program define power_cmd_poisson, rclass
	
	syntax anything [, ///
		Alpha(string) ///
		Power(string) Beta(string) ///
		n(string) ///
		phi(real 1) ///
		mu_t(real 1) ///
		pi(real 0.5) ///
		ONESIDed DIRection(string) ///
		opt_lower(real 0) opt_upper(real 0) trace ///
	]
	
//	it currently ignores the onesided & direction options! why??
	
	* in all specifications, either 2 or 3 command arguments required
	local anything_count : word count `anything'
	if (`anything_count'>3) {
		display as error "too many command arguments specified"
		error 198
		exit
	}		
	else if (`anything_count'<2) {
		display as error "too few command arguments specified"
		error 198
		exit
	}
	
	* figure out if we're estimating sample size, power, or effect size
	if ( `anything_count'==3 & missing("`n'") ) {
		local mode "sample size"
	}
	else if ( `anything_count'==3 & missing("`power'") & missing("`beta'") & !missing("`n'") ) {
		local mode "power"
	}
	else if ( `anything_count'==2 & (!missing("`power'") | !missing("`beta'")) & !missing("`n'") ) {
		local mode "effect size"
	}
	else {
		display as error "please choose sample size, power, or effect size syntax"
		error 198
		exit
	}
	
	* extract a, b0, and (if needed) ba
	local a : word 1 of `anything'
	local b0 : word 2 of `anything'
	if (inlist("`mode'","sample size","power")) {
		local ba : word 3 of `anything'
		
		/* direction option not being passed to the program ... why?!
		if (!missing("`direction'")) {
			display as error "direction only valid for effect size syntax"
			error 198
		}
		*/
	}
	
	if (inlist("`mode'","sample size","effect size")) {
		if (!missing("`power'")) {
			local beta = 1-`power'
		}

		if (!missing("`beta'")) {
			local power = 1-`beta'
		}
	}

	local z_alpha = invnormal(1-(`alpha'/2))
	//display "z_alpha: `z_alpha'"
	local null_part = `z_alpha'*sqrt(1/(1-`pi')+(1/`pi')*exp(-`b0'))
	//display "null part: `null_part'"
	
	if ("`mode'"=="sample size") {

		local z_beta = invnormal(1-`beta')
		//display "z_beta: `z_beta'"
		local alt_part = `z_beta'*sqrt(1/(1-`pi')+(1/`pi')*exp(-`ba'))
		//display "alt part: `alt_part'"
		
		//display "pfx part"
		//display exp(-`a')*`phi'*(`mu_t'^-1)*((`ba'-`b0')^-2)
	
		local n = exp(-`a')*`phi'*(`mu_t'^-1)*((`ba'-`b0')^-2)*((`null_part'+`alt_part')^2)
		//display "n: `n'"
		
	}
	else if ("`mode'"=="power") {
		
		local numerator = sqrt((`n'*`mu_t'*((`ba'-`b0')^2))/(`phi'*exp(-`a'))) - `null_part'
		local denominator = sqrt(1/(1-`pi')+(1/`pi')*exp(-`ba'))
		
		local power = normal(`numerator'/`denominator')
		local beta = 1-`power'
		
	}
	else if ("`mode'"=="effect size") {
		
		/* USELESS BECAUSE STATA WON'T PASS direction OPTION TO COMMAND
		if (missing("`direction'") | "`direction'"=="upper") {
			local opt_lower = cond(`opt_lower'==0,.0001,`opt_lower')
			local opt_upper = cond(`opt_upper'==0,1000,`opt_upper')
		}
		else {
			local opt_lower = cond(`opt_lower'==0,-1000,`opt_lower')
			local opt_upper = cond(`opt_upper'==0,-.0001,`opt_upper')
		}
		*/

		* set optimizer extents
		
		local opt_lower = cond(`opt_lower'==0,.0001,`opt_lower')
		local opt_upper = cond(`opt_upper'==0,10,`opt_upper')
		
		* estimate effect size for ba > b0
		local range_start = `b0'+`opt_lower'
		local range_end = `b0'+`opt_upper'			
		minbound _poisson_minbound, ///
			range(`range_start' `range_end') ///
			arguments(`a' `b0' `alpha' `beta' `pi' `phi' `mu_t' `n') ///
			`trace'
		local ba_upper = r(x)
		
		* estimate effect size for ba < b0
		local range_start = `b0'-`opt_upper'
		local range_end = `b0'-`opt_lower'		
		minbound _poisson_minbound, ///
			range(`range_start' `range_end') from(`b0'-`ba_upper') ///
			arguments(`a' `b0' `alpha' `beta' `pi' `phi' `mu_t' `n') ///
			`trace'
		local ba_lower = r(x)
		
	}
	
	return scalar b0 = `b0'
	
	if ("`mode'"=="effect size") {
		return scalar ba_lower = `ba_lower'		
		return scalar ba_upper = `ba_upper'		
	}
	else {
		return scalar ba = `ba'
	}
	
	return scalar power = `power'
	return scalar beta = `beta'
	return scalar N = `n'
	return scalar alpha = `alpha'

end

// optimizer program for estimating effect size
capture program drop _poisson_minbound
program define _poisson_minbound, rclass
	
	syntax anything
	
	args ba a b0 alpha beta pi phi mu_t n
	
	local z_alpha = invnormal(1-(`alpha'/2))
	local null_part = `z_alpha'*sqrt(1/(1-`pi')+(1/`pi')*exp(-`b0'))

	local z_beta = invnormal(1-`beta')
	local alt_part = `z_beta'*sqrt(1/(1-`pi')+(1/`pi')*exp(-`ba'))
	
	local eqn = `n' - exp(-`a')*`phi'*(`mu_t'^-1)*((`ba'-`b0')^-2)*((`null_part'+`alt_part')^2)
	return scalar fx = (`eqn')^2
	
	
end

