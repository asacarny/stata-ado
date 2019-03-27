* 1.0 Adam Sacarny March 2019
* randomization inference confidence intervals

***** USE AT YOUR OWN RISK, THIS CODE HAS NOT BEEN TESTED FOR GENERAL USE ******

* usage: rici , ritest(string) reg(string) treatvar(varname) outcome(varname)
*
* ritest( ) - full ritest command before the colon, including the command name
* e.g.: ritest treat _b[treat], reps(100)
* 
* reg( ) - regression command to pass to ritest after the colon
* e.g.: regress y treat
*
* treatvar( ) - treatment variable of interest to permute, e.g. treat
*
* outcome( ) - outcome variable, e.g. y

program define rici, rclass
	syntax [if], ritest(string) reg(string) treatvar(varname) outcome(varname)
	
	// run regression to collect point estimate and asymptotic s.e.
	`reg'
	
	local b = _b[`treatvar']
	local se = _se[`treatvar']
	
	// the region of non-rejection
	local ci_upper = .
	local ci_lower = .
	
	// the innermost region of rejection (we know the CI does not go to here)
	local ci_upper_limit = .
	local ci_lower_limit = .
	
	// search for upper and lower CI limit
	
	foreach direction in up down {
	
		display "*** going `direction' ***"
		
		if ("`direction'"=="up") {
			local plusminus "+"
			local arithmetic "Adding"
		}
		else {
			local plusminus "-"
			local arithmetic "Subtracting"
		}
	
		// for up: upwards search until we get a rejected test
		// for down: downwards search until we get a rejected test
	
		local s_outer = `b' `plusminus' 2*`se'
		local pval = 1
	
		display "Initializing search with s_outer `s_outer'"
	
		while (`pval'>0.05) {
			display "Searching with s_outer `s_outer'"
			`ritest' null(`outcome' `s_outer') : `reg'
			mat pmat = r(p)
			local pval = pmat[1,1]
		
			if (`pval'>=0.05) {
				display "Failed to reject with p=" string(`pval')
				display "`arithmetic' " string(2*`se') " to s_outer"
				local s_outer = `s_outer' `plusminus' 2*`se'
			}
			else {
				display "Rejected with p=" string(`pval')
			}
		}
	
		display "Found valid s_outer: " string(`s_outer')
	
		// bisection search
		
		// we progressively divide up an interval (inner, outer)
		// inner is a value we know does not reject
		// outer is a value we know rejects
		// we initialize the algorithm with inner = b, outer=s_outer
		// then we test the midpoint of the interval to see if it rejects
		
		// ganong and jaeger def of converged:
		// 2 consecutive midpoints < 1/10 the asymptotic CI
		
		// our approach:
		// interval range (where we know outer rejects and inner does not)
		// is < 1/10 the asymptotic s.e.
		
		local TOL = `se'/10
	
		local interval_inner = `b'
		local interval_outer = `s_outer'
		local pval = 1
	
		while ( abs(`interval_outer' - `interval_inner') > `TOL') {
			display "Searching with interval range (" string(`interval_inner') ", " string(`interval_outer') ")"

			local midpoint = (`interval_inner' + `interval_outer')/2
			display "Testing at midpoint " string(`midpoint')

			`ritest' null(`outcome' `midpoint') : `reg'
			mat pmat = r(p)
			local pval = pmat[1,1]
		
			if (`pval'<0.05) {
				display "Rejected at midpoint, midpoint is new outer end"
				local interval_outer = `midpoint'
			}
			else {
				display "Failed to reject at midpoint, midpoint is new inner end"
				local interval_inner = `midpoint'
			}

			display "New interval range (" string(`interval_inner') ", " string(`interval_outer') ")"

		}
		
		display "Exited with interval range (" string(`interval_inner') ", " string(`interval_outer') ")"
		display "CI limit: " string(`interval_inner')
		
		if ("`direction'"=="up") {
			// upper end of CI
			
			// highest value above b that we know fails to reject
			local ci_upper `interval_inner'
			// lowest value above b that we know rejects
			local ci_upper_limit `interval_outer'
		}
		else {
			// lower end of CI
		
			// lowest value below b that we know fails to reject
			local ci_lower `interval_inner'
			// highest value below b that we know rejects
			local ci_lower_limit `interval_outer'
		}
	}
	
	display "Finished!"
	display "Beta: " string(`b')
	display "CI: (" string(`ci_lower') ", " string(`ci_upper') ")"
	
	return scalar beta = `b'
	return scalar ci_lower = `ci_lower'
	return scalar ci_upper = `ci_upper'
	return scalar ci_lower_limit = `ci_lower_limit'
	return scalar ci_upper_limit = `ci_upper_limit'
	
end
