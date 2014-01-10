*! version 6.01  4oct2013  Michael Stepner, stepner@mit.edu

/* CC0 license information:
To the extent possible under law, the author has dedicated all copyright and related and neighboring rights
to this software to the public domain worldwide. This software is distributed without any warranty.

This code is licensed under the CC0 1.0 Universal license.  The full legal text as well as a
human-readable summary can be accessed at http://creativecommons.org/publicdomain/zero/1.0/
*/

* Why did I include a formal license? Jeff Atwood gives good reasons: http://www.codinghorror.com/blog/2007/04/pick-a-license-any-license.html


program define binscatter, eclass
	version 11
	
	syntax varlist(min=2 numeric) [if] [in] [aweight fweight pweight], [by(varname) ///
		Nquantiles(integer 20) GENxq(name) discrete xq(varname numeric) ///
		CONTROLs(varlist numeric ts fv) absorb(varname) NOAddmean ///
		LINEtype(string) rd(numlist ascending) reportreg ///
		COLors(string) MColors(string) LColors(string) Msymbols(string) ///
		savegraph(string) savedata(string) replace ///
		nofastxtile randvar(varname numeric) randcut(real 1) randn(integer -1) ///
		/* LEGACY OPTIONS */ nbins(integer 20) create_xq x_q(varname numeric) symbols(string) method(string) unique(string) ///
		*]

	set more off

	* Create convenient weight local
	if ("`weight'"!="") local wt [`weight'`exp']
	
	
	***** Begin legacy option compatibility code
	
	if (`nbins'!=20) {
		if (`nquantiles'!=20) {
			di as error "Cannot specify both nquantiles() and nbins(): both are the same option, nbins is supported only for backward compatibility."
			exit
		}
		di as text "NOTE: legacy option nbins() has been renamed nquantiles(), and is supported only for backward compatibility."
		local nquantiles=`nbins'
	}
	
	if ("`create_xq'"!="") {
		if ("`genxq'"!="") {
			di as error "Cannot specify both genxq() and create_xq: both are the same option, create_xq is supported only for backward compatibility."
			exit
		}
		di as text "NOTE: legacy option create_xq has been renamed genxq(), and is supported only for backward compatibility."
		local genxq="q_"+word("`varlist'",-1)
	}
	
	if ("`x_q'"!="") {
		if ("`xq'"!="") {
			di as error "Cannot specify both xq() and x_q(): both are the same option, x_q() is supported only for backward compatibility."
			exit
		}
		di as text "NOTE: legacy option x_q() has been renamed xq(), and is supported only for backward compatibility."
		local xq `x_q'
	}
	
	if ("`symbols'"!="") {
		if ("`msymbols'"!="") {
			di as error "Cannot specify both msymbols() and symbols(): both are the same option, symbols() is supported only for backward compatibility."
			exit
		}
		di as text "NOTE: legacy option symbols() has been renamed msymbols(), and is supported only for backward compatibility."
		local msymbols `symbols'
	}
	
	if ("`linetype'"=="noline") {
		di as text "NOTE: legacy line type 'noline' has been renamed 'none', and is supported only for backward compatibility."
		local linetype none
	}
	
	if ("`method'"!="") {
		di as text "NOTE: method() is no longer a recognized option, and will be ignored. binscatter now uses method(log) but without a need for two instances"
	}
	
	if ("`unique'"!="") {
		di as text "NOTE: unique() is no longer a recognized option, and will be ignored. binscatter now considers the x-variable discrete if it has fewer unique values than nquantiles()"
	}
		
	***** End legacy option capatibility code

	*** Perform checks

	* Check that binscatter isn't being run quietly
	if !c(noisily) {
		di as error "binscatter cannot be run quietly"
		exit
	}

	* Set default linetype and check valid
	if ("`linetype'"=="") local linetype lfit
	else if !inlist("`linetype'","connect","lfit","qfit","none") {
		di as error "linetype() must either be connect, lfit, qfit, or none"
		exit
	}
	
	* Check that nofastxtile isn't combined with fastxtile-only options
	if "`fastxtile'"=="nofastxtile" & ("`randvar'"!="" | `randcut'!=1 | `randn'!=-1) {
		di as error "Cannot combine randvar, randcut or randn with nofastxtile"
		exit
	}

	* Misc checks
	if ("`genxq'"!="" & ("`xq'"!="" | "`discrete'"!="")) | ("`xq'"!="" & "`discrete'"!="") {
		di as error "Cannot specify more than one of genxq(), xq(), and discrete simultaneously."
		exit
	}
	if ("`genxq'"!="") confirm new variable `genxq'
	if ("`xq'"!="") {
		capture assert `xq'==int(`xq') & `xq'>0
		if _rc!=0 {
			di as error "xq() must contain only positive integers."
			exit
		}
	}
	if `nquantiles'!=20 & ("`xq'"!="" | "`discrete'"!="") {
		di as error "Cannot specify nquantiles in combination with discrete or an xq variable."
		exit
	}
	if "`reportreg'"!="" & !inlist("`linetype'","lfit","qfit") {
		di as error "Cannot specify 'reportreg' when no fit line is being created."
		exit
	}
	if "`replace'"=="" {
		if `"`savegraph'"'!="" {
			if regexm(`"`savegraph'"',"\.[a-zA-Z0-9]+$") confirm new file `"`savegraph'"'
			else confirm new file `"`savegraph'.gph"'
		}
		if `"`savedata'"'!="" {
			confirm new file `"`savedata'.csv"'
			confirm new file `"`savedata'.do"'
		}
	}

	* Mark sample (reflects the if/in conditions, and includes only nonmissing observations)
	marksample touse
	markout `touse' `by' `xq' `controls' `absorb', strok

	* Parse varlist into y-vars and x-var
	local x_var=word("`varlist'",-1)
	local y_vars=regexr("`varlist'"," `x_var'$","")
	local ynum=wordcount("`y_vars'")

	* Check number of unique byvals & create local storing byvals
	if "`by'"!="" {
		capture confirm numeric variable `by'
		if !_rc { 
			* by-variable is numeric => use tab
			tempname byvalmatrix
			qui tab `by' if `touse', nofreq matrow(`byvalmatrix')
			
			local bynum=r(r)
			forvalues i=1/`bynum' {
				local byvals `byvals' `=`byvalmatrix'[`i',1]'
			}
		}
		else {
			* by-variable is string => use levelsof (which is much slower)
			quietly levelsof `by' if `touse'
			local byvals `"`r(levels)'"'
			local bynum : word count `"`byvals'"'
		}
	}
	else local bynum=1

	*** Create residuals (if needed)
	
	if (`"`controls'`absorb'"'!="") quietly {
	
		* Parse absorb to define the type of regression to be used
		if `"`absorb'"'!="" {
			local regtype "areg"
			local absorb "absorb(`absorb')"
		}
		else {
			local regtype "reg"
		}
	
		* Generate residuals
		
		local firstloop=1
		foreach var of varlist `x_var' `y_vars' {
			tempvar residvar
			`regtype' `var' `controls' `wt' if `touse', `absorb'
			predict `residvar' if e(sample), residuals
			label variable `residvar' "`var'"
			if ("`noaddmean'"=="") {
				summarize `var' `wt' if `touse', meanonly
				replace `residvar'=`residvar'+r(mean)
			}
			
			if `firstloop'==1 {
				local x_r `residvar'
				local firstloop=0
			}
			else local y_vars_r `y_vars_r' `residvar'
		}
		
	}
	else { 	/*absorb and controls both empty, no need for regression*/
		local x_r `x_var'
		local y_vars_r `y_vars'
	}


	*** Regressions for fit lines
	if ("`reportreg'"=="") local reg_verbosity "quietly"

	if inlist("`linetype'","lfit","qfit") `reg_verbosity' {

		* If doing a quadratic fit, generate a quadratic term in x
		if "`linetype'"=="qfit" {
			tempvar x_r2
			gen `x_r2'=`x_r'^2
		}
		
		* Create matrices to hold regression results
		tempname e_b_temp
		forvalues i=1/`ynum' {
			tempname y`i'_coefs
		}
		
		* LOOP over by-vars
		local counter_by=1
		if ("`by'"=="") local noby="noby"
		foreach byval in `byvals' `noby' {
		
			* LOOP over rd intervals
			tokenize  "`rd'"
			local counter_rd=1	
				
			while ("`1'"!="" | `counter_rd'==1) {
			
				* display text headers
				if "`reportreg'"!="" {
					di "{txt}{hline}"
					if ("`by'"!="") di "-> `by' = `byval'"
					if ("`rd'"!="") {
						if (`counter_rd'==1) di "RD: x<=`1'"
						else if ("`2'"!="") di "RD: x>`1' & x<=`2'"
						else di "RD: x>`1'"
					}
				}
				
				* set conditions on reg
				local conds `touse'
				
				if ("`by'"!="") local conds `conds' & `by'==`byval'
				
				if ("`rd'"!="") {
					if (`counter_rd'==1) local conds `conds' & `x_r'<=`1'
					else if ("`2'"!="") local conds `conds' & `x_r'>`1' & `x_r'<=`2'
					else local conds `conds' & `x_r'>`1'
				}

				* LOOP over y-vars
				local counter_depvar=1
				foreach depvar of varlist `y_vars_r' {
				
					* display text headers
					if (`ynum'>1) {
						if ("`controls'`absorb'"!="") local depvar_name : var label `depvar'
						else local depvar_name `depvar'
						di as text "{bf:y_var = `depvar_name'}"
					}
					
					* perform regression
					if ("`reg_verbosity'"=="quietly") capture reg `depvar' `x_r2' `x_r' `wt' if `conds'
					else capture noisily reg `depvar' `x_r2' `x_r' `wt' if `conds'
					
					* store results
					if (_rc==0) matrix e_b_temp=e(b)
					else if (_rc==2000) {
						if ("`reg_verbosity'"=="quietly") di as error "no observations for one of the fit lines. add 'reportreg' for more info."
						
						if ("`linetype'"=="lfit") matrix e_b_temp=.,.
						else /*("`linetype'"=="qfit")*/ matrix e_b_temp=.,.,.
					}
					else {
						error _rc
						exit _rc
					}
					
					* relabel matrix row			
					if ("`by'"!="") matrix roweq e_b_temp = "by`counter_by'"
					if ("`rd'"!="") matrix rownames e_b_temp = "rd`counter_rd'"
					else matrix rownames e_b_temp = "="
					
					* save to y_var matrix
					if (`counter_by'==1 & `counter_rd'==1) matrix `y`counter_depvar'_coefs'=e_b_temp
					else matrix `y`counter_depvar'_coefs'=`y`counter_depvar'_coefs' \ e_b_temp
					
					* increment depvar counter
					local ++counter_depvar
				}
			
				* increment rd counter
				if (`counter_rd'!=1) mac shift
				local ++counter_rd
				
			}
			
			* increment by counter
			local ++counter_by
			
		}
	
		* relabel matrix column names
		forvalues i=1/`ynum' {
			if ("`linetype'"=="lfit") matrix colnames `y`i'_coefs' = "`x_var'" "_cons"
			else if ("`linetype'"=="qfit") matrix colnames `y`i'_coefs' = "`x_var'^2" "`x_var'" "_cons"
		}
	
	}

	*** Generate variable containing bins of x-var
	
	local force_discrete=0
	
	* Specify and/or create the xq var, as necessary
	if "`xq'"=="" {
	
		if "`discrete'"=="" { /* xq() and discrete are not specified */
			
			* Check whether the number of unique values > nquantiles, or <= nquantiles
			capture tab `x_r' if `touse', nofreq
			
			if (_rc==0 & r(r)<=`nquantiles') { /* number of unique values <= nquantiles, set to discrete */
				local discrete discrete
				local force_discrete=1
				if "`genxq'"!="" {
					di as text `"note: the x-variable has fewer unique values than the number of bins specified (`nquantiles').  It will therefore be treated as discrete, and genxq() will be ignored"'
				}
			}
			else if (_rc==134 | r(r)>`nquantiles') { /* number of unique values > nquantiles, perform binning */
				if ("`genxq'"!="") local xq `genxq'
				else tempvar xq
	
				if ("`fastxtile'"!="nofastxtile") fastxtile `xq' = `x_r' `wt' if `touse', nq(`nquantiles') randvar(`randvar') randcut(`randcut') randn(`randn')
				else xtile `xq' = `x_r' `wt' if `touse', nq(`nquantiles')
			}
			else {
				error _rc
			}
			
		}
		
		if "`discrete'"!="" { /* discrete is specified, xq() & genxq() are not. note that we don't use 'else' here because discrete could be set in the previous conditional */
		
			if ("`controls'`absorb'"!="" & `force_discrete'!=1) di as text "warning: discrete is specified in combination with controls() or absorb(). note that binning takes places after residualization, so the residualized x-variable may contain many more unique values."
		
			local xq `x_var'
			tempname xbin_means
			* set nquantiles var
			qui tab `xq' if `touse', nofreq matrow(`xbin_means')
			local nquantiles=r(r)
		}
	}
	else {
		if ("`controls'`absorb'"!="") di as text "warning: xq() is specified in combination with controls() or absorb(). note that binning takes places after residualization, so the xq variable should contain bins of the residuals."
		
		* set nquantiles var
		qui sum `xq' if `touse', meanonly
		local nquantiles=r(max)
	}


	**************************************
	
	* if not discrete, create a vector with the x-var mean for each bin
	if "`discrete'"=="" {
		__tab_sum_parse `xq' if `touse' `wt', sum(`x_r') means wrap nolabel noobs
	
		tempname xbin_means	
		if r(r)==`nquantiles' {
			matrix `xbin_means'=r(yval)
			local xbin_min=`xbin_means'[1,1]
			local xbin_max=`xbin_means'[`nquantiles',1]
		}
		else {
			* there were fewer quantiles than `nquantiles'.  this can happen with xtile.  ex: "sysuse auto; xtile test=mpg, nq(20); tab test"
			tempname r_xval r_yval
			matrix `r_xval'=r(xval)
			matrix `r_yval'=r(yval)
			
			local xbin_min=`r_yval'[1,1]
			local xbin_max=`r_yval'[`r(r)',1]
			
			matrix `xbin_means'=J(`nquantiles',1,.)
			forvalues i=1/`r(r)' {
				matrix `xbin_means'[`r_xval'[`i',1],1]=`r_yval'[`i',1]
			}
			
		}

	}
	
	*** Create matrices containing scatter points for each y-var
	
	* LOOP over y-vars
	local counter_depvar=0
	foreach depvar of varlist `y_vars_r' {
		local ++counter_depvar
	
		* LOOP over by-vars
		local counter_by=0
		if ("`by'"=="") local noby="noby"
		foreach byval in `byvals' `noby' {
			local ++counter_by
		
			* set conditions
			local conds `touse'
			if ("`by'"!="") local conds `conds' & `by'==`byval'
			
			* compute tab
			__tab_sum_parse `xq' if `conds' `wt', sum(`depvar') means wrap nolabel noobs	
	
			* store to matrix
			if (`counter_by'==1) {
				tempname y`counter_depvar'_scatterpts
				matrix `y`counter_depvar'_scatterpts' = r(xval),r(yval)
			}
			else {
				* make matrices conformable before right appending			
				local rowdiff=rowsof(`y`counter_depvar'_scatterpts')-rowsof(r(xval))
				if (`rowdiff'==0) matrix `y`counter_depvar'_scatterpts' = `y`counter_depvar'_scatterpts',r(xval),r(yval)
				if (`rowdiff'>0) matrix `y`counter_depvar'_scatterpts' = `y`counter_depvar'_scatterpts', ( (r(xval),r(yval)) \ J(`rowdiff',2,.) )
				if (`rowdiff'<0) matrix `y`counter_depvar'_scatterpts' = ( `y`counter_depvar'_scatterpts' \ J(-`rowdiff',colsof(`y`counter_depvar'_scatterpts'),.) ) ,r(xval),r(yval)
			}
			
		}
		
	}
	

	*********** Perform Graphing ***********

	* If rd is specified, prepare xline parameters
	if "`rd'"!="" {
		foreach xval in "`rd'" {
			local xlines `xlines' xline(`xval', lpattern(dash) lcolor(gs8))
		}
	}

	* Fill colors if missing
	if `"`colors'"'=="" local colors ///
		navy maroon forest_green dkorange teal cranberry lavender ///
		khaki sienna emidblue emerald brown erose gold bluishgray ///
		/* lime magenta cyan pink blue */
	if `"`mcolors'"'=="" {
		if (`ynum'==1 & `bynum'==1 & "`linetype'"!="connect") local mcolors `: word 1 of `colors''
		else local mcolors `colors'
	}
	if `"`lcolors'"'=="" {
		if (`ynum'==1 & `bynum'==1 & "`linetype'"!="connect") local lcolors `: word 2 of `colors''
		else local lcolors `colors'
	}
	local num_mcolor=wordcount(`"`mcolors'"')
	local num_lcolor=wordcount(`"`lcolors'"')


	* Prepare connect & msymbol options
	if ("`linetype'"=="connect") local connect "c(l)"
	if "`msymbols'"!="" {
		local symbol_prefix "msymbol("
		local symbol_suffix ")"
	}
	
	*** Prepare scatters
	
	* c indexes which color is to be used
	local c=0
	
	local counter_series=0
	
	* LOOP over by-vars
	local counter_by=0
	if ("`by'"=="") local noby="noby"
	foreach byval in `byvals' `noby' {
		local ++counter_by
		
		local xind=`counter_by'*2-1
		local yind=`counter_by'*2

		* LOOP over y-vars
		local counter_depvar=0
		foreach depvar of varlist `y_vars' {
			local ++counter_depvar
			local ++c
			
			* LOOP over rows (each row contains a coordinate pair)
			local row=1
			local xval=`y`counter_depvar'_scatterpts'[`row',`xind']
			local yval=`y`counter_depvar'_scatterpts'[`row',`yind']
			
			if !missing(`xval',`yval') {
				local ++counter_series
				local scatters `scatters' (scatteri
				if ("`savedata'"!="") {
					if ("`by'"=="") local savedata_scatters `savedata_scatters' (scatter `depvar' `x_var'
					else local savedata_scatters `savedata_scatters' (scatter `depvar'_by`counter_by' `x_var'
				}
			}
			else {
				* skip the rest of this loop iteration
				continue
			}
			
			while (`xval'!=. & `yval'!=.) {
			
				* If the x-var isn't discrete, `xval' currently contains the bin #, and we need to fetch the corresponding mean
				if (`"`discrete'"'=="") local xval=`xbin_means'[`xval',1]
				
				local scatters `scatters' `yval' `xval'
			
				local ++row
				local xval=`y`counter_depvar'_scatterpts'[`row',`xind']
				local yval=`y`counter_depvar'_scatterpts'[`row',`yind']
			}
			
			* Add options
			local scatter_options `connect' mcolor(`: word `c' of `mcolors'') lcolor(`: word `c' of `lcolors'') `symbol_prefix'`: word `c' of `msymbols''`symbol_suffix'
			local scatters `scatters', `scatter_options')
			if ("`savedata'"!="") local savedata_scatters `savedata_scatters', `scatter_options')
		

			* Add legend
			if "`by'"=="" {
				if (`ynum'==1) local legend_labels off
				else local legend_labels `legend_labels' lab(`counter_series' `depvar')
			}
			else {
				if (`ynum'==1) local legend_labels `legend_labels' lab(`counter_series' `by'=`byval')			
				else local legend_labels `legend_labels' lab(`counter_series' `depvar': `by'=`byval')
			}
			if ("`by'"!="" | `ynum'>1) local order `order' `counter_series'
			
		}
		
	}
	
	*** Fit lines
		
	if inlist(`"`linetype'"',"lfit","qfit") {
	
		* c indexes which color is to be used
		local c=0
		
		local rdnum=wordcount("`rd'")+1
		
		tempname fitline_bounds
		if ("`discrete'"=="") {
			if ("`rd'"=="") matrix `fitline_bounds'=`xbin_min', `xbin_max'
			else matrix `fitline_bounds'=`xbin_min', `=subinstr("`rd'"," ",",",.)', `xbin_max'
		}
		else {
			quietly sum `xq' if `touse', meanonly
			
			if ("`rd'"=="") matrix `fitline_bounds'=r(min), r(max)
			else matrix `fitline_bounds'=r(min), `=subinstr("`rd'"," ",",",.)', r(max)
		}

		* LOOP over by-vars
		local counter_by=0
		if ("`by'"=="") local noby="noby"
		foreach byval in `byvals' `noby' {
			local ++counter_by
			
			* Set the row to start seeking from
			*     note: each time we seek a coeff, it should be from row (rd_num)(counter_by-1)+counter_rd
			local row0=( `rdnum' ) * (`counter_by' - 1)
			
			
			* LOOP over y-vars
			local counter_depvar=0
			foreach depvar of varlist `y_vars_r' {
				local ++counter_depvar
				local ++c
		
				* LOOP over rd intervals
				forvalues counter_rd=1/`rdnum' {
					
					if (`"`linetype'"'=="lfit") {
						local coef_quad=0
						local coef_lin=`y`counter_depvar'_coefs'[`row0'+`counter_rd',1]
						local coef_cons=`y`counter_depvar'_coefs'[`row0'+`counter_rd',2]
					}
					else if (`"`linetype'"'=="qfit") {
						local coef_quad=`y`counter_depvar'_coefs'[`row0'+`counter_rd',1]
						local coef_lin=`y`counter_depvar'_coefs'[`row0'+`counter_rd',2]
						local coef_cons=`y`counter_depvar'_coefs'[`row0'+`counter_rd',3]
					}
					
					if !missing(`coef_quad',`coef_lin',`coef_cons') {
						local leftbound=`fitline_bounds'[1,`counter_rd']
						local rightbound=`fitline_bounds'[1,`counter_rd'+1]
					
						local fits `fits' (function `coef_quad'*x^2+`coef_lin'*x+`coef_cons', range(`leftbound' `rightbound') lcolor(`: word `c' of `lcolors''))
					}
				}
			}
		}
	}
	
	* Prepare y-axis title
	if (`ynum'==1) local ytitle `y_vars'
	else if (`ynum'==2) local ytitle : subinstr local y_vars " " " and "
	else local ytitle : subinstr local y_vars " " "; ", all

	* Display graph
	local graphcmd twoway `scatters' `fits', graphregion(fcolor(white)) `xlines' xtitle(`x_var') ytitle(`ytitle') legend(`legend_labels' order(`order')) `options'
	if ("`savedata'"!="") local savedata_graphcmd twoway `savedata_scatters' `fits', graphregion(fcolor(white)) `xlines' xtitle(`x_var') ytitle(`ytitle') legend(`legend_labels' order(`order')) `options'
	`graphcmd'
	
	* Save graph
	if `"`savegraph'"'!="" {
		* check file extension using a regular expression
		if regexm(`"`savegraph'"',"\.[a-zA-Z0-9]+$") local graphextension=regexs(0)
		
		if inlist(`"`graphextension'"',".gph","") graph save `"`savegraph'"', `replace'
		else graph export `"`savegraph'"', `replace'
	}

	* Save data
	if ("`savedata'"!="") {
	
		*** Save a CSV containing the scatter points
	
		file open __savedata using `"`savedata'.csv"', write text `replace'
		
		* LOOP over rows
		forvalues row=0/`nquantiles' {
		
			* x-var
			if (`row'==0) file write __savedata "`x_var'"
			else {
				local cur_xval=`xbin_means'[`row',1]
				file write __savedata (`cur_xval')
			}
			
			* LOOP over y-vars
			local counter_depvar=0
			foreach depvar of varlist `y_vars' {
				local ++counter_depvar

				* LOOP over by-vals
				forvalues counter_by=1/`bynum' {
				
					if (`row'==0) {
						if "`by'"=="" file write __savedata ",`depvar'"
						else file write __savedata ",`depvar'_by`counter_by'"
						local shift_y`counter_depvar'_by`counter_by'=0
					}
					else {
						local cur_row=`row'+`shift_y`counter_depvar'_by`counter_by''

						if (`cur_row'<=`=rowsof(`y`counter_depvar'_scatterpts')') {
					
							* Check if x-value matches the current one being processed
							local xval=`y`counter_depvar'_scatterpts'[`cur_row',`counter_by'*2-1]
							if ("`discrete'"=="" & `xval'==`row') | ("`discrete'"!="" & `xval'==`cur_xval') {
								* x-value is correct, write the y-value to the csv
								file write __savedata "," (`y`counter_depvar'_scatterpts'[`cur_row',`counter_by'*2])
							}
							else {
								* x-value is incorrect, shift the index so that it is tried again
								local --shift_y`counter_depvar'_by`counter_by'
								file write __savedata ",."
							}
						}
						else file write __savedata ",."
					}
					
				} /* end by-val loop */
				
			} /* end y-var loop */
			
			file write __savedata _n
			
		} /* end row loop */

		file close __savedata
		di as text `"(file `savedata'.csv written containing saved data)"'
		
		
		
		*** Save a do-file with the commands to generate a nicely labeled dataset and re-create the binscatter graph
		
		file open __savedata using `"`savedata'.do"', write text `replace'
		
		file write __savedata `"insheet using `savedata'.csv"' _n _n
		
		if "`by'"!="" {
			foreach depvar of varlist `y_vars' {
				local counter_by=0
				foreach byval in `byvals' {
					local ++counter_by
					file write __savedata `"label variable `depvar'_by`counter_by' "`depvar'; `by'==`byval'""' _n
				}
			}
			file write __savedata _n
		}
		
		file write __savedata `"`savedata_graphcmd'"' _n
		
		file close __savedata
		di as text `"(file `savedata'.do written containing commands to process saved data)"'
		
	}

	*** Return items
	ereturn post, esample(`touse')
	
	qui count if e(sample)
	ereturn scalar N = r(N)
	
	ereturn local graphcmd `"`graphcmd'"'
	if inlist("`linetype'","lfit","qfit") {
		forvalues yi=`ynum'(-1)1 {
			ereturn matrix y`yi'_coefs=`y`yi'_coefs'
		}
	}
	
	if ("`rd'"!="") {
		tempname rdintervals
		matrix `rdintervals' = (. \ `=subinstr("`rd'"," ","\",.)' ) , ( `=subinstr("`rd'"," ","\",.)' \ .)

		forvalues i=1/`=rowsof(`rdintervals')' {
			local rdintervals_labels `rdintervals_labels' rd`i'
		}
		matrix rownames `rdintervals' = `rdintervals_labels'
		matrix colnames `rdintervals' = gt lt_eq
		ereturn matrix rdintervals=`rdintervals'
	}
	
	capture confirm matrix `byvalmatrix'
	if (_rc==0) {
		forvalues i=1/`=rowsof(`byvalmatrix')' {
			local byvalmatrix_labels `byvalmatrix_labels' by`i'
		}
		matrix rownames `byvalmatrix' = `byvalmatrix_labels'
		matrix colnames `byvalmatrix' = `by'
		ereturn matrix byvalues=`byvalmatrix'
	}
	
end


**********************************

* Helper programs

program define __tab_sum_parse, rclass
	version 11
	
	syntax anything(everything equalok name=cmd id="command"), [*]
	
	* Find a working tempfile to log to (fixes a Windows-only crash)
	local no_inf_loop=0
		
	while `no_inf_loop'==0 | (_rc!=0 & `no_inf_loop'<100) {
		tempfile log_file
		capture confirm new file `log_file'
		
		local ++no_inf_loop
	}
	if _rc!=0 {
		di as error "could not open a tempfile to save log"
		exit
	}
	
	* Start the log
	capture log close __tab_log
	quietly log using `log_file', name(__tab_log) replace text
	
	* Run -tab, sum()- command
	noisily tab `cmd', `options'
	
	* Close the log
	qui log close __tab_log
	
	* Open the log to be read
	capture file close __tab_log
	file open __tab_log using `log_file', read
	file read __tab_log line
	
	if "`line'"=="no observations" {
		tempname nullmatrix
		matrix `nullmatrix'=.
		
		return scalar r=0
		return matrix yval=`nullmatrix', copy
		return matrix xval=`nullmatrix'
	}
	else {
		* Find the top crossline of the table
		tokenize `"`line'"', parse("|")
		while (`"`1'"'!="------------+------------" | `"`2'"'!="") & r(eof)==0 {
			file read __tab_log line
			tokenize `"`line'"', parse("|")
		}
	
		* Load the table into the matrices until we reach the bottom crossline of the table
		tempname xval yval
		
		local firstloop=1
		
	
		file read __tab_log line
		tokenize `"`line'"', parse("|")
	
		while (`"`1'"'!="------------+------------" | `"`2'"'!="") & r(eof)==0 {
		
			if (`firstloop'==1) {
				matrix `xval'=real(subinstr("`1'",",","",.))
				matrix `yval'=real(subinstr("`3'",",","",.))
			}
			else {
				matrix `xval'=`xval' \ real(subinstr("`1'",",","",.))
				matrix `yval'=`yval' \ real(subinstr("`3'",",","",.))
			}
	
			file read __tab_log line
			tokenize `"`line'"', parse("|")
			
			local firstloop=0
		}
		
		* Return the parsed matrices with -tab, sum() results-
		return scalar r = rowsof(`xval')
		return matrix yval=`yval'
		return matrix xval=`xval'
	}
	
end

*** copy of: version 1.20  7sep2013  Michael Stepner, stepner@mit.edu
program define fastxtile, rclass
	version 11

	* Parse weights, if any
	_parsewt "aweight fweight pweight" `0' 
	local 0  "`s(newcmd)'" /* command minus weight statement */
	local wt "`s(weight)'"  /* contains [weight=exp] or nothing */

	* Extract parameters
	syntax newvarname=/exp [if] [in] [,Nquantiles(integer 2) Cutpoints(varname numeric) ALTdef ///
		CUTValues(numlist ascending) randvar(varname numeric) randcut(real 1) randn(integer -1)]

	* Mark observations which will be placed in quantiles
	marksample touse, novarlist
	markout `touse' `exp'
	qui count if `touse'
	local popsize=r(N)

	if "`cutpoints'"=="" & "`cutvalues'"=="" { /***** NQUANTILES *****/
		if `"`wt'"'!="" & "`altdef'"!="" {
			di as error "altdef option cannot be used with weights"
			exit 198
		}
		
		if `randn'!=-1 {
			if `randcut'!=1 {
				di as error "cannot specify both randcut() and randn()"
				exit 198
			}
			else if `randn'<1 {
				di as error "randn() must be a positive integer"
				exit 198
			}
			else if `randn'>`popsize' {
				di as text "randn() is larger than the population. using the full population."
				local randvar=""
			}
			else {
				local randcut=`randn'/`popsize'
				
				if "`randvar'"!="" {
					qui sum `randvar', meanonly
					if r(min)<0 | r(max)>1 {
						di as error "with randn(), the randvar specified must be in [0,1] and ought to be uniformly distributed"
						exit 198
					}
				}
			}
		}

		* Check if need to gen a temporary uniform random var
		if "`randvar'"=="" {
			if (`randcut'<1 & `randcut'>0) { 
				tempvar randvar
				gen `randvar'=runiform()
			}
			* randcut sanity check
			else if `randcut'!=1 {
				di as error "if randcut() is specified without randvar(), a uniform r.v. will be generated and randcut() must be in (0,1)"
				exit 198
			}
		}

		* Mark observations used to calculate quantile boundaries
		if ("`randvar'"!="") {
			tempvar randsample
			mark `randsample' `wt' if `touse' & `randvar'<=`randcut'
		}
		else {
			local randsample `touse'
		}

		* Error checks
		qui count if `randsample'
		local samplesize=r(N)
		if (`nquantiles' > r(N) + 1) {
			if ("`randvar'"=="") di as error "nquantiles() must be less than or equal to the number of observations [`r(N)'] plus one"
			else di as error "nquantiles() must be less than or equal to the number of sampled observations [`r(N)'] plus one"
			exit 198
		}
		else if (`nquantiles' < 2) {
			di as error "nquantiles() must be greater than or equal to 2"
			exit 198
		}

		* Compute quantile boundaries
		_pctile `exp' if `randsample' `wt', nq(`nquantiles') `altdef'

		* Store quantile boundaries in list
		forvalues i=1/`=`nquantiles'-1' {
			local cutvallist `cutvallist',`r(r`i')'
		}
	}
	else if "`cutpoints'"!="" { /***** CUTPOINTS *****/
	
		* Parameter checks
		if "`cutvalues'"!="" {
			di as error "cannot specify both cutpoints() and cutvalues()"
			exit 198
		}		
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 | `randn'!=-1 {
			di as error "cutpoints() cannot be used with nquantiles(), altdef, randvar(), randcut(), randn() or weights"
			exit 198
		}

		tempname cutvals
		qui tab `cutpoints', matrow(`cutvals')
		
		if r(r)==0 {
			di as error "cutpoints() all missing"
			exit 2000
		}
		else {
			local nquantiles = r(r) + 1
			
			forvalues i=1/`r(r)' {
				local cutvallist `cutvallist',`=`cutvals'[`i',1]'
			}
		}
	}
	else { /***** CUTVALUES *****/
		if "`wt'"!="" | "`randvar'"!="" | "`ALTdef'"!="" | `randcut'!=1 | `nquantiles'!=2 | `randn'!=-1 {
			di as error "cutvalues() cannot be used with nquantiles(), altdef, randvar(), randcut(), randn() or weights"
			exit 198
		}
		
		* parse numlist
		numlist "`cutvalues'"
		local nquantiles=wordcount(`"`r(numlist)'"')+1
		
		* put commas between each value
		local cutvallist `",`r(numlist)'"'
		local cutvallist : subinstr local cutvallist " " ",", all
	
	}

	* Create quantile variable
	qui gen `varlist'=1+irecode(`exp'`cutvallist') if `touse'
	label var `varlist' "`nquantiles' quantiles of `exp'"
	
	* Return values
	if ("`samplesize'"!="") return scalar n = `samplesize'
	else return scalar n = .
	
	return scalar N = `popsize'
	
	local cutvallist : subinstr local cutvallist "," " ", all
	tokenize `"`cutvallist'"'
	forvalues i=`=`nquantiles'-1'(-1)1 {
		return scalar r`i' = ``i''
	}

end
