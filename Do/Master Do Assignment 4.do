* Causal Inference and Research Design
* Assignment 4
* Autor: Diana Perez

clear all
set more off
cap log close
cls
graph set window fontface "Calibri" 


*********************************************************************************
** ASSIGNMENT 3																   **
*********************************************************************************
global path "C:\Users\Diana\Documents\GitHub\RDD"
cd "${path}"

* I. Github repo and summary
* -------------------------------------------------------------------------------

*) 1. Saving data
use 																			///
"https://github.com/scunning1975/causal-inference-class/raw/master/hansen_dwi", ///
clear

compress
save "${path}\Data\Hansen_dwi.dta", replace

* II. Replication 
* -------------------------------------------------------------------------------
use "${path}\Data\Hansen_dwi.dta", clear

*) 3. Eligibility variable
gen eligibility=bac1>=0.08 if !missing(bac1)
gen bac1_ajust=bac1-0.08
gen bac2_ajust=bac1_ajust*bac1_ajust

	* Labels
	label var male "Male"
	label var white "White"
	label var aged "Age"
	label var acc "Accident"
	label var bac1_ajust "DUI"

*) 4. Testing Manipulation on the RV

	* Packages
	net install rddensity, 														///
		from("https://sites.google.com/site/rdpackages/rddensity/stata") replace
	net install lpdensity, 														///
		from("https://sites.google.com/site/nppackages/lpdensity/stata") replace

	* Cattaneo et al.
	rddensity bac1_ajust /*P-value: 0.0276 */

	* Figure 1
	hist bac1, freq bin(450) bc(gs9) lc(gs9) graphregion(fcolor(white)) 		///
		ti("BAC histogram", color(black) size(vlarge) lwidth(vvthick)) 			///
		xti("BAC", size(vlarge) lwidth(vvthick)) 								///
		yti("Frequency", size(medium) lwidth(vvthick)) 							///
		addplot(pci  0 0.08 2000 0.08, lc(black)) 								///
		xlabel(0(0.1)0.4) xvarformat(%2.1f) yvarformat(%9.0gc) 					///
		ylabel(, angle(0)) legend(off)
		
	gr export "${path}\Figures\Figure 01.pdf", replace as(pdf)

*) 5. Table 2. Covariance continuity

	* Editing e(N)
	cap program drop changeN
	program define changeN, eclass
		/* This program edits the e(N). It replace it for any scalar named nobs.*/
		
		ereturn scalar N = nobs
	end

	* Regressions and table
	local bw=0.05
	local n=0
	
	global covs male white aged acc
	
	foreach var of varlist $covs {

		local ++n
		local vlab: variable label `var'
		
		* Estimation
		rdrobust `var' bac1_ajust, kernel(uniform) h(`bw' `bw') p(1) vce(hc0)
		est store rdrob
		* Output
		if "`n'"=="1" local comp="replace"
		else local comp="append"
		
		if "`n'"=="3" local j=1
		else local j=3
		
		scalar nobs=e(N_h_l)+e(N_h_r)
		local nobs=string(e(N_h_l)+e(N_h_r),"%9.0gc")
		changeN
		
		qui sum `var' if bac1_ajust>=-`bw' & bac1_ajust<0
		local mu=string(r(mean), "%5.`j'f")

		outreg2 using "${path}\Tables\Table2.tex", `comp' nocons nor2 decm(.) 	///
			dec(3)																///
			addstat(Mean at (0.079), `mu') 										///
			addtext(Controls,No) 												///
			label nonotes 														///
			addn("Standard errors are in parentheses." 							///
			"*** Significant at the 1 percent level." 							///
			"** Significant at the 5 percent level." 							///
			"* Significant at the 10 percent level.")

	}

*) 6. Figure 2.
local cut=0.08
local cut2=0.15
	
local bw=0.05
local bw2=`cut2'-`cut'
	
local bwl=`cut'-`bw'
local bwu=`cut'+`bw'
local bwl2=`cut2'-`bw'
local bwu2=`cut2'+`bw'
	
local nbinl=35
local nbinr=40
	
global covs acc male aged white 

foreach j of numlist 1/2{	

	local n=0
	foreach var of  varlist $covs {

		local ++n
		local vlab: variable label `var'

		if "`var'"=="acc" {
			local lti="Panel A. Accident at scene"
			local yax="ylabel(0.05(0.05)0.25) yscale(range(0.03 0.26))"
		}
		else if "`var'"=="male" {
			local lti="Panel B. Male"
			local yax="ylabel(0.74(0.02)0.82) yscale(range(0.73 0.83))"
		}
		else if "`var'"=="aged" {
			local lti="Panel C. Age"
			local yax="ylabel(34(1)39) yscale(range(33.5 40))"
			*local yax=" "
		}
		else if "`var'"=="white" {
			local lti="Panel D. White"
			local yax="ylabel(0.8(0.02)0.9) yscale(range(0.79 0.91))"
		}
		
		if "`n'"=="3" local l=0
		else local l=2
		
		* 2nd cutoff and above
		cap drop rdplot_*
		rdplot `var' bac1 if inrange(bac1,`bwl2',`bwu2'), binselect(es) 		///
			c(`cut2') genvars p(`j')  kernel(uniform) h(`bw' `bw2') 			///
			nbins(`nbinr' `nbinr')
		cap drop v2rdplot_mean_y v2rdplot_mean_x v2rdplot_hat_y
		rename (rdplot_mean_y rdplot_mean_x rdplot_hat_y)						///
			(v2rdplot_mean_y v2rdplot_mean_x v2rdplot_hat_y)
		
		* origin to 2nd cutoff
		cap drop rdplot_*
		rdplot `var' bac1 if inrange(bac1,`bwl',`cut2'), binselect(es) 			///
			c(`cut') genvars p(`j')  kernel(uniform) h(`bw' `bw') 				///
			nbins(`nbinl' `nbinr')

		replace rdplot_mean_y=v2rdplot_mean_y if inrange(bac1,`cut2',`bwu2')
		replace rdplot_mean_x=v2rdplot_mean_x if inrange(bac1,`cut2',`bwu2')
		replace rdplot_hat_y=v2rdplot_hat_y if inrange(bac1,`cut2',`bwu2')

		* Confidence intervals
		cap drop ep2
		cap drop sd
		cap drop rdplot_ci_l
		cap drop rdplot_ci_r
		
		gen ep2=(rdplot_mean_y-rdplot_hat_y)^2
		sum ep2 if inrange(bac1,`bwl',`cut')
		gen sd=r(mean)*r(N)/(r(N)-2) if inrange(bac1,`bwl',`cut')
		sum ep2 if inrange(bac1,`cut',`cut2')
		replace sd=r(mean)*r(N)/(r(N)-2) if inrange(bac1,`cut',`cut2')
		sum ep2 if inrange(bac1,`cut2',`bwu2')
		replace sd=r(mean)*r(N)/(r(N)-2) if inrange(bac1,`cut2',`bwu2')
		
		gen rdplot_ci_l=rdplot_hat_y-1.96*sd
		gen rdplot_ci_r=rdplot_hat_y+1.96*sd
		
		* Graphs
		if `j'==1 local q="l"
		else local q="q"

		tw (`q'fit rdplot_hat_y rdplot_mean_x if inrange(bac1,`bwl',`cut'), 	///
			lcolor(black) lwidth(medthick) lpattern(solid) xvarformat(%2.1f) 	///
			yvarformat(%3.`l'f))												///
			(`q'fit rdplot_ci_r rdplot_mean_x if inrange(bac1,`bwl',`cut'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_l rdplot_mean_x if inrange(bac1,`bwl',`cut'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_hat_y rdplot_mean_x if inrange(bac1,`cut',`cut2'),	///
			lcolor(black) lwidth(medthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_r rdplot_mean_x if inrange(bac1,`cut',`cut2'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_l rdplot_mean_x if inrange(bac1,`cut',`cut2'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_hat_y rdplot_mean_x if inrange(bac1,`cut2',`bwu2'),	///
			lcolor(black) lwidth(medthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_r rdplot_mean_x if inrange(bac1,`cut2',`bwu2'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_l rdplot_mean_x if inrange(bac1,`cut2',`bwu2'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(scatter rdplot_mean_y rdplot_mean_x, msymbol(circle_hollow) 		///
			mcolor(gs9%80)), 													///
			xline(`cut', lcolor(black) lw(medium))								///
			xline(`cut2', lcolor(black) lw(medium))								///
			xti("BAC", size(medlarge) lwidth(vvthick))							///
			yti(" ", size(medium))												///
			ti("`lti'", color(black) size(medlarge) lwidth(vvthick) 			///
			j(left) placement(nwest))											///
			xlabel(0.05 0.1 0.15 0.2) xscale(range(0.03 0.21)) `yax'			///
			graphregion(fcolor(white)) 											///
			legend(order (1 10) label (1 "Fitted") label (10 "`vlab'")) 		///
			name(G`n'_`j', replace) xsize(7) ysize(7)	
			
	}
	
	gr combine G1_`j' G2_`j' G3_`j' G4_`j', col(2) graphregion(color(white)) 	///
		ysize(14) xsize(12)
	gr export "${path}\Figures\Figure 02_P`j'.pdf", replace as(pdf)
}

*) Table 3
mat MAT_tabla=J(8,6,.)
mat MAT_tabla_s=J(8,6,0)

global controls aged white i.year male

local n=1
foreach j of numlist 0.05 0.03{
	
	* BAC
	qui reg recid eligibility bac1_ajust 										///
		${controls} if inrange(bac1_ajust,-`j',`j'), r
	
	mat MAT_tabla[`n'+1,1]=_b[eligibility]
	mat MAT_tabla[`n'+1,2]=_se[eligibility]
	local N=string(e(N),"%5.0gc")
	mat MAT_tabla[`n'+3,1]=`N'
	
	local p=string(ttail(e(df_r),abs(_b[eligibility]/_se[eligibility]))*2,		///
		"%6.4f")
	matrix MAT_tabla_s[`n'+1,1] = (`p' <= 0.1) + (`p' <= 0.05) + (`p' <= 0.01)
		
	
	qui sum recid if e(sample)
	mat MAT_tabla[`n'+2,1]=r(mean)
	
	* BAC x Elegibility
	qui reg recid eligibility bac1_ajust eligibility#c.bac1_ajust				///
		${controls} if inrange(bac1_ajust,-`j',`j'), r
		
	mat MAT_tabla[`n'+1,3]=_b[eligibility]
	mat MAT_tabla[`n'+1,4]=_se[eligibility]
	local N=string(e(N),"%5.0gc")
	mat MAT_tabla[`n'+3,3]=`N'
	
	local p=string(ttail(e(df_r),abs(_b[eligibility]/_se[eligibility]))*2,		///
		"%6.4f")
	matrix MAT_tabla_s[`n'+1,3] = (`p' <= 0.1) + (`p' <= 0.05) + (`p' <= 0.01)

	qui sum recid if e(sample)
	mat MAT_tabla[`n'+2,3]=r(mean)

	* BAC2 x Eligibility
	qui reg recid eligibility bac1_ajust bac2_ajust eligibility#c.bac1_ajust 	///
		eligibility#c.bac2_ajust 												///
		${controls} if inrange(bac1_ajust,-`j',`j'), r
		
	mat MAT_tabla[`n'+1,5]=_b[eligibility]
	mat MAT_tabla[`n'+1,6]=_se[eligibility]
	local N=string(e(N),"%5.0gc")
	mat MAT_tabla[`n'+3,5]=`N'
	
	local p=string(ttail(e(df_r),abs(_b[eligibility]/_se[eligibility]))*2,		///
		"%6.4f")
	matrix MAT_tabla_s[`n'+1,5] = (`p' <= 0.1) + (`p' <= 0.05) + (`p' <= 0.01)
	
	qui sum recid if e(sample)
	mat MAT_tabla[`n'+2,5]=r(mean)
	
	local n=`n'+4
	
}

frmttable using "${path}\Tables\Table3", replace tex 							///
	statmat(MAT_tabla) annotate(MAT_tabla_s) asymbol("*","**","***")			///
	substat(1) noblankrows														///
	ct("","Linear","Linear diferenciated","Quadratic diferenciated")			///
	rt("{\i Panel A}. BAC $\in$ [0.03,0.13]"\""\"DUI"\""\"Mean"\""\				///
	"Observations"\""\															///
	"{\i Panel B}. BAC $\in$ [0.055,0.105]"\""\"DUI"\""\"Mean"\""\				///
	"Observations")																///
	sdec(3,3,3\3,3,3\3,3,3\3,3,3\3,3,3\3,3,3\0,0,0\								///
	3,3,3\3,3,3\3,3,3\3,3,3\3,3,3\3,3,3\0,0,0)
	
*) Figure 3	
local cut=0.08
local cut2=0.15
	
local bw=0.05
local bw2=`cut2'-`cut'
	
local bwl=`cut'-`bw'
local bwu=`cut'+`bw'
	
local nbinl=35
local nbinr=40

foreach j of numlist 1/2{	

		if "`j'"=="1" local lti="{it: Panel A.} Linear polynomial"
		else if "`j'"=="2" local lti="{it: Panel B.} Quadratic polynomial"
		local yax="ylabel(0.08(0.02)0.16) yscale(range(0.07 0.17))"
						
		* origin to 2nd cutoff
		cap drop rdplot_*
		rdplot recid bac1 if inrange(bac1,`bwl',`cut2'), binselect(es) 			///
			c(`cut') genvars p(`j')  kernel(uniform) h(`bw' `bw') 				///
			nbins(`nbinl' `nbinr')

		* Confidence intervals
		cap drop ep2
		cap drop sd
		cap drop rdplot_ci_l
		cap drop rdplot_ci_r
		
		gen ep2=(rdplot_mean_y-rdplot_hat_y)^2
		sum ep2 if inrange(bac1,`bwl',`cut')
		gen sd=r(mean)*r(N)/(r(N)-2) if inrange(bac1,`bwl',`cut')
		sum ep2 if inrange(bac1,`cut',`cut2')
		replace sd=r(mean)*r(N)/(r(N)-2) if inrange(bac1,`cut',`cut2')
		
		gen rdplot_ci_l=rdplot_hat_y-1.96*sd
		gen rdplot_ci_r=rdplot_hat_y+1.96*sd
		
		* Graphs
		if `j'==1 local q="l"
		else local q="q"

		tw (`q'fit rdplot_hat_y rdplot_mean_x if inrange(bac1,`bwl',`cut'), 	///
			lcolor(black) lwidth(medthick) lpattern(solid) xvarformat(%2.1f) 	///
			yvarformat(%3.2f))												///
			(`q'fit rdplot_ci_r rdplot_mean_x if inrange(bac1,`bwl',`cut'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_l rdplot_mean_x if inrange(bac1,`bwl',`cut'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_hat_y rdplot_mean_x if inrange(bac1,`cut',`cut2'),	///
			lcolor(black) lwidth(medthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_r rdplot_mean_x if inrange(bac1,`cut',`cut2'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(`q'fit rdplot_ci_l rdplot_mean_x if inrange(bac1,`cut',`cut2'), 	///
			lcolor(black%15) lwidth(vvthick) lpattern(solid)) 					///
			(scatter rdplot_mean_y rdplot_mean_x, msymbol(circle_hollow) 		///
			mcolor(gs9%80)), 													///
			xline(`cut', lcolor(black) lw(medium))								///
			xti("BAC", size(medlarge) lwidth(vvthick))							///
			yti(" ", size(medium))												///
			ti("`lti'", color(black) size(medlarge) lwidth(vvthick) 			///
			j(left) placement(nwest))											///
			xlabel(0.05 0.1 0.15) xscale(range(0.03 0.16)) `yax'				///
			graphregion(fcolor(white)) 											///
			legend(off) 														///
			name(G_`j', replace) xsize(8) ysize(5)		

	gr export "${path}\Figures\Figure 03_P`j'.pdf", replace as(pdf)
}

	
	
	
	
	



