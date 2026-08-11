 *********************************
* Runge-Kutta- Fehlberg RKF45 method. 
*Application to Population data.P0, P1, P2 at t0, t1,t2
*Initializations of Variables and path
 ******************************
clear 
local data_path "D:\popProj\rkf45"
local eps=0.5
use "`data_path'\input_data.dta" 
************************************
* Variables Corresponding to Projection Years
************************************
 local a=1991
 local b=2031
 local yrInt= (`b'-`a')/5
 gen prYr=`a'+5
 gen recid = _n
 order recid, first 
 
 gen r1=ln(pop01/pop91)/10 
 gen r2=ln(pop11/pop01)/10
 gen r=(r1+r2)/2 
*****************************************	 
* Step Size, *Initial population 
*****************************************
	gen h=(`b'-`a')/8
	gen Yk=pop91

*****************************************
*Carrying capacity (k) which is maximum sustainable population of an area.
*****************************************
	gen k=pop11*1.5
****************************************
* Equations for Generating Coefficients
*Compute k's for first step, k1-k6
****************************************
	gen f1=r*Yk*(1-(Yk/k))
	gen k1=h*f1

	gen f2=Yk+((1/4)*k1)
	gen k2=h*f2*r*(1-(f2/k)) 

	gen f3=Yk+(3/32*k1)+(9/32*k2)
	gen k3=h*f3*r*(1-(f3/k)) 

	gen f4=Yk+((1932/2197)*k1)-((7200/2197)*k2)+((7296/2197)*k3)
	gen k4=h*f4*r*(1-(f4/k))

	gen f5=(Yk+((439/216)*k1))-(8*k2)+((3680/513)*k3)-((845/4104)*k4)
	gen k5=h*f5*r*(1-(f5/k))

	gen f6=(Yk -((8/27)*k1) +(2*k2) - ((3544/2565)*k3) + ((1859/4104) *k4) -((11/40) * k5))
	gen k6=h*f6*r*(1-(f6/k))
*************************************************
*Compute 4th order RKF45  estimates 
*************************************************
	gen yk1=Yk+(25/216*k1)+(1408/2565*k3)+(2197/4101*k4)-(1/5*k5)
*************************************************
*Compute 5th order RKF45 estimates
************************************************
	gen zk1=Yk+(16/135*k1)+(6656/12825*k3)+(28561/56430*k4)-(9/50*k5)+(2/55*k6)
***********************************************
*Error Estimation between RK4 & RK5
***********************************************
	gen R=(1/h)*abs(yk1-zk1) 
	
**********************************************
*Step size adjustment factor
**********************************************
	gen dlt=0.84*(`eps'/R)^(1/4)
	gen iter=1

save "`data_path'\kf45_computeddata.dta" , replace

keep if R>`eps'
save "`data_path'\kf45_intermediate.dta" , replace
clear
use "`data_path'\kf45_intermediate.dta" 

if R>`eps'{
		local n = 0
****************************************		
* Repeat the Procedure with Adjusted h
****************************************
while `n' >= 0{
		clear
		use "`data_path'\kf45_intermediate.dta" 
		replace Yk = zk1 
		keep recid name pop91 pop01 pop11 k r1 r2 r h dlt Yk iter prYr
		****************************************
		*Adjusted h, New optimal step size
		****************************************
		replace h=h*dlt
		drop dlt
		save "`data_path'\kf45_residuals.dta" , replace
		clear
		use "`data_path'\kf45_residuals.dta" 
		****************************************
		* Equations for Generating Coefficients
		*Compute k's for first step, k1-k6
		****************************************
		gen f1=r*Yk*(1-(Yk/k))
		gen k1=h*f1

		gen f2=Yk+((1/4)*k1)
		gen k2=h*f2*r*(1-(f2/k)) 

		gen f3=Yk+(3/32*k1)+(9/32*k2)
		gen k3=h*f3*r*(1-(f3/k)) 

		gen f4=Yk+((1932/2197)*k1)-((7200/2197)*k2)+((7296/2197)*k3)
		gen k4=h*f4*r*(1-(f4/k))

		gen f5=(Yk+((439/216)*k1))-(8*k2)+((3680/513)*k3)-((845/4104)*k4)
		gen k5=h*f5*r*(1-(f5/k))

		gen f6=(Yk -((8/27)*k1) +(2*k2) - ((3544/2565)*k3) + ((1859/4104) *k4) -((11/40) * k5))
		gen k6=h*f6*r*(1-(f6/k))
		****************************************
		*Compute 4th order RKF45  estimates
		****************************************
		gen yk1=Yk+(25/216*k1)+(1408/2565*k3)+(2197/4101*k4)-(1/5*k5)
		****************************************
		*Compute 5th order RKF45 estimates
		****************************************
		gen zk1=Yk+(16/135*k1)+(6656/12825*k3)+(28561/56430*k4)-(9/50*k5)+(2/55*k6)
		****************************************
		*Error Estimation between RK4 & RK5	
		****************************************
		gen R=(1/h)*abs(yk1-zk1)
		****************************************
		*Step size adjustment factor
		****************************************
		gen dlt=0.84*(`eps'/R)^(1/4)
		replace iter=iter+1
		save "`data_path'\kf45_residuals.dta" , replace

		use "`data_path'\kf45_computeddata.dta", clear
		append using "`data_path'\kf45_residuals.dta"
		save "`data_path'\kf45_computeddata.dta" , replace
		clear
		
		use "`data_path'\kf45_residuals.dta"
		****************************************
		*Further extending loop
		****************************************
		keep if R>`eps'
		save "`data_path'\kf45_intermediate.dta" , replace
		clear
		use "`data_path'\kf45_intermediate.dta" 
		count if !missing(recid)   // condition for existing records
		local n = r(N)

		display "Number of records found: `n'"
		
		
		if `n' == 0 {
			sleep 1000   // wait 1 second before checking again
			
			 continue, break
		}
	}	


display "Records now exist in the dataset."
	
}
use "`data_path'\kf45_computeddata.dta", clear
sort recid iter
order iter, after(recid)

save "`data_path'\kf45_finaloutput.dta", replace

*table name, statistic(max iter)

bysort recid: egen maxiter = max(iter)
keep if iter == maxiter
drop maxiter
save "`data_path'\kf45_projection1.dta", replace
save "`data_path'\kf45_collection_projections.dta" , replace

**************************************************
* Determine the next projection year.
**************************************************

local next5yr = `a'+10

use "`data_path'\kf45_projection1.dta", clear

keep recid name pop91 pop01 pop11 k r1 r2 r  zk1

save "`data_path'\kf45_projection_loop.dta", replace

forvalues year = `next5yr'(5)`b' {
	
use "`data_path'\kf45_projection_loop.dta", clear
	
gen prYr=`year'
*gen prYr=`next5yr'
gen Yk= zk1
gen h=(`b'-`a')/8

drop zk1
****************************************
* Equations for Generating Coefficients
*Compute k's for first step, k1-k6
****************************************
gen f1=r*Yk*(1-(Yk/k))
gen k1=h*f1

gen f2=Yk+((1/4)*k1)
gen k2=h*f2*r*(1-(f2/k)) 

gen f3=Yk+(3/32*k1)+(9/32*k2)
gen k3=h*f3*r*(1-(f3/k)) 

gen f4=Yk+((1932/2197)*k1)-((7200/2197)*k2)+((7296/2197)*k3)
gen k4=h*f4*r*(1-(f4/k))

gen f5=(Yk+((439/216)*k1))-(8*k2)+((3680/513)*k3)-((845/4104)*k4)
gen k5=h*f5*r*(1-(f5/k))

gen f6=(Yk -((8/27)*k1) +(2*k2) - ((3544/2565)*k3) + ((1859/4104) *k4) -((11/40) * k5))
gen k6=h*f6*r*(1-(f6/k))

gen yk1=Yk+(25/216*k1)+(1408/2565*k3)+(2197/4101*k4)-(1/5*k5)
gen zk1=Yk+(16/135*k1)+(6656/12825*k3)+(28561/56430*k4)-(9/50*k5)+(2/55*k6)

gen R=(1/h)*abs(yk1-zk1)
gen dlt=0.84*(`eps'/R)^(1/4)
gen iter=1
save "`data_path'\kf45_computeddata_np.dta" , replace

keep if R>`eps'
save "`data_path'\kf45_intermediate_np.dta" , replace
clear
use "`data_path'\kf45_intermediate_np.dta" 

if R>`eps'{
		local n = 0
while `n' >= 0{
		clear
		use "`data_path'\kf45_intermediate_np.dta" 
		replace Yk = zk1 
		keep recid name pop91 pop01 pop11 k r1 r2 r h dlt Yk iter prYr
		replace h=h*dlt
		drop dlt
		save "`data_path'\kf45_residuals_np.dta" , replace
		clear
		use "`data_path'\kf45_residuals_np.dta" 
****************************************
* Equations for Generating Coefficients
*Compute k's for first step, k1-k6
****************************************	
		gen f1=r*Yk*(1-(Yk/k))
		gen k1=h*f1

		gen f2=Yk+((1/4)*k1)
		gen k2=h*f2*r*(1-(f2/k)) 

		gen f3=Yk+(3/32*k1)+(9/32*k2)
		gen k3=h*f3*r*(1-(f3/k)) 

		gen f4=Yk+((1932/2197)*k1)-((7200/2197)*k2)+((7296/2197)*k3)
		gen k4=h*f4*r*(1-(f4/k))

		gen f5=(Yk+((439/216)*k1))-(8*k2)+((3680/513)*k3)-((845/4104)*k4)
		gen k5=h*f5*r*(1-(f5/k))

		gen f6=(Yk -((8/27)*k1) +(2*k2) - ((3544/2565)*k3) + ((1859/4104) *k4) -((11/40) * k5))
		gen k6=h*f6*r*(1-(f6/k))

		gen yk1=Yk+(25/216*k1)+(1408/2565*k3)+(2197/4101*k4)-(1/5*k5)
		gen zk1=Yk+(16/135*k1)+(6656/12825*k3)+(28561/56430*k4)-(9/50*k5)+(2/55*k6)
		
		gen R=(1/h)*abs(yk1-zk1)
		gen dlt=0.84*(`eps'/R)^(1/4)
		replace iter=iter+1
		save "`data_path'\kf45_residuals_np.dta" , replace

		use "`data_path'\kf45_computeddata_np.dta", clear
		append using "`data_path'\kf45_residuals_np.dta"
		save "`data_path'\kf45_computeddata_np.dta" , replace
		clear
		
		use "`data_path'\kf45_residuals_np.dta"
		****************************************
		*Further extending loop
		****************************************
		keep if R>`eps'
		save "`data_path'\kf45_intermediate_np.dta" , replace
		clear
		use "`data_path'\kf45_intermediate_np.dta" 
		count if !missing(recid)   // condition for existing records
		local n = r(N)

		display "Number of records found: `n'"
		
		
		if `n' == 0 {
			sleep 1000   // wait 1 seconds before checking again
			
			 continue, break
		}
	}	

display "Records now exist in the dataset."
	
}
use "`data_path'\kf45_computeddata_np.dta", clear
sort recid iter
order iter, after(recid)

****************************************************
*Final output saved in file
****************************************************
save "`data_path'\kf45_finaloutput_np.dta", replace

*****************************************************
* working files saved ,append, for additional parameters
******************************************************
*table name, statistic(max iter)

bysort recid: egen maxiter = max(iter)
keep if iter == maxiter
drop maxiter
save "`data_path'\kf45_projection2.dta", replace

use "`data_path'\kf45_projection2.dta",clear
append using "`data_path'\kf45_collection_projections.dta"
		save "`data_path'\kf45_collection_projections.dta" , replace
		clear
use "`data_path'\kf45_projection2.dta", clear	
keep recid name pop91 pop01 pop11 k r1 r2 r  zk1
save "`data_path'\kf45_projection_loop.dta", replace
		
    display `year'
}

use "`data_path'\kf45_collection_projections.dta", clear 
sort recid prYr
save "`data_path'\kf45_collection_projections.dta" , replace


*****************************************************
* Generate Excel files of the output.
******************************************************
export excel using "`data_path'\kf45_output_collection_projections.xlsx", firstrow(variables) replace

keep recid  name prYr zk1 
reshape wide zk1, i(recid  name)  j(prYr)
browse
export excel using "`data_path'\kf45_output_pop_projections9631.xlsx", firstrow(variables) replace
