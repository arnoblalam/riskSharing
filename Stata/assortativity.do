**** Project title: The structure of risk-sharing networks
**** Authors: H. Henderson and A. Alam
**** Last updated: 1/1/2021
**** Purpose: Creates results for Section 3.3 of the paper

* Open Tanzania data
clear all
cd "/Users/hendersonhl/Documents/Articles/Risk-Sharing-Networks/Code"
insheet using "Nyakatoke.csv"

* Generate link variables for the three models
sort hh1 hh2
gen overreport = 0
replace overreport=1 if willingness_link1==1 & willingness_link2==1
gen underreport = 0 
replace underreport=1 if willingness_link1==1 | willingness_link2==1
order hh1 hh2 overreport underreport

*** Heat map for underreporting model
preserve
drop if underreport==0
by hh1, sort: egen d1_underreport = total(underreport) 
by hh2, sort: egen d2_underreport = total(underreport)
replace d1_underreport = d1_underreport - 1 // Excess degrees
replace d2_underreport = d2_underreport - 1
by d1_underreport d2_underreport, sort: egen jp = count(d1_underreport) 
replace jp = jp/_N  // Divide by the number of links
pwcorr d1_underreport d2_underreport, sig // Degree correlation coefficient
collapse jp, by(d1_underreport d2_underreport)
tempfile temp_under
save "`temp_under'"
clear

* Missing to zero
set obs 1024  // Equals max. degree squared
gen d1_underreport=.
gen d2_underreport=.
local k=1
forvalues i=1/32 {
  forvalues j=1/32 {
    qui replace d1_underreport = `i' - 1 if _n==`k'
    qui replace d2_underreport = `j' - 1 if _n==`k'
	local k = `k' + 1
  }
}
merge 1:1 d1_underreport d2_underreport using "`temp_under'" // Merge data in
drop _merge
replace jp=0 if jp==.  // Missing to zero
sum jp

* Plot
twoway contour jp d1_underreport d2_underreport, heatmap ccuts(0(0.002)0.016) ///
xtitle("") ytitle("") ztitle("") title("{bf: (a) Underreporting}") plotr(m(zero)) ///
xlabels(#6) ylabels(#6) scheme(s1mono) zlabels(0(0.002)0.016, format(%04.3f)) ///
name(assort1, replace) nodraw
restore

*** Heat map for overreporting model
preserve
drop if overreport==0
by hh1, sort: egen d1_overreport = total(overreport) 
by hh2, sort: egen d2_overreport = total(overreport)
replace d1_overreport = d1_overreport - 1 // Excess degrees
replace d2_overreport = d2_overreport - 1
by d1_overreport d2_overreport, sort: egen jp = count(d1_overreport) 
replace jp = jp/_N  // Divide by the number of links
pwcorr d1_overreport d2_overreport, sig // Degree correlation coefficient
collapse jp, by(d1_overreport d2_overreport)
tempfile temp_over
save "`temp_over'"
clear

* Missing to zero
set obs 100  // Equals max. degree squared
gen d1_overreport=.
gen d2_overreport=.
local k=1
forvalues i=1/10 {
  forvalues j=1/10 {
    qui replace d1_overreport = `i' - 1 if _n==`k'
    qui replace d2_overreport = `j' - 1 if _n==`k'
	local k = `k' + 1
  }
}
merge 1:1 d1_overreport d2_overreport using "`temp_over'" // Merge data in
drop _merge
replace jp=0 if jp==.  // Missing to zero
sum jp

* Plot
twoway contour jp d1_overreport d2_overreport, heatmap scheme(s1mono) ///
plotr(m(zero)) xtitle("") ytitle("") ztitle("") title("{bf: (b) Overreporting}") ///
xlabels(#9) ylabels(#9) ccuts(0(0.01)0.07) zlabels(0(0.01)0.07, format(%04.3f)) ///
name(assort2, replace) nodraw
restore

**** Heat map for desire-to-link model
* Note: A couple households do not report links. Under the desire-to-link
* assumption, households can thus have a negative excess out-degrees.
* Because of this, we simply use out-degree here.
preserve // First calculate in-degree correlation coefficient
by hh1, sort: egen d1_desire = total(willingness_link2) // In-degrees
by hh2, sort: egen d2_desire = total(willingness_link1)
drop if willingness_link1==0
pwcorr d1_desire d2_desire, sig // In-degree correlation coefficient
restore
preserve
by hh1, sort: egen d1_desire = total(willingness_link1) // Out-degrees
by hh2, sort: egen d2_desire = total(willingness_link2)
drop if willingness_link1==0
by d1_desire d2_desire, sort: egen jp = count(d1_desire) 
replace jp = jp/_N  // Divide by the number of links
pwcorr d1_desire d2_desire, sig // Out-degree correlation coefficient
collapse jp, by(d1_desire d2_desire)
tempfile temp_desire
save "`temp_desire'"
clear

* Missing to zero
set obs 441  // Equals max. degree squared
gen d1_desire=.
gen d2_desire=.
local k=1
forvalues i=1/21 {
  forvalues j=1/21 {
    qui replace d1_desire = `i' - 1 if _n==`k'
    qui replace d2_desire = `j' - 1 if _n==`k'
	local k = `k' + 1
  }
}
merge 1:1 d1_desire d2_desire using "`temp_desire'" // Merge data in
drop _merge
replace jp=0 if jp==.  // Missing to zero
sum jp

* Plot
twoway contour jp d1_desire d2_desire, heatmap scheme(s1mono) ///
plotr(m(zero)) xtitle("") ytitle("") ztitle("") title("{bf: (c) Desire-to-Link}") ///
xlabels(0 5 10  15 20) ylabels(0 5 10 15 20) ccuts(0(0.005)0.04) ///
zlabels(0(0.005)0.04, format(%04.3f)) name(assort3, replace) nodraw
restore

* Combine graphs
* Note: The aspect ratio is off with the default options. The easiest way
* to fix this is to generate the graph, and then use the graph editor to
* change the size. Here the height is set to 8 inches.
graph combine assort1 assort2 assort3, rows(3) scheme(s1mono)













