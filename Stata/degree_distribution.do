**** Project title: The structure of risk-sharing networks
**** Authors: H. Henderson and A. Alam
**** Last updated: 1/1/2021
**** Purpose: Creates results for Section 3.2 of the paper

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
gen desire = willingness_link2 // In-degree for desire-to-link model
order hh1 hh2 overreport underreport desire

* Create dataset that lists the degree for each hh
by hh1, sort: egen degree_overreport = total(overreport)
by hh1, sort: egen degree_underreport = total(underreport)
by hh1, sort: egen degree_desire = total(desire)
collapse degree_overreport degree_underreport degree_desire, by(hh1) 

* Generate frequency distribution for each model
by degree_overreport, sort: egen freq_overreport = count(degree_overreport)
replace freq_overreport = freq_overreport/119 
by degree_underreport, sort: egen freq_underreport = count(degree_underreport)
replace freq_underreport = freq_underreport/119 
by degree_desire, sort: egen freq_desire = count(degree_desire)
replace freq_desire = freq_desire/119 

**** Create graphs for underreporting model
preserve 
collapse freq_underreport, by(degree_underreport)
gen cum_freq=sum(freq_underreport) // Cumulative frequency
gen comp_cum_freq = 1 - cum_freq // Complementary cumulative frequency
gen log_comp_cum_freq = ln(comp_cum_freq + 0.01) // Log of CCF

* Calculate r
local r "1"  // Starting value
local m "4.12" // Average degree for underreporting model (half of avg. degree)
local d "10"   // Arbitrary starting point
local tol "0.0001"
while `d' > `tol' { 
    gen rhs = ln(degree + `r'*`m')
    quietly reg log_comp_cum_freq rhs 
    local d = abs(`r' - (- _b[rhs] - 1)) // Difference between r values
    local r = - _b[rhs] - 1  // Update r
    disp `r' 
    drop rhs
    }
 
* Generate variables for plots
gen log_degree = ln(degree_underreport)
gen log_pois_comp_cum_freq = ln(1 - poisson(2*`m', degree_underreport))
gen rhs = ln(degree_underreport + `r'*`m')  // Uses final r
quietly reg log_comp_cum_freq rhs 
predict log_pred_comp_cum_freq
    
* Compare actual vs. Poisson
twoway (line log_comp_cum_freq log_degree, lpattern(dash_dot) lcolor(black) lwidth(medthick)) ///
(line log_pois_comp_cum_freq log_degree, lcolor(black) lwidth(medthick)), ///
legend(label(1 "Nyakatoke") label (2 "Random")) xtitle("Log Degree") ///
ytitle("Log CCDF") title("{bf: (a) Underreporting}")  ///
scheme(s1mono) name(under1, replace) nodraw  

* Compare actual vs. Jackson-Rogers 
twoway (line log_comp_cum_freq log_degree, lpattern(dash_dot) lcolor(black) lwidth(medthick)) ///
(line log_pred_comp_cum_freq log_degree, lcolor(black) lwidth(medthick)), ///
legend(label(1 "Nyakatoke") label (2 "Jackson-Rogers")) xtitle("Log Degree") ///
ytitle("Log CCDF") title("{bf: (a) Underreporting}")  ///
scheme(s1mono) name(under2, replace) nodraw 
restore

**** Create graphs for overreporting model
preserve 
collapse freq_overreport, by(degree_overreport)
gen cum_freq=sum(freq_overreport) // Cumulative frequency
gen comp_cum_freq = 1 - cum_freq // Complementary cumulative frequency
gen log_comp_cum_freq = ln(comp_cum_freq + 0.01) // Log of CCF

* Calculate r
local r "1"  // Starting value
local m "1.175" // Average degree for overreporting model (half of avg. degree)
local d "10"   // Arbitrary starting point
local tol "0.0001"
while `d' > `tol' { 
    gen rhs = ln(degree + `r'*`m')
    quietly reg log_comp_cum_freq rhs 
    local d = abs(`r' - (- _b[rhs] - 1)) // Difference between r values
    local r = - _b[rhs] - 1  // Update r
    disp `r' 
    drop rhs
    }
    
* Generate variables for plots
gen log_degree = ln(degree_overreport)
gen log_pois_comp_cum_freq = ln(1 - poisson(2*`m', degree_overreport))
gen rhs = ln(degree_overreport + `r'*`m')  // Uses final r
quietly reg log_comp_cum_freq rhs 
predict log_pred_comp_cum_freq 
    
* Compare actual vs. Poisson
twoway (line log_comp_cum_freq log_degree, lpattern(dash_dot) lcolor(black) lwidth(medthick)) ///
(line log_pois_comp_cum_freq log_degree, lcolor(black) lwidth(medthick)), ///
legend(label(1 "Nyakatoke") label (2 "Random")) xtitle("Log Degree") ///
ytitle("Log CCDF") title("{bf: (b) Overreporting}")  ///
scheme(s1mono) name(over1, replace) nodraw 

* Compare actual vs. Jackson-Rogers 
twoway (line log_comp_cum_freq log_degree, lpattern(dash_dot) lcolor(black) lwidth(medthick)) ///
(line log_pred_comp_cum_freq log_degree, lcolor(black) lwidth(medthick)), ///
legend(label(1 "Nyakatoke") label (2 "Jackson-Rogers")) xtitle("Log Degree") ///
ytitle("Log CCDF") title("{bf: (b) Overreporting}")  ///
scheme(s1mono) name(over2, replace) nodraw 
restore

**** Create graphs for desire-to-link model
preserve 
collapse freq_desire, by(degree_desire)
gen cum_freq=sum(freq_desire) // Cumulative frequency
gen comp_cum_freq = 1 - cum_freq // Complementary cumulative frequency
gen log_comp_cum_freq = ln(comp_cum_freq + 0.01) // Log of CCF

* Calculate r
local r "1"  // Starting value
local m "5.29" // Average degree for desire-to-link model (equal to avg. degree)
local d "10"   // Arbitrary starting point
local tol "0.0001"
while `d' > `tol' { 
    gen rhs = ln(degree + `r'*`m')
    quietly reg log_comp_cum_freq rhs 
    local d = abs(`r' - (- _b[rhs] - 1)) // Difference between r values
    local r = - _b[rhs] - 1  // Update r
    disp `r' 
    drop rhs
    }
    
* Generate variables for plots
gen log_degree = ln(degree_desire)
gen log_pois_comp_cum_freq = ln(1 - poisson(`m', degree_desire))
gen rhs = ln(degree_desire + `r'*`m')  // Uses final r
quietly reg log_comp_cum_freq rhs 
predict log_pred_comp_cum_freq
    
* Compare actual vs. Poisson
twoway (line log_comp_cum_freq log_degree, lpattern(dash_dot) lcolor(black) lwidth(medthick)) ///
(line log_pois_comp_cum_freq log_degree, lcolor(black) lwidth(medthick)), ///
legend(label(1 "Nyakatoke") label (2 "Random")) xtitle("Log Degree") ///
ytitle("Log CCDF") title("{bf: (c) Desire-to-Link}")  ///
scheme(s1mono) name(desire1, replace) nodraw 

* Compare actual vs. Jackson-Rogers 
twoway (line log_comp_cum_freq log_degree, lpattern(dash_dot) lcolor(black) lwidth(medthick)) ///
(line log_pred_comp_cum_freq log_degree, lcolor(black) lwidth(medthick)), ///
legend(label(1 "Nyakatoke") label (2 "Jackson-Rogers")) xtitle("Log Degree") ///
ytitle("Log CCDF") title("{bf: (c) Desire-to-Link}")  ///
scheme(s1mono) name(desire2, replace) nodraw 
restore

* Combine graphs
* Note: grc1leg is to combine the graphs into one, and create a common legend.
* The aspect ratio is off with the default options. The easiest way to fix this
* is to generate the graph, and then use the graph editor to change the size.
* Here the height is set to 8 inches.
grc1leg under1 over1 desire1, legendfrom(under1) rows(3) scheme(s1mono) 
grc1leg under2 over2 desire2, legendfrom(under2) rows(3) scheme(s1mono) 












