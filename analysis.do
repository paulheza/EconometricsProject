********************************************************************************
* SILVER PRICES & THE ENERGY TRANSITION
* Full Analysis Do-File
********************************************************************************

clear all
set more off
cd "/Users/admin/Documents/College/Junior Sophister/SEM2/Econometrics/Project/Data"

import delimited "master_dataset.csv", clear varnames(1)

********************************************************************************
* PHASE 1: DATA PREPARATION
********************************************************************************

*Fix date and set time series
gen t = date(date, "YMD")
gen t2 = mofd(t)
format t2 %tm
tsset t2
drop t
rename t2 t

*Verify range
sum t

*Check for gaps
xtset t

*Generate interaction terms
gen indprod_post = dlog_indpro * post2016
gen solar_post   = d_solar_mw  * post2016

*Install estout
ssc install estout, replace

*Summary statistics
estpost summarize dlog_silver dlog_gold dlog_indpro ///
    d_solar_mw d_treasury d_dollar post2016 covid

esttab using "01_summary_stats.rtf", ///
    cells("mean(fmt(4)) sd(fmt(4)) min(fmt(4)) max(fmt(4)) count(fmt(0))") ///
    title("Table 1: Summary Statistics") replace

save "master_dataset.dta", replace

********************************************************************************
* PHASE 2: UNIT ROOT TESTS
********************************************************************************

tsset t

*Log levels - expect to FAIL to reject
dfuller log_silver,   lags(12) trend
dfuller log_gold,     lags(12) trend
dfuller log_indpro,   lags(12) trend
dfuller log_copper,   lags(12) trend
dfuller treasury_10y, lags(12) trend
dfuller dollar_index, lags(12) trend
dfuller solar_mw,     lags(12) trend

*First differences - expect to REJECT
dfuller dlog_silver,  lags(12)
dfuller dlog_gold,    lags(12)
dfuller dlog_indpro,  lags(12)
dfuller dlog_copper,  lags(12)
dfuller d_treasury,   lags(12)
dfuller d_dollar,     lags(12)
dfuller d_solar_mw,   lags(12)

*Phillips-Perron robustness
pperron log_silver,   lags(4) trend
pperron dlog_silver,  lags(4)
pperron log_indpro,   lags(4) trend
pperron dlog_indpro,  lags(4)
pperron solar_mw,     lags(4) trend
pperron d_solar_mw,   lags(4)

*Additional tests for problematic variables
pperron dlog_gold, lags(4)

*Solar stationarity checks
dfuller d_solar_mw, lags(6)

*Generate second difference of solar (I(2) confirmed)
gen d2_solar_mw = D.d_solar_mw
dfuller d2_solar_mw, lags(12)

********************************************************************************
* PHASE 3: STRUCTURAL BREAK TESTS
********************************************************************************

tsset t

*--- Step 3.1: Verify post2016 dummy ---
tab post2016
list t post2016 if t == monthly("2015m12","tm") | t == monthly("2016m1","tm")

*--- Generate interactions for Phase 3 ---
gen d2_solar_post    = d2_solar_mw * post2016
gen d2_solar_mw_lag6 = L6.d2_solar_mw
gen d2_solar_post_l6 = d2_solar_mw_lag6 * post2016
gen gold_post        = dlog_gold  * post2016
gen int_post         = d_treasury * post2016
gen dollar_post      = d_dollar   * post2016

*--- Step 3.2: Chow Test ---
reg dlog_silver dlog_indpro d2_solar_mw dlog_gold d_treasury d_dollar ///
    post2016 indprod_post d2_solar_post gold_post int_post dollar_post, ///
    robust

*Joint significance test = Chow Test
testparm post2016 indprod_post d2_solar_post gold_post int_post dollar_post

*--- Step 3.3: Save before rolling regressions ---
save "master_dataset.dta", replace

*--- Step 3.4: Rolling Window Regressions ---


rolling _b _se, window(36) saving("rolling_indpro.dta", replace): ///
    reg dlog_silver dlog_indpro d2_solar_mw dlog_gold d_treasury d_dollar

*Plot industrial production coefficient
use "rolling_indpro.dta", clear
tsset end

tsline _b_dlog_indpro, ///
    title("Rolling 36-Month: IndProd Coefficient on Silver") ///
    ytitle("Coefficient") xtitle("Window End Date") ///
    yline(0, lcolor(red) lpattern(dash)) ///
    xline(672, lcolor(blue) lpattern(dash)) ///
    note("Blue dashed = Jan 2016. Red dashed = zero.") ///
    scheme(s2color)

graph export "03_rolling_indpro.png", replace

*Plot solar coefficient
tsline _b_d2_solar_mw, ///
    title("Rolling 36-Month: Solar Coefficient on Silver") ///
    ytitle("Coefficient") xtitle("Window End Date") ///
    yline(0, lcolor(red) lpattern(dash)) ///
    xline(672, lcolor(blue) lpattern(dash)) ///
    note("Blue dashed = Jan 2016. Red dashed = zero.") ///
    scheme(s2color)

graph export "03_rolling_solar.png", replace

*Reload main dataset
use "master_dataset.dta", clear
tsset t

********************************************************************************
* PHASE 4: MAIN REGRESSION MODELS
********************************************************************************

tsset t

*--- Model 1: Post-2016 Level Shift ---
reg dlog_silver post2016 dlog_gold d_treasury d_dollar, robust

estimates store model1

*--- Model 2: Industrial Production Interaction (H1 Test) ---
reg dlog_silver dlog_indpro post2016 indprod_post ///
    dlog_gold d_treasury d_dollar, robust

estimates store model2

*Test H1: is the interaction term significant?
test indprod_post

*--- Model 3: Solar-Specific Channel (H2 Test) ---
reg dlog_silver dlog_indpro d2_solar_mw post2016 ///
    indprod_post d2_solar_post ///
    dlog_gold d_treasury d_dollar, robust

estimates store model3

*Test H2: solar interaction significant?
test d2_solar_post

*Has indprod_post shrunk vs Model 2?
test indprod_post

*Joint test of both interactions
test indprod_post d2_solar_post

*--- Export all three models ---
esttab model1 model2 model3 using "04_main_results.rtf", ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    b(4) se(4) ///
    title("Table 2: Main Regression Results") ///
    mtitles("Model 1: Level Shift" "Model 2: H1 Test" "Model 3: H2 Test") ///
    scalars("r2 R-squared" "N Observations") ///
    addnote("Robust standard errors in parentheses. * p<0.10 ** p<0.05 *** p<0.01") ///
    replace

********************************************************************************
* PHASE 5: OLS DIAGNOSTICS
********************************************************************************

tsset t

*--- Re-run Model 3 without robust SEs for diagnostic tests ---
reg dlog_silver dlog_indpro d2_solar_mw post2016 ///
    indprod_post d2_solar_post ///
    dlog_gold d_treasury d_dollar

*--- Step 5.1: Serial Correlation (Breusch-Godfrey) ---
estat bgodfrey, lags(1 2 3 6 12)

*--- Step 5.2: Heteroskedasticity (simplified) ---
estat hettest

*--- Step 5.3: ARCH Effects ---
predict resid3, residuals
estat archlm, lags(1 2 3 6)

*--- Step 5.4: Multicollinearity (VIF) ---
vif

*--- Step 5.5: Normality of residuals ---
sktest resid3
histogram resid3, normal ///
    title("Model 3 Residuals") ///
    scheme(s2color)
graph export "05_residuals_hist.png", replace

********************************************************************************
* PHASE 6: ROBUSTNESS CHECKS
********************************************************************************

tsset t

*--- Step 6.1: Granger Causality ---
var dlog_silver d2_solar_mw, lags(1/6)
vargranger

*--- Step 6.2: 6-Month Lag Robustness ---
reg dlog_silver dlog_indpro d2_solar_mw_lag6 post2016 ///
    indprod_post d2_solar_post_l6 ///
    dlog_gold d_treasury d_dollar, robust

estimates store model3_lag6

*--- Step 6.3: Copper Falsification Test ---
reg dlog_copper dlog_indpro d2_solar_mw post2016 ///
    indprod_post d2_solar_post ///
    dlog_gold d_treasury d_dollar, robust

estimates store copper_test

*Test solar interaction for copper
test d2_solar_post

*--- Step 6.4: COVID Exclusion ---
reg dlog_silver dlog_indpro d2_solar_mw post2016 ///
    indprod_post d2_solar_post ///
    dlog_gold d_treasury d_dollar ///
    if covid == 0, robust

estimates store model3_nocovid

*--- Step 6.5 Fixed: Alternative Break Dates ---

*Regenerate break date variables
gen post2015       = (t >= monthly("2015m1","tm"))
gen post2017       = (t >= monthly("2017m1","tm"))

gen indprod_post15 = dlog_indpro * post2015
gen solar_post15   = d2_solar_mw * post2015

gen indprod_post17 = dlog_indpro * post2017
gen solar_post17   = d2_solar_mw * post2017

*Test 2015 break
reg dlog_silver dlog_indpro d2_solar_mw post2015 ///
    indprod_post15 solar_post15 ///
    dlog_gold d_treasury d_dollar, robust
estimates store break2015

*Test 2016 break (comparable specification)
reg dlog_silver dlog_indpro d2_solar_mw post2016 ///
    indprod_post d2_solar_post ///
    dlog_gold d_treasury d_dollar, robust
estimates store break2016

*Test 2017 break
reg dlog_silver dlog_indpro d2_solar_mw post2017 ///
    indprod_post17 solar_post17 ///
    dlog_gold d_treasury d_dollar, robust
estimates store break2017

*Compare AIC/BIC across all three break dates
estimates stats break2015 break2016 break2017

*--- Export Robustness Table ---
esttab model3 model3_lag6 copper_test model3_nocovid using "06_robustness.rtf", ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    b(4) se(4) ///
    title("Table 3: Robustness Checks") ///
    mtitles("Baseline" "6M Lag Solar" "Copper Falsif." "No COVID") ///
    scalars("r2 R-squared" "N Observations") ///
    addnote("Robust standard errors in parentheses. Copper test: DV = dlog_copper.") ///
    replace
	
********************************************************************************
* PHASE 7: FINAL OUTPUTS
********************************************************************************

tsset t
use "master_dataset.dta", clear
tsset t

*--- Table 1: Unit Root Summary (manual — export to rtf) ---
*This is a formatted summary of your Phase 2 results
*You will fill this in manually in Word from your Phase 2 output

*--- Table 2: Main Results (already exported as 04_main_results.rtf) ---

*--- Table 3: Robustness (already exported as 06_robustness.rtf) ---

*--- Figure 1: Silver Price and Solar Capacity Over Time ---
twoway (tsline real_silver, yaxis(1) lcolor(navy) lwidth(medium)) ///
       (tsline solar_mw, yaxis(2) lcolor(orange) lwidth(medium) lpattern(dash)), ///
    title("Real Silver Price and Global Solar PV Capacity") ///
    ytitle("Real Silver Price (USD/oz)", axis(1)) ///
    ytitle("Cumulative Solar Capacity (MW)", axis(2)) ///
    xtitle("") ///
    xline(672, lcolor(red) lpattern(dash)) ///
    legend(label(1 "Real Silver Price") label(2 "Solar Capacity (MW)")) ///
    note("Red dashed = January 2016 structural break.") ///
    scheme(s2color)
graph export "07_fig1_silver_solar.png", replace

*--- Figure 2: Residual Plot over Time ---
reg dlog_silver dlog_indpro d2_solar_mw post2016 ///
    indprod_post d2_solar_post ///
    dlog_gold d_treasury d_dollar, robust
predict resid_final, residuals

tsline resid_final, ///
    title("Model 3 Residuals Over Time") ///
    ytitle("Residual") xtitle("") ///
    yline(0, lcolor(red) lpattern(dash)) ///
    xline(672, lcolor(blue) lpattern(dash)) ///
    note("Blue dashed = January 2016. Red dashed = zero.") ///
    scheme(s2color)
graph export "07_fig2_residuals.png", replace

*--- Figure 3: Actual vs Fitted ---
predict fitted_final
twoway (tsline dlog_silver, lcolor(navy) lwidth(medium)) ///
       (tsline fitted_final, lcolor(orange) lpattern(dash) lwidth(medium)), ///
    title("Actual vs Fitted: Model 3") ///
    ytitle("Monthly Log Return") xtitle("") ///
    legend(label(1 "Actual") label(2 "Fitted")) ///
    scheme(s2color)
graph export "07_fig3_actual_fitted.png", replace

*--- Final combined regression table: Models 1-3 with controls ---
esttab model1 model2 model3 using "07_final_table.rtf", ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    b(4) se(4) ///
    title("Main Regression Results") ///
    mtitles("Model 1" "Model 2" "Model 3") ///
    varlabels(dlog_indpro "Industrial Production (Δlog)" ///
              d2_solar_mw "Solar Capacity (Δ²MW)" ///
              post2016 "Post-2016 Dummy" ///
              indprod_post "IndProd × Post-2016" ///
              d2_solar_post "Solar × Post-2016" ///
              dlog_gold "Gold Price (Δlog)" ///
              d_treasury "10Y Treasury (Δ)" ///
              d_dollar "Dollar Index (Δ)" ///
              _cons "Constant") ///
    scalars("r2 R-squared" "N Observations") ///
    addnote("Robust standard errors in parentheses." ///
            "* p<0.10 ** p<0.05 *** p<0.01" ///
            "All variables in first differences. Solar in second differences (I(2)).") ///
    replace

*--- Save final dataset ---
save "master_dataset_final.dta", replace
