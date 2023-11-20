clear

use ".\datav4.0\Data\2. Demographics and Outcomes\individual_characteristics.dta"
drop if resp_status==2 | resp_status==3
duplicates drop hhid, force
save ".\datav4.0\Data\2. Demographics and Outcomes\individual_characteristics_headonly.dta", replace

use ".\datav4.0\Data\2. Demographics and Outcomes\household_characteristics.dta"

g caste=1 if castesubcaste=="SCHEDULE CASTE"
replace caste=2 if castesubcaste=="SCHEDULE TRIBE"
replace caste=3 if castesubcaste=="OBC"
replace caste=4 if castesubcaste=="GENERAL"
replace caste=5 if castesubcaste=="MINORITY"

merge 1:1 hhid using ".\datav4.0\Data\2. Demographics and Outcomes\individual_characteristics_headonly.dta", keepusing(caste) update replace
replace caste=5 if caste==-999

drop _merge
rename hohreligion religion
merge 1:1 hhid using ".\datav4.0\Data\2. Demographics and Outcomes\individual_characteristics_headonly.dta", keepusing(religion) update replace

drop _merge castesubcaste
drop hhSurveyed adjmatrix_key


// Convert categorical variables to binary variables


// Caste. Base: OBC
g caste_gen = 1 if caste == 4
replace caste_gen=0 if missing(caste_gen)

g caste_min = 1 if caste == 5
replace caste_min=0 if missing(caste_min)

g caste_sc = 1 if caste == 1
replace caste_sc=0 if missing(caste_sc)

g caste_st = 1 if caste == 2
replace caste_st=0 if missing(caste_st)

g caste_obc = 1 if caste == 3
replace caste_obc=0 if missing(caste_obc)



// Religion. Base: Hinduism
g rel_isl = 1 if religion == 2
replace rel_isl = 0 if missing(rel_isl)

g rel_hin = 1 if religion == 1
replace rel_hin = 0 if missing(rel_hin)

drop if religion == 3 // Too little sample for christians
drop religion

// Latrine. Base: anything other than no latrine
g lat_none = 1 if latrine == 3
replace lat_none = 0 if missing(lat_none)

drop latrine

// Roof type. Base: everything except tile

drop rooftype1 rooftype3 rooftype4 rooftype5 rooftypeoth

// Rent status. Base: owned
g notowned = 1 if ownrent != 1
replace notowned = 0 if missing(notowned)
drop ownrent

// Electricity. Base: everything except no electricity

g noelectricity = 1 if electricity != 1
replace noelectricity = 0 if missing(noelectricity)
drop electricity

// Summary statistics 1
summarize

// Drop unused villages

drop if village == 5 | village == 7 | village == 8 | village == 10 | village == 11 | village == 13 | village == 14 | village == 16 | village == 17 | village == 18 | village == 22 | village == 26 | village == 27 | village == 28 | village == 30 | village == 34 | village == 35 | village == 37 | village == 38 | village == 40 | village == 41 | village == 44 | village == 49 | village == 53 | village == 54 | village == 56 | village == 58 | village == 61 | village == 63 | village == 66 | village == 69 | village == 74 | village == 76 | village == 77

summarize

// Drop the villages without caste information
// Summary Statistics 3

drop if missing(caste)
summarize

drop caste
drop rel_hin
order caste_obc, last // the obc variable will be used only for a network specification

export delimited using ".\processed\hh_cov_cleaned.csv", nolabel replace