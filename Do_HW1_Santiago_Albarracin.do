*Quantitative Macroeconomics
*Homework 1
*Santiago Albarracin

clear
cd "/Users/santiago/Documents/OneDrive/Documents/QEM Joint Degree/3. Universitat Autònoma de Barcelona/3rd Semester/Quantitative Macroeconomics/Homework/"
use "cps_00001.dta"

*Question 1



*	1.1
*Generate time variable
sort year month
gen time = ym(year, month) 
format time %tm 

*Generate employment situation
gen unemp=0
replace unemp=1 if empstat==21 | empstat==22

gen laborforce=0
replace laborforce=1 if labforce==2

*Generate unemployment
egen tot_unemp = sum(unemp), by(time)
egen pea = sum(laborforce), by(time)

gen unemployment=tot_unemp/pea
gen employment=1-unemployment

*Graph
line employment time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export employment.png, replace



*	1.2
*Generate unemployent by education group
replace educ=. if educ==1

*Less than HighSchool
gen less_hs=0
replace less_hs=1 if educ <= 71
egen less_hs_unemp = sum(unemp)  if less_hs==1,  by(time)
egen less_hs_pea = sum(laborforce) if less_hs==1, by(time)
gen unemployment_less_hs= less_hs_unemp/less_hs_pea
gen employment_less_hs= 1-unemployment_less_hs

*High School
gen hs=0
replace hs=1 if educ >= 72 &  educ<= 110
egen hs_unemp = sum(unemp)  if hs==1,  by(time)
egen hs_pea = sum(laborforce) if hs==1, by(time)
gen unemployment_hs= hs_unemp/hs_pea
gen employment_hs= 1-unemployment_hs

*College
gen coll=0
replace coll=1 if educ >= 111 &  educ<= 122
egen coll_unemp = sum(unemp)  if coll==1,  by(time)
egen coll_pea = sum(laborforce) if coll==1, by(time)
gen unemployment_coll= coll_unemp/coll_pea
gen employment_coll= 1-unemployment_coll

*More than college
gen more_coll=0
replace more_coll=1 if educ>= 123 &  educ<= 125
egen more_coll_unemp = sum(unemp)  if more_coll==1,  by(time)
egen more_coll_pea = sum(laborforce) if more_coll==1, by(time)
gen unemployment_more_coll= more_coll_unemp/more_coll_pea 
gen employment_more_coll= 1-unemployment_more_coll

*Graph
line employment_less_hs employment_hs employment_coll employment_more_coll time, legend(size(medsmall))  note("Source: IPSUM, CPS. Own made")
graph export employment_by_edu.png, replace



*	1.3
*Generate unemployent by industry
replace ind=. if ind==0

*Physical work
gen phy_work=0 
replace phy_work=1 if ind<=6390 | ind>=8560 & ind<=9180
egen phy_unemp=sum(unemp) if phy_work==1, by(time)
egen phy_pea=sum(laborforce) if phy_work==1, by(time)
gen unemployment_phy=phy_unemp/phy_pea
gen employment_phy=1-unemployment_phy

*Intellectual work
gen int_work=0
replace int_work=1 if ind>=9190 | ind<=8470 & ind>=6470
egen int_unemp=sum(unemp) if int_work==1, by(time)
egen int_pea=sum(laborforce) if int_work==1, by(time)
gen unemployment_int=int_unemp/int_pea
gen employment_int=1-unemployment_int

line employment_phy employment_int time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export employment_by_industry.png, replace



*	1.4
*Generate unemployent by occupation
drop if occ==0

*Private
gen private=0
replace private=1 if occ<= 3540
egen private_unemp = sum(unemp)  if private==1,  by(time)
egen private_pea = sum(laborforce) if private==1, by(time)
gen unemployment_priv= private_unemp/private_pea
gen employment_priv= 1-unemployment_priv

*Public
gen public=0
replace public=1 if occ>= 3600 & occ<= 5940
egen public_unemp = sum(unemp)  if public==1,  by(time)
egen public_pea = sum(laborforce) if public==1, by(time)
gen unemployment_pub= public_unemp/public_pea
gen employment_pub= 1-unemployment_pub

*Field
gen field=0
replace field=1 if occ>= 6000
egen field_unemp= sum(unemp)  if field==1,  by(time)
egen field_pea = sum(laborforce) if field==1, by(time)
gen unemployment_field= field_unemp/field_pea
gen employment_field= 1-unemployment_field

*Graph
line employment_priv employment_pub employment_field time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export employment_by_occupation.png, replace



*	2.1
*ssc install asgen

replace uhrsworkt=. if uhrsworkt==997
replace uhrsworkt=. if uhrsworkt==999

egen wh=mean(uhrsworkt), by(time)
*Graph
line wh time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export weekly_hours.png, replace



*	2.2
egen less_hs_wh = mean(uhrsworkt)  if less_hs==1,  by(time)
egen hs_wh = mean(uhrsworkt)  if hs==1,  by(time)
egen coll_wh = mean(uhrsworkt)  if coll==1,  by(time)
egen more_coll_wh = mean(uhrsworkt)  if more_coll==1,  by(time)

*Graph
line less_hs_wh hs_wh coll_wh more_coll_wh time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export weekly_hours_by_edu.png, replace



*	2.3
egen phy_wh=mean(uhrsworkt) if phy_work==1, by(time)
egen int_wh=mean(uhrsworkt) if int_work==1, by(time)


*Graph
line phy_wh int_wh time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export weekly_hours_by_industry.png, replace


*	2.4
egen private_wh = mean(uhrsworkt)  if private==1,  by(time)
egen public_wh = mean(uhrsworkt)  if public==1,  by(time)
egen field_wh = mean(uhrsworkt)  if field==1,  by(time)

*Graph
line private_wh public_wh field_wh time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export weekly_hours_by_occupation.png, replace


*	3
*Aggregate hours






*	4.1
replace hourwage=. if hourwage== 999.99

egen wage=mean(hourwage), by(time)

*Graph
line wage time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export wage.png, replace



*	4.2
egen less_hs_wage = mean(hourwage)  if less_hs==1,  by(time)
egen hs_wage = mean(hourwage)  if hs==1,  by(time)
egen coll_wage = mean(hourwage)  if coll==1,  by(time)
egen more_coll_wage = mean(hourwage)  if more_coll==1,  by(time)

*Graph
line less_hs_wage hs_wage coll_wage more_coll_wage time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export wage_by_edu.png, replace



*	4.3
egen phy_wage=mean(hourwage) if phy_work==1, by(time)
egen int_wage=mean(hourwage) if int_work==1, by(time)

*Graph
line phy_wage int_wage time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export wage_by_industry.png, replace



*	4.4
egen private_wage = mean(hourwage)  if private==1,  by(time)
egen public_wage = mean(hourwage)  if public==1,  by(time)
egen field_wage = mean(hourwage)  if field==1,  by(time)

*Graph
line private_wage public_wage field_wage time, legend(size(medsmall)) note("Source: IPSUM, CPS. Own made")
graph export wage_by_occupation.png, replace



