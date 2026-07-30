## ArielREU_kmt20250918.R
Code for analyses and visuals

## coccidia_longdata.csv
Data file for analysis in long form

## coccidia_widedata.csv
Data file for analysis in wide form

## Column headers for longdata  
band: individual ID  
sex: 1=female, 0=male  
treatment: 1=ivermectin, 0=control  
count: coccidia oocysts per gram feces   
mass: junco mass (grams)
plasmodium_exp: role in Plasmodium experiment prior to Ivermectin study (none = bird was not involved; control = bird inoculated with uninfected blood; exp = bird was inoculated with Plasmodium infected blood)  
plasmodium_treat: whether bird received a Plasmodium inoculation (0 = no, 1 = yes)

## Column headers for widedata  
band: individual ID  
sex: 1=female, 0=male  
treatment: 1=ivermectin, 0=control  
count: coccidia oocysts per gram feces   
mass: junco mass (grams)
plasmodium_exp: role in Plasmodium experiment prior to Ivermectin study (none = bird was not involved; control = bird inoculated with uninfected blood; exp = bird was inoculated with Plasmodium infected blood)  
plasmodium_treat: whether bird received a Plasmodium inoculation (0 = no, 1 = yes)
mass_base: junco mass at baseline (grams)
mass_postdose: junco mass following final ivermectin dose (grams)
count_base: coccidia oocyst count per grams feces at baseline
count_postdose: coccidia oocyst count per grams feces following final ivermectin dose (grams)
dcount: change in coccidia oocyst count from baseline to postdose (dcount = count_postdose - count_base)
