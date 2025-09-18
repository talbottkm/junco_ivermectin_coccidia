# junco_ivermectin_coccidia
Code and data for project investigating sex differences in the impact of ivermectin on coccidia shedding in dark-eyed juncos

## ArielREU_kmt20250918.R
Code for analyses and visuals

## master_juncococcidia.xlsx
Data file for analysis; two tabs labeled 'long' and 'longchange'

Column headers for 'longchange'  
band: individual ID  
sex: 1=female, 0=male  
treatment: 1=ivermectin, 0=control  
group: sex + treatment  
dcount: change in oocyst count   
change: timeframe of change (note "pd" means "post-dose")  
malaria: treatment (if any) during an experimental project that occurred previously  
mtreat: malaria treatment; 1=Plasmodium inoculation, 0=no Plasmodium inoculation   
treat2: 'group' + 'mtreat'   

Column headers for 'long'  
band: individual ID  
sex: 1=female, 0=male  
treatment: 1=ivermectin, 0=control  
group: sex + treatment  
mass: junco mass in grams  
count_addone: mean number of oocysts plus a value of one  
sample: point during experiment when sample was collected (note 'base2' is baseline)  
count: raw mean oocyst count  
