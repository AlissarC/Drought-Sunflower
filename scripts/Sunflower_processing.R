
##################### packages ##############################
install.packages("LeafArea")
install.packages("rJava")
###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(LeafArea)
library(rJava)


###############################################################################
## Load data : leaves area and leaves disks 
###############################################################################
leaves_area <- read.csv("../input/leaves_area.csv")
leaves_disk <- read.csv("../input/leaves_disk.csv")
nb_disk_plant_time <- read.csv("../input/nb_disk_plant_time.csv")

####################### open files: Join to the file containing the number and biomass of disks for each plant, trt and time: nb_disk_plant_time ###################### 
## In this file the biomass is for the dried disks  after removing the fresh ones for chlorophyll measurements ######
## Except time 1: there was no chlorophyll measurements, the biomass is for all disks 

nb_disk_plant_time$plant <- as.character(nb_disk_plant_time$plant)
leaves_disk$plant <- as.character(leaves_disk$plant)

leaves_disk_area_chl = left_join(nb_disk_plant_time,leaves_disk, by = c("ID_trt", "n_trt", "drought", "plant", "time" )) 
names(leaves_disk_area_chl)
leaves_disk_area_chl$total.leaf.area

#### approximation to deduce the total.leaf area for time 1 since the disks were not scanned by error
# Calculate the mean of other times when disks == 10
subset_data <- subset(leaves_disk_area_chl, nb_disks_scanned == 10)
nrow(subset_data)#121
mean_value_leaf_area <- mean(subset_data$total.leaf.area, na.rm = TRUE) # 6.265537

######## calculate the surface of leaves disks at t1 depending on their numbers based on the mean leaf area value for 10 disks
leaves_disk_area_chl$total.leaf.area[leaves_disk_area_chl$time == "t1"] <- 
  mean_value_leaf_area *leaves_disk_area_chl$nb_disks_scanned[leaves_disk_area_chl$time == "t1"]/10 

#### four disks used to chl extraction: calculate the surface of 4 disks
leaves_disk_area_chl$chl_disk_leaf_area =(leaves_disk_area_chl$total.leaf.area*leaves_disk_area_chl$nb_disks_chlorophyll)/leaves_disk_area_chl$nb_disks_scanned 
names(leaves_disk_area_chl)
nrow(leaves_disk_area_chl)

#### Calculate SLA ############
leaves_disk_area_chl$nb_disk_SLA = leaves_disk_area_chl$nb_disks_scanned - leaves_disk_area_chl$nb_disks_chlorophyll 
leaves_disk_area_chl$leaf_area_for_SLA = (leaves_disk_area_chl$nb_disk_SLA*leaves_disk_area_chl$total.leaf.area)/leaves_disk_area_chl$nb_disks_scanned
leaves_disk_area_chl$LMA = (leaves_disk_area_chl$leaves_disk_mg_noChl/leaves_disk_area_chl$leaf_area_for_SLA)*10 ## *10 to convert from mg/cm2 to g/m2
leaves_disk_area_chl$SLA = 1/leaves_disk_area_chl$LMA #m2/g
names(leaves_disk_area_chl)
nrow(leaves_disk_area_chl)

max(leaves_disk_area_chl$LMA) # 101.23
min(leaves_disk_area_chl$LMA) # 29.72

# open file with chlorophyll measurements ###################
Chl_data <- read.csv("../input/Chlorophyll_all.csv")
names(Chl_data)

Chl_data <- separate(Chl_data, col = Name, into=c("ID_trt", "n_trt", "drought", "plant", "time"), sep= "_", remove=FALSE)
Chl_data$plant <- gsub("P", "", Chl_data$plant) ### select the number of the plant
Chl_data$plant <- as.character(Chl_data$plant)
view(Chl_data) # Verification 
names(Chl_data)

## calculate the mean of technical repetitions 
Chl_data_mean <- Chl_data %>% 
  group_by(ID_trt, n_trt, drought, plant, time) %>% 
  summarize(Abs_649_mean = mean (Blank.chlorophyll_absorbance.649, na.rm=TRUE), Abs_665_mean = mean (Blank.chlorophyll_absorbance.665, na.rm=TRUE), 
            Abs_649_sd= sd (Blank.chlorophyll_absorbance.649, na.rm=TRUE),  Abs_665_sd= sd (Blank.chlorophyll_absorbance.665, na.rm=TRUE))

hist(Chl_data_mean$Abs_649_mean)
hist(Chl_data_mean$Abs_665_mean)

hist(Chl_data_mean$Abs_649_sd)
hist(Chl_data_mean$Abs_665_sd)

view(Chl_data_mean)
## join chlorophyll to leaves disk surface area 
chl_area = left_join(leaves_disk_area_chl,Chl_data_mean, by = c("ID_trt", "n_trt", "drought", "plant", "time" )) 
chl_area$Abs_649_mean
names(chl_area)
nrow(chl_area)
view(chl_area)

########## Estimate chl in mmole/g #######

### Chl A and chl B ug per ml 
chl_area$chlA_ugml = 12.47 * chl_area$Abs_665_mean - 3.62 * chl_area$Abs_649_mean
chl_area$chlB_ugml = 25.06 * chl_area$Abs_649_mean - 6.5 * chl_area$Abs_665_mean

### Chl A and chl B g per ml 
chl_area$chlA_gml = chl_area$chlA_ugml / 1000000
chl_area$chlB_gml = chl_area$chlB_ugml / 1000000

### Chl A and chl B g per 10 ml (the volume of the extractant)
chl_area$chlA_g = chl_area$chlA_gml * 10 # extracted in 10mL DMSO
chl_area$chlB_g = chl_area$chlB_gml * 10 # extracted in 10mL DMSO

### Chl A and chl B g per m2 : chl_disk_leaf_area: leaf area for the 4 disks used to extract cholrophyll
chl_area$chlA_gm2 = chl_area$chlA_g / (chl_area$chl_disk_leaf_area / 10000)
chl_area$chlB_gm2 = chl_area$chlB_g / (chl_area$chl_disk_leaf_area / 10000)

### Chl A and chl B mmol per m2
chl_area$chlA_mmolm2 = chl_area$chlA_gm2 / 893.51 * 1000
chl_area$chlB_mmolm2 = chl_area$chlB_gm2 / 907.47 * 1000

### Chl A and chl B mmol per g
chl_area$chlA_mmolg = chl_area$chlA_mmolm2/chl_area$LMA
chl_area$chlB_mmolg = chl_area$chlB_mmolm2/chl_area$LMA

#### total chlorophyll #chl_LMA_disk$chl_mmolg = chl_LMA_disk$chlA_mmolg  + chl_LMA_disk$chlB_mmolg
chl_area$chl_mmolm2 = chl_area$chlA_mmolm2  + chl_area$chlB_mmolm2

#### total chlorophyll #chl_LMA_disk$chl_mmolg = chl_LMA_disk$chlA_mmolg  + chl_LMA_disk$chlB_mmolg
chl_area$chl_gm2 = chl_area$chlA_gm2  + chl_area$chlB_gm2
chl_area$chlmass= chl_area$chl_gm2/chl_area$LMA
chl_area$chl_mmolg = chl_area$chlA_mmolg  + chl_area$chlB_mmolg

## ratio chla/chlb
chl_area$ratio_chlA_B = chl_area$chlA_mmolm2 / chl_area$chlB_mmolm2
nrow(chl_area)

### Isotopes, leaf nitrogen, A-Ci and biomass data ################
A_Ci_N_isotopes_biomass <-  read.csv("../input/A_Ci_N_isotopes_biomass.csv") 
names(A_Ci_N_isotopes_biomass)
nrow(A_Ci_N_isotopes_biomass) #192

table(A_Ci_N_isotopes_biomass$time)
table(A_Ci_N_isotopes_biomass$block)
table(A_Ci_N_isotopes_biomass$drought)

################ Big Delta calculus #####################################
delta_air = -8
A_Ci_N_isotopes_biomass$big_D13 = (delta_air - A_Ci_N_isotopes_biomass$deltaC13)/(1+A_Ci_N_isotopes_biomass$deltaC13*0.001)
hist(A_Ci_N_isotopes_biomass$big_D13)

######################## Chi calculus ###############################
a = 4.4
b = 28
A_Ci_N_isotopes_biomass$chi= (A_Ci_N_isotopes_biomass$big_D13 -a)/(b-a)
hist(A_Ci_N_isotopes_biomass$chi)
names(A_Ci_N_isotopes_biomass)

################ Merge Iso_N_Aci and chlorophyll  data ################## 
A_Ci_N_isotopes_biomass$plant <- as.character(A_Ci_N_isotopes_biomass$plant)

Iso_N_Aci_chl = left_join(A_Ci_N_isotopes_biomass,chl_area, by = c("ID_trt", "n_trt", "drought", "plant",
                                                                        "time", "date")) 
names(Iso_N_Aci_chl)
nrow(Iso_N_Aci_chl)
View(Iso_N_Aci_chl)


######## calculate Nmass 
Iso_N_Aci_chl$Nmass= Iso_N_Aci_chl$percentage_N/100
hist(Iso_N_Aci_chl$Nmass)

######## calculate Narea 
Iso_N_Aci_chl$Narea= Iso_N_Aci_chl$Nmass * Iso_N_Aci_chl$LMA
hist(Iso_N_Aci_chl$Narea)

### calculate the relative proportion of nitrogen in Rubisco (ρrubisco; gN gN−1) Waring et al., 2023)
Nr = 0.16 # the amount of nitrogen in Rubisco, assumed to be 0.16 gN (g Rubisco)−1
Vcr = 20.5 # is the specific activity of Rubisco, assumed to be 20.5 µmol CO2 (g Rubisco)−1 s−1 at 25 °C
Iso_N_Aci_chl$nrubisco = (Iso_N_Aci_chl$vcmax_tleaf * Nr)/(Iso_N_Aci_chl$Narea*Vcr)
hist(Iso_N_Aci_chl$nrubisco)

### calculate the relative proportion of nitrogen in bioenergetics: ρbioenergetics; gN gN−1 (Waring et al., 2023)
Nb = 0.1240695 # the amount of nitrogen in cytochrome f, assumed to be 0.1240695 g N (µmol cytochrome f)−1
Jmc = 156 # a the capacity of electron transport per cytochrome f, set to 156 μmol electron (μmol cytochrome f)−1 s−1 

Iso_N_Aci_chl$nbioe = (Iso_N_Aci_chl$jmax_tleaf * Nb)/(Iso_N_Aci_chl$Narea*Jmc)
hist(Iso_N_Aci_chl$nbioe)

### calculate the relative proportion of nitrogen in light harvesting (ρlightharvesting; gN gN−1), Waring et al., 2023)
Cb = 2.75 #  the chlorophyll binding of the thylakoid protein complexes, assumed to be 2.75 mmol chlorophyll (g chlorophyll N)−1

Iso_N_Aci_chl$nlightharvesting = (Iso_N_Aci_chl$chl_mmolm2)/(Iso_N_Aci_chl$Narea*Cb)
hist(Iso_N_Aci_chl$nlightharvesting)

############### Total nitrogen allocated to photosynthesis #########
Iso_N_Aci_chl$nphoto = Iso_N_Aci_chl$nrubisco + Iso_N_Aci_chl$nbioe + Iso_N_Aci_chl$nlightharvesting

### Ncw : N allocated to the cell wall
Iso_N_Aci_chl$Ncw= 0.000355*Iso_N_Aci_chl$LMA^1.39

############## N structure ##############
Iso_N_Aci_chl$Nstrucutre = Iso_N_Aci_chl$Ncw/ Iso_N_Aci_chl$Narea
#Iso_N_Aci_chl$Nstrucutre = (10^-2.67) * (Iso_N_Aci_chl$lma ^ 0.99)
names(Iso_N_Aci_chl)

############################# call functions to calculate beta and chi #################################################
install.packages ("R.utils") # to read source function
library (R.utils)

getwd()
sourceDirectory("../scripts/functions", modifiedOnly = FALSE, verbose = TRUE)
################ calculate atmospheric pressure (Pa) from elevation (m)#####################
Iso_N_Aci_chl$z = 976
calc_patm = function(z) {
  
  kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
  kTo = 298.15   # base temperature, K (Prentice, unpublished)
  kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
  kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
  kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
  kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  
  patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
  
  patm
}

Iso_N_Aci_chl$patm<- calc_patm(Iso_N_Aci_chl$z)

################ calculate  nstar (unitless relative viscosity of h2o at temperature relative to 25°C) ####################
calc_nstar = function(temp, z){ # temp in °C and z in m
  
  patm = calc_patm(z)
  
  # viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
  ns      = calc_viscosity_h2o( temp, z )  # Pa s 
  ns25    = calc_viscosity_h2o( 25, z )  # Pa s 
  nstar = ns / ns25                       # (unitless)
  
  nstar
  
}
# Apparently it call automatically calc_viscositt_h2o 
Iso_N_Aci_chl$nstar<- calc_nstar(temp = Iso_N_Aci_chl$tmp, z = Iso_N_Aci_chl$z)

########################## calculate gammastar (Pa) ###############################################
calc_gammastar_pa = function(temp, z) {
  
  patm = calc_patm(z)
  rat = calc_patm(z) / calc_patm(0)
  gammastar25 = 4.332 * rat  # Pa
  Hgm=37830 # J mol-1
  R = 8.314        # J K-1 mol-1
  O2 = 2.09476e5 # ppm
  O2_0 = O2 * 1e-6 * calc_patm(0)
  O2_z = O2 * 1e-6 * calc_patm(z)
  
  temp_k = 273.15+ temp
  
  gStar_pa = gammastar25*exp((Hgm/R)*(1/298.15-1/temp_k))
  
  gStar_pa
  
}
Iso_N_Aci_chl$gammastar_pa <-calc_gammastar_pa(temp = Iso_N_Aci_chl$tmp, z = Iso_N_Aci_chl$z)

######## calcualte the Michaelis-Menton coefficient (Pa) for Rubisco from temperature##############################
calc_km_pa = function(temp, z) {
  
  patm = calc_patm(z) 
  rat = patm / calc_patm(0)
  
  R = 8.314        
  O2 = 2.09476e5      
  Kc25 = 41.03 * rat 
  Ko25 = 28210 * rat 
  Hkc = 79430  
  Hko = 36380 
  
  temp_k = 273.15 + temp
  
  Kc_pa =Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  Ko_pa =Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  O2_pa = O2 * (1e-6) * patm 
  
  Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
  
  Km_pa 
  
}

Iso_N_Aci_chl$Km<-calc_km_pa(temp = Iso_N_Aci_chl$tmp, z = Iso_N_Aci_chl$z)

names(Iso_N_Aci_chl)
##################### Beta calculus ######################################
Iso_N_Aci_chl$D <- Iso_N_Aci_chl$VPD*1000
Iso_N_Aci_chl$ca <- 420
Iso_N_Aci_chl$beta_num <- 1.6 * Iso_N_Aci_chl$nstar * Iso_N_Aci_chl$D * ((Iso_N_Aci_chl$chi-Iso_N_Aci_chl$gammastar_pa/Iso_N_Aci_chl$ca) ^ 2)
Iso_N_Aci_chl$beta_denom <- ((1 - Iso_N_Aci_chl$chi)^2) * (Iso_N_Aci_chl$Km + Iso_N_Aci_chl$gammastar_pa)

Iso_N_Aci_chl$betaleaf <- Iso_N_Aci_chl$beta_num / Iso_N_Aci_chl$beta_denom

summary(Iso_N_Aci_chl$beta_num)
hist(Iso_N_Aci_chl$beta_num)
plot(Iso_N_Aci_chl$beta_num, Iso_N_Aci_chl$chi)
plot(Iso_N_Aci_chl$beta_denom, Iso_N_Aci_chl$chi)
plot(Iso_N_Aci_chl$chi, Iso_N_Aci_chl$betaleaf)
plot(Iso_N_Aci_chl$chi, log(Iso_N_Aci_chl$betaleaf))


##################### Morphology ##################################################
morpho_all_sum_plant <- read.csv("../input/morpho_all_sum_plant.csv")

######### Combine with Iso_N_Aci_chl
morpho_all_sum_plant$plant <- gsub("P", "", morpho_all_sum_plant$plant) ### select the number of the plant
names(morpho_all_sum_plant)
names(Iso_N_Aci_chl)
sunflower_all_data = left_join(Iso_N_Aci_chl, morpho_all_sum_plant, by = c("ID_trt", "n_trt", "drought", "plant", "time"))
Iso_N_Aci_chl$time
morpho_all_sum_plant$time
names(sunflower_all_data)
View(sunflower_all_data)
#### calculate sum biomass ###########
sunflower_all_data$Abg_all = sunflower_all_data$total_leaves+sunflower_all_data$stem+sunflower_all_data$rep_org+sunflower_all_data$senescence
sunflower_all_data$Abg_t4= sunflower_all_data$total_leaves+sunflower_all_data$stem+sunflower_all_data$rep_org
sunflower_all_data$biomass_t4= sunflower_all_data$Abg_t4+sunflower_all_data$roots
sunflower_all_data$biomass_all= sunflower_all_data$Abg_all+sunflower_all_data$roots
sunflower_all_data$RS_t4 = sunflower_all_data$roots/sunflower_all_data$Abg_t4
sunflower_all_data$RS_all = sunflower_all_data$roots/sunflower_all_data$Abg_all
sunflower_all_data$root_total= sunflower_all_data$roots/sunflower_all_data$biomass_t4
sunflower_all_data$leaves_stem_ratio = sunflower_all_data$total_leaves/sunflower_all_data$stem
sunflower_all_data$Abg_flowers_ratio = sunflower_all_data$rep_org /sunflower_all_data$Abg_t4
sunflower_all_data$stem_and_roots = sunflower_all_data$stem + sunflower_all_data$roots

########### calculation of the cost of transpiration ###################
######### roots and stems 
sunflower_all_data$cost_trans = sunflower_all_data$stem_and_roots/sunflower_all_data$transpiration
hist(sunflower_all_data$cost_trans)
######## roots
sunflower_all_data$Ecost = sunflower_all_data$roots/sunflower_all_data$transpiration
hist(sunflower_all_data$Ecost)

sunflower_all_data$ratio_J_V = sunflower_all_data$jmax_tleaf/sunflower_all_data$vcmax_tleaf
hist(sunflower_all_data$ratio_J_V )

names(sunflower_all_data)
sunflower_all_data$N_total_leaves = sunflower_all_data$Nmass * sunflower_all_data$total_leaves
sunflower_all_data$Ncost = sunflower_all_data$roots/sunflower_all_data$N_total_leaves
sunflower_all_data$WUE = sunflower_all_data$Abg_all/sunflower_all_data$transpiration
sunflower_all_data$NUE = sunflower_all_data$Abg_all/sunflower_all_data$Cumul_N

sunflower_all_data$betaplant = sunflower_all_data$Ncost/sunflower_all_data$cost_trans


########## calculate the standardized cumulative drought severity (Cumul-DD)  ########
sunflower_all_data$Sdt_DD = sunflower_all_data$Cumul_DD/sunflower_all_data$DAD

sunflower_all_data$n_trt
sunflower_all_data$N_level[sunflower_all_data$n_trt == 'LN'] <- '70 ppm'
sunflower_all_data$N_level[sunflower_all_data$n_trt == 'HN'] <- '210 ppm'
sunflower_all_data$N_level <- as.factor(sunflower_all_data$N_level)
sunflower_all_data$block <- factor(sunflower_all_data$block)
sunflower_all_data$plant <- factor(sunflower_all_data$plant)
sunflower_all_data$n_trt <- factor(sunflower_all_data$n_trt)
sunflower_all_data$drought <- factor(sunflower_all_data$drought)
sunflower_all_data$time <- factor(sunflower_all_data$time)
sunflower_all_data$Cumul_N <- as.numeric(sunflower_all_data$Cumul_N)
sunflower_all_data$Cumul_FC <- as.numeric(sunflower_all_data$Cumul_FC)
sunflower_all_data$Cumul_DD <- as.numeric(sunflower_all_data$Cumul_DD)
sunflower_all_data$treatment <- interaction(sunflower_all_data$n_trt, sunflower_all_data$drought, sep = " x ")

sunflower_all_data$Sdt_DD[is.na(sunflower_all_data$Sdt_DD)] <- 0

#############################################################################################################
#############################################################################################################
write.csv(sunflower_all_data, "../output/sunflower.csv", row.names = FALSE)
names(sunflower_all_data)

