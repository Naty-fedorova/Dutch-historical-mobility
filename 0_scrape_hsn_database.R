

# dependencies
library(foreign)
library(viridis)
library(rethinking)
library(tidyverse)
library(testthat)


## reading in data

# BEVOP - research persons 
rp <- read.dbf("./Data_files/BEVOP.DBF")
new_col_names_bevop <- c("rp_id", "region", "cohort", "rel_to_h", "f_rel_to_h", "m_rel_to_h", "cert_files", "birth_d", "birth_m", "birth_y", "birth_date", "death_d", "death_m", "death_y", "death_date", "death_source", "marriage_d", "marriage_m", "marriage_y", "marriage_date", "marriage_obsr_d", "marriage_obsr_m", "marriage_obsr_y", "marriage_obsr_date", "s_rel_to_h", "sex", "obsr_d", "obsr_m", "obsr_y", "obsr_date", "obse_end_d", "obsr_end_m", "obsr_end_y", "obsr_end_date", "gap_n", "no_rp", "release", "correction_date")
names(rp) <- new_col_names_bevop

# BEVHUISH - households
hh <- read.dbf("./Data_files/BEVHUISH.DBF")
new_col_names_bevhuish <- c("rp_id", "household", "household_pk", "link", "link_h", "obsr_d", "obsr_m", "obsr_y", "obsr_date", "obsr_end_d", "obsr_end_m", "obsr_end_y", "obsr_end_date", "release", "correction_date")
names(hh) <- new_col_names_bevhuish

# BEVADRES - household addresses
ha <- read.dbf("./Data_files/BEVADRES.DBF")
new_col_names_bevadres <- c("rp_id", "household", "household_rk", "version", "address_start_d", "address_start_m", "address_start_y", "address_start_date", "address_order", "address_end_d", "address_end_m", "address_end_y", "address_end_date", "pers_idf", "street", "house_n", "house_add", "district", "district_house_n", "district_house_add", "address_name_change", "owner_board_house", "institution", "boat", "boat_loc", "loc_municipality", "municipality", "release", "correction_date")
names(ha) <- new_col_names_bevadres

# BEVSTATP - static data (on all individuals)
sd<- read.dbf("./Data_files/BEVSTATP.DBF")
new_col_names_bevstatp <- c("rp_id", "household", "household_rk", "n_pph", "static_pk", "link", "person_id", "father_id", "father_cert", "mother_id", "mother_cert", "link_gen", "obsr_d", "obsr_m", "obsr_y", "obsr_date", "obsr_end_d", "obsr_end_m", "obsr_end_y", "obsr_end_date", "role", "surname", "name", "sex", "birth_d", "birth_m", "birth_y", "birth_date", "birth_municipality", "death_d", "death_m", "death_y", "death_date", "death_municipality", "release", "correction_date")
names(sd) <- new_col_names_bevstatp


## new ID creation

set.seed(1)

id_maker <- function(n, reserved='', seed=NA, nchars=NA){
  my_let <- letters 
  my_num <- 0:9 
  if(is.na(seed) | !is.numeric(seed)) set.seed(as.numeric(as.POSIXlt(Sys.time())))
  if(!is.na(seed) & is.numeric(seed)) set.seed(seed)
  output <- replicate(n, paste(sample(c(my_let, my_num), nchars, replace=TRUE), 
                               collapse=''))
  rejected <- duplicated(output) | output %in% reserved | substr(output, 1, 1) %in% my_num
  while(any(rejected)){
    output <- output[-which(rejected)]
    remaining <- n-length(output)
    output <- c(output, replicate(remaining, paste(sample(c(my_let, my_num), nchars, 
                                                          replace=TRUE), collapse='')))
    rejected <- duplicated(output) | output %in% reserved | substr(output, 1, 1) %in% my_num
  }
  output
}

expect_true(all(rp$rp_id %in% sd$rp_id)) 
expect_true(all(sd$rp_id %in% rp$rp_id)) 

sd$person_id <- paste(sd$rp_id, sd$person_id, sep = "-")
sd$father_id <- paste(sd$rp_id, sd$father_id, sep = "-")
sd$mother_id <- paste(sd$rp_id, sd$mother_id, sep = "-")

sd$person_id_fam <- sd$person_id
sd$mother_id_fam <- sd$mother_id
sd$father_id_fam <- sd$father_id

# we think 2018-us did this b/c "-0" is the code for missingness
sd$father_id[grep("-0$", sd$father_id)] <- NA
sd$mother_id[grep("-0$", sd$mother_id)] <- NA

unique_people <- unique(sd$person_id)
new_ids <- id_maker(n = length(unique_people), nchars = 5, seed = 1)
sd$person_id <- new_ids[match(sd$person_id, unique_people)]

sd$mother_id <- sd$person_id[match(sd$mother_id_fam, sd$person_id_fam)]
sd$father_id <- sd$person_id[match(sd$father_id_fam, sd$person_id_fam)]

# 'OP' indicates they are research person, and 'ALLEENST' indicates they are the only person in the household so they must be the research person!
sd$is_rp <- sd$role %in% c("HOOFD=OP", "VROUW=OP", "OP", "ALLEENST")

# check to make sure this makes sense
sd[which(sd$rp_id == sample(sd$rp_id, 1)), c("rp_id", "household", "sex", "person_id", "father_id", "mother_id", "obsr_date", "obsr_end_date", "is_rp")]

## ID matching

# in order to match person_id with specific persons we can use rp_id and role. Role of OP(dutch for rp) is given by role is hoofd=op, vrouw=op, op, and alleenst and captured in the variable is_rp calculated in chunk above

# index of rps
rp_index <- which(sd$role %in% c("HOOFD=OP", "VROUW=OP", "OP", "ALLEENST"))

# check - make sure indexing works
# nrow d should be the same as research_persons
t <- unique(sd[rp_index, c("rp_id", "person_id")])
expect_true(nrow(t) == nrow(rp))
expect_true(nrow(t) == 37173)

t <- sd[rp_index, c("rp_id", "person_id")]

# check - need to make sure each person_id is associated with only 1 rp_id
# nothing should print

for(i in 1:nrow(t)){
  temp <- t[which(t$rp_id == t[i, 1]),]
  temp2 <- unique(temp$person_id) 
  
  if(length(temp2) > 1){
    print(temp)  }
}


## match rps with person id
# create working files that will then be written to the folder, so that original data files don't need to be used

# sd already done
w_sd <- sd

# match rps with person id
rp$person_id <- unique(w_sd$person_id[rp_index])
w_rp <- rp

# matching for household addresses (noted per rp)
expect_true(sum(is.na(ha$rp_id)) == 0)
expect_true(sum(is.na(w_rp$rp_id)) == 0)
ha$person_id <- w_rp$person_id[match(ha$rp_id, w_rp$rp_id)]
w_ha <- ha  

mean(ha$household_rk %in% hh$household_pk, na.rm = TRUE)

## ID maker for households

# make correct number of unique ids
unique_households <- unique(hh$household_pk)
new_household_ids <- id_maker(n = length(unique_households), nchars = 4, seed = 2)

# match ids to households in hh
hh$household_id <- new_household_ids[match(hh$household_pk, unique_households)]
w_hh <- hh

# these should be the same, because hh is row per household
expect_true(length(unique(hh$household_id)) == 86777)
expect_true(nrow(hh) == 86777)

# matching for household addresses
expect_true(sum(is.na(hh$household_pk)) == 0)
expect_true(sum(is.na(ha$household_rk)) == 0)
ha$household_id  <- hh$household_id[match(ha$household_rk, hh$household_pk)]
w_ha <- ha

# need to check that household_rk is indeed a subset of household_pk:
mean(ha$household_rk %in% hh$household_pk, na.rm = TRUE)

# match with sd
w_sd$household_id <- w_hh$household_id[match(sd$household_rk, w_hh$household_pk)]

# streamlining
w_sd <- subset(w_sd, select = -c(household_rk, link, father_cert, mother_cert, link_gen, surname, name, release, correction_date))
w_rp <- subset(w_rp, select = -c(region, cohort, cert_files, release, correction_date))
w_hh <- subset(w_hh, select = -c(household_pk, link, link_h, release, correction_date))
w_ha <- subset(w_ha, select = -c(household_rk, version, pers_idf, release, correction_date))


## saving working files
dir.create("Working_files")
setwd("./Working_files")

# w_rp
write.csv(w_rp, file ="Working_data_rp.csv")

# w_sd
write.csv(w_sd, file ="Working_data_sd.csv")

# w_hh
write.csv(w_hh, file ="Working_data_hh.csv")

# w_ha
write.csv(w_ha, file ="Working_data_ha.csv")
