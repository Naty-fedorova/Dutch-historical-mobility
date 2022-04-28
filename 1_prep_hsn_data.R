# code for raw data cleaning and checks, and preparation of working data to be used in analyses


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




# address data

# load in household address data
w_ha <- read.csv("./Working_files/Working_data_ha.csv", stringsAsFactors = FALSE)

# streamline to a dataset of person_id, household_id, start date, address, municipality

w_ha <- w_ha[ , c("person_id", "household_id",  "address_start_d",  "address_start_m", "address_start_y","address_end_d",  "address_end_m", "address_end_y", "municipality")]

expect_true(length(unique(w_ha$person_id)) == 36706)
expect_true(nrow(w_ha) == 338766)

# birth and death data
w_rp <- read.csv("./Working_files/Working_data_rp.csv", stringsAsFactors = FALSE)

# streamline
w_rp <- w_rp[ , c("person_id", "birth_y", "obsr_y", "death_y", "obsr_end_y", "sex")]

expect_true(nrow(w_rp) == 37173)
expect_true(!any(duplicated(w_rp$person_id)))

# birth y occuring after death y 
invalid_rps <- w_rp$person_id[which(w_rp$birth_y > w_rp$death_y)]
expect_true(length(invalid_rps) == 10)

# these individuals most likely represent errors in rounding and are removed
w_rp <- subset(w_rp, !w_rp$person_id %in% invalid_rps)
expect_true(nrow(w_rp) == 37163)
w_ha <- subset(w_ha, !w_ha$person_id %in% invalid_rps)
expect_true(nrow(w_ha) == 338679)

# unique rp ids
# check these against rp list
expect_true(all(w_ha$person_id %in% w_rp$person_id))

# merge
# create w_df set
w_df <- merge(w_rp, w_ha, by = "person_id")
expect_true(nrow(w_df) == 338679)

# dealing with age at move

# order moved by year and month, with each person_id
w_df <- w_df[with(w_df, order(person_id, address_start_y, address_start_m)), ]
w_df$age_at_move <- w_df$address_start_y - w_df$birth_y

# the two moves at -1868 represent an address start year of 0 - i.e. missing data
age_errors <- w_df[which(w_df$age_at_move < 0 ),]

# the errors account for 
expect_true(length(unique(age_errors$person_id)) == 29174)  # n of rps with suspicious ages
expect_true(nrow(age_errors) == 43738)    # n of suspicious cases

age_start_first_residence <- tapply(w_df$age_at_move, w_df$person_id, min)
neg_max_age_first_residence <- tapply(w_df$age_at_move, w_df$person_id, function(z) z[max(which(z <= 0))])
# this gives warnings b/c its NA for anyone who *doesn't* have an age at first residence

age_start_first_residence[which(!is.na(neg_max_age_first_residence))] <- neg_max_age_first_residence[which(!is.na(neg_max_age_first_residence))]

# create a `to_remove` flag for entries that occur before birth in the first residence
w_df$to_remove <- w_df$age_at_move < age_start_first_residence[w_df$person_id]
targets <- which(w_df$age_at_move == age_start_first_residence[w_df$person_id] & w_df$age_at_move < 0)
w_df$address_start_y[targets] <- w_df$birth_y[targets]
w_df$age_at_move[targets] <- 0


# remove error terms
w_df <- w_df[which(w_df$to_remove == FALSE),]

# recalculate reg event ages
w_df$age_at_move <- w_df$address_start_y - w_df$birth_y

expect_true(nrow(w_df[which(w_df$age_at_move < 0 ),]) == 0)
expect_true(nrow(w_df) == 320933)

# given this messes up the move order, we need to reorder again
# order moved by year and month, with each person_id
w_df <- w_df[with(w_df, order(person_id, address_start_y, address_start_m)), ]

# we also check the presence of reg events after death or end of obsr period
end_errors <- w_df[which((w_df$address_start_y > w_df$death_y) | (w_df$address_start_y > w_df$obsr_end_y) ),]

expect_true(nrow(end_errors) == 11882)
expect_true(length(unique(end_errors$person_id)) == 4224)

# we postulate that the errors have the same source - addresses for the hh are logged even after rps death, and kept in rps log
# although it is difficult to decide whether address logs or death/observation end years are more accurate
# we remove end_errors
# for additional conservatism, we use observation end years for this removal

w_df <- w_df[which(w_df$address_start_y <= w_df$obsr_end_y),]

expect_true(nrow(w_df[which((w_df$address_start_y > w_df$death_y) | (w_df$address_start_y > w_df$obsr_end_y) ),]) == 0)

# add order of move

w_df$nmove <- NA

prev_id <- w_df$person_id[1]
prev_rolling_count <- 0

for(i in 1:nrow(w_df)){
  if(w_df$person_id[i] == prev_id){
    w_df$nmove[i] <- prev_rolling_count + 1
    prev_rolling_count <- prev_rolling_count + 1
  }
  else {
    prev_id <- w_df$person_id[i]
    prev_rolling_count <- 1
    w_df$nmove[i] <- prev_rolling_count
  }
  if (i %% 1000 == 0) print(i)
}


# move categorization loop

prev_id <- ""
last_mun <- ""
natal_mun <- ""
prev_muns <- c()

# w_df is ordered by person_id and nmove
for(i in 1:nrow(w_df)){
  
  # since we haven't remove any of the missing data, will just skip over the empty ones here
  if(is.na(w_df$municipality[i])){
    w_df$mun_move_category[i] <- NA
  } else {
    
    if(w_df$person_id[i] != prev_id){
      # first registered address
      w_df$mun_move_category[i] <- 'first'
      prev_id <- w_df$person_id[i]
      last_mun <- w_df$municipality[i]
      natal_mun <- w_df$municipality[i]
      prev_muns <- c(w_df$municipality[i])
    } else {
      if(w_df$municipality[i] == last_mun){
        # RP moved to this address from within the same municipality
        w_df$mun_move_category[i] <- 'within'
      } else {
        last_mun <- w_df$municipality[i]
        if(!(w_df$municipality[i] %in% prev_muns)){
          # RP moved here from a different mun, and this is the first time they have come here
          w_df$mun_move_category[i] <- 'new mun'
          prev_muns <- c(prev_muns, w_df$municipality[i])
        } else {
          if(w_df$municipality[i] == natal_mun){
            # RPs natal mun, they have come back to it
            w_df$mun_move_category[i] <- 'natal return'
          } else {
            # RPs return to a mun they have lived in in the past
            w_df$mun_move_category[i] <- 'non-natal return'
          }
        }
      } 
    }
  }
  if (i %% 1000 == 0) print(i)
}


# we can now save w_df as the set that can be used for models
write.csv(w_df, file = "w_df.csv", row.names = FALSE)

## creating df for model 

df <- read.csv(file = "w_df.csv", stringsAsFactors = FALSE)

expect_true(nrow(df) == 309051)
expect_true(sum(is.na(df$obsr_end_y)) == 0)

df %>%
  group_by(person_id, address_start_y) %>%
  summarize(
    n_moves = n(),
    age = mean(age_at_move),
    b_y = first(birth_y),
    obs_end = first(obsr_end_y),
    gender = first(sex)
  ) %>%
  complete(
    address_start_y = seq(min(address_start_y), max(obs_end), by = 1),    
    nesting(person_id), fill = list(n_moves = 0)
  )  %>%
  mutate(age=age[1] + 1*(0:(length(age)-1))) %>%
  fill(b_y) %>%
  fill(obs_end) %>%
  fill(gender) %>%
  as.data.frame(stringsAsFactors = FALSE) -> s_person_year

s_person_year %>%
  mutate(gender = recode(gender,v = 1, m = 2)) %>%
  as.data.frame(stringsAsFactors = FALSE) -> s_person_year
  


# remove first address if age is 0 (i.e. first reg at birth) to turn table from reg events to move events
v <- which(s_person_year$age == 0)
s_person_year$n_moves[v] <- s_person_year$n_moves[v] - 1

expect_true(all(s_person_year$person_id %in% df$person_id))
expect_true(nrow(s_person_year) == 1078279)
expect_true(length(unique(s_person_year$person_id)) == 36595)
expect_true(length(unique(s_person_year$person_id)) == length(unique(df$person_id)))

# save person year table for use in model
write.csv(s_person_year, "s_person_year_df.csv", row.names = FALSE)


## creating a pure lifecourse set

w_df <- read.csv(file = "w_df.csv", stringsAsFactors = FALSE)

# remove RPs whose first registered address does not match the year of birth

s <- w_df[which((w_df$nmove == 1) & (w_df$address_start_y == w_df$birth_y)), ] 
s <- subset(w_df, w_df$person_id %in% s$person_id)

expect_true(nrow(s) == 285286)
expect_true(length(unique(s$person_id)) == 33820)

# s now contains the registration events for length(unique(s$person_id)) RPS, who are tracked from the birth year onwards

# because many death years are missing and there is a lot of variation between individuals with death years and without them, we consider only individuals with both birth and death years here

# rps with death year
temp_s <- subset(s, s$death_y != "NA")

# rps without death y 
temp_ss <- s[is.na(s$death_y),]

# compare distribution of death and obsvr end years
par(mfrow=c(2,1), mar = c(4, 4, 2, 2) )

plot(table(temp_s$death_y - temp_s$birth_y), ylab = "death y present")  # note! 33 cases where death year is before birth year - these will need to be removed
plot(table(temp_ss$obsr_end_y - temp_ss$birth_y), ylab = "death y missing")

# remove rps without death year
s <- s[which(!is.na(s$death_y)),]

expect_true(sum(is.na(s$death_y)) == 0)

# save s, the lifecourse set 
write.csv(s, "w_s.csv")


