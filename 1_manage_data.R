# code for raw data cleaning and checks, and preparation of working data to be used in analyses

# dependencies
library(foreign)
library(viridis)
library(rethinking)
library(tidyverse)
library(testthat)

## data cleaning and creation of datasets for analysis

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
  mutate(sex = recode(sex,v = 1, m = 2)) %>%  
  group_by(person_id, address_start_y) %>%
  summarize(
    n_moves = n(),
    age = mean(age_at_move),
    b_y = first(birth_y),
    obs_end = first(obsr_end_y)
  ) %>%
  complete(
    address_start_y = seq(min(address_start_y), max(obs_end), by = 1),    
    nesting(person_id), fill = list(n_moves = 0)
  )  %>%
  mutate(age=age[1] + 1*(0:(length(age)-1))) %>%
  fill(b_y) %>%
  fill(obs_end) %>%
  as.data.frame(stringsAsFactors = FALSE) -> s_person_year


# remove first address if age is 0 (i.e. first reg at birth) to turn table from reg events to move events
v <- which(s_person_year$age == 0)
s_person_year$n_moves[v] <- s_person_year$n_moves[v] - 1

expect_true(all(s_person_year$person_id %in% df$person_id))
expect_true(nrow(s_person_year) == 1078279)
expect_true(length(unique(s_person_year$person_id)) == 36595)
expect_true(length(unique(s_person_year$person_id)) == length(unique(df$person_id)))

# save person year table for use in model
write.csv(s_person_year, "s_person_year_df.csv")

## creating a pure lifecourse set

# remove RPs whose first registered address does not match the year of birth

s <- w_df[which((w_df$nmove == 1) & (w_df$address_start_y == w_df$birth_y)), ] 
s <- subset(w_df, w_df$person_id %in% s$person_id)

expect_true(nrow(s) == 285286)
expect_true(length(unique(s$person_id)) == 33820)

# remove rps without death year
s <- s[which(!is.na(s$death_y)),]

expect_true(sum(is.na(s$death_y)) == 0)

# save s, the lifecourse set 
write.csv(s, "w_s.csv")
