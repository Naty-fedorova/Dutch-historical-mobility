# code to simulate data for use in analysis

# dependencies 
library(tidyverse)

# simulating w_df 
n <- 100    # number of RPs to simulate for, note that for the cohort regressions to work, you will need atleast 10000 rps
ID <- 1:n
birth_y <- sample(1850:1920, n, replace = TRUE)  # sample birth years, not reproducing birth year proportions of HSN data
df_identities <- data.frame(ID, birth_y)
columns <- c("birth_y", "person_id", "nmove", "address_start_y")
df_sim <- data.frame(matrix(ncol = length(columns), nrow = 0))
colnames(df_sim) <- columns

for(i in 1:n){
  n_move <- ceiling(rexp(1, rate = 1/5))
  
  for (nmove in 1:n_move){
    if (nmove == 1) {
      address_start_y <- df_identities$birth_y[i] + abs(round(rnorm(1)))
    }
    else {
      address_start_y <- df_sim$address_start_y[nrow(df_sim)] + round(runif(1, min = 1, max = 10))
    }
    df_sim[nrow(df_sim) + 1, ] <- list(df_identities$birth_y[i], df_identities$ID[i], nmove, address_start_y)
  }
}


# add obsr_end_year
df_sim$obsr_end_y <- 0

for(id in 1:n){
  sub <- df_sim[which(df_sim$person_id == id),]
  df_sim$obsr_end_y[which(df_sim$person_id == id)] <- max(sub$address_start_y)  + abs(round(rnorm(1, 5, 2)))
}

# add age at move
df_sim$age_at_move <- df_sim$address_start_y - df_sim$birth_y 


# creating s_person_year as per modified code from data_management.R

df_sim %>%
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
  as.data.frame(stringsAsFactors = FALSE) -> s_person_year_sim

write.csv(s_person_year_sim, "s_person_year_sim.csv")


