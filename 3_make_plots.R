# dependencies
library(rethinking)
library(viridis)
library(testthat)
library(tidyverse)
# for color palette
library(MetBrewer)

# create folder to contain plots
dir.create("Figures")

# load in data if necessary
# simulated data
d <- read.csv("d_sim.csv", stringsAsFactors = FALSE)

# when working with real data:
# df used for model
#d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

# df of birth-death lifecourses
#s <- read.csv(file = "w_s.csv", stringsAsFactors = FALSE)

# load in samples if necessary
# load("post_pois.RData")

# set plotting colors
cols <- met.brewer(name = "Tam", 8)[c(5,8)]  
col_neut <- met.brewer(name = "Tam", 8)[3]
cols_int <- adjustcolor(cols, 0.6)
cols_int_neut <- adjustcolor(col_neut, 0.4)

# subset to plot for
person_ids <- sort(unique(d$person_id)) 
n_rp <- 36595  #36595 for full set
rp_sub <- sample(person_ids, size = n_rp)
dm_sim <- subset(d, d$person_id %in% rp_sub)
uniq_ages <- length(unique(dm_sim$age))

#-----------------------------------------------------------------------------------------------------------------------
### plotting total residential moves over the lifetime (only possible with HSN data)

#check first reg address is at age 0
expect_true(nrow(s[which((s$mun_move_category == "first") & (s$age_at_move != 0)),]) == 0)   # first reg address is always at age 0

ss <- s[which((s$death_y - s$birth_y) >= 20),]

ss %>%
  mutate(gender = recode(sex,v = 1, m = 2)) -> ss  # female (vrou) = 1, male (mann) = 2

ss %>%
  group_by(person_id) %>%
  summarize(
    total_moves = length(nmove),
    gender = first(gender)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> ss_total_moves

table(ss_total_moves$total_moves) # this shows us the number of registered addresses logged per person
table(ss_total_moves$total_moves, ss_total_moves$gender) # and split by gender

dist <- table(ss_total_moves$total_moves - 1) # this shows us the number of moves, as for the individuals with just one registered address, they are now logged as having zero moves

# for range and IQR
summary(dist)

# percent of individuals that have more than 30 moves
sum(dist[31:82])/sum(dist)*100

#life_mob_gen <- table(ss_total_moves$gender, (ss_total_moves$total_moves - 1))
#temp <- as.data.frame.matrix(temp)  


# Variable am to factor
gender <- factor(ss_total_moves$gender)
n_moves <- factor(ss_total_moves$total_moves)  


# Change factor levels
levels(gender) <- c("1", "2")
levels(n_moves) <- sort(unique(ss_total_moves$total_moves))


# Table cylinder - transmission type
life_mob_gen <- table(gender, n_moves)


# stacked bar graph to visualize mobility over age with split by gender

png("Figures/hist_s_total_moves_gen.png",type = "cairo", res = 300, height = 15, width = 20, units = "cm")
blank()
plot(x = 0:130,
     y = 0:130,
     xlim = c(0,130),
     ylim = c(0,250),
     xaxt = "n",
     main = "Total number of moves recorded over lifetime, by gender",
     xlab = "Total lifetime moves",
     ylab = "Frequency of total move category",
     bty = "n",
     type = "n")

for(i in 1:131){
  x <- i
  
  if(x %in% unique(ss_total_moves$total_moves)){
    y <- which(colnames(life_mob_gen) == x)
    
    # x - 1 to change from registrations to moves  
    lines(x = c(x-1, x-1) ,y = c(0, life_mob_gen[1,y]), col = cols[1], lwd = 5)
    lines(x = c(x-1, x-1) ,y = c(life_mob_gen[1,y], life_mob_gen[1,y]+life_mob_gen[2,y]), col = cols[2], lwd = 5)
  } else{
    lines(x = c(x-1, x-1) ,y = c(0, 0), col = cols[1])
    lines(x = c(x-1, x-1) ,y = c(0, 0), col = cols[2])
  }
}

# female median
abline(v = median(ss_total_moves$total_moves[which(ss_total_moves$gender == 1)]-1), col = "mistyrose3", lty = 2, lwd = 2)

# male median
abline(v = median(ss_total_moves$total_moves[which(ss_total_moves$gender == 2)]-1), col = "thistle4", lty = 2, lwd = 2)

axis(side = 1, seq(from = 0, to = 130, by = 10))
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### Simulating and plotting predicted moves |post-stratification|

## post-stratification females only

dm_sim_f <- dm_sim[which(dm_sim$gender == "1"),]

uniq_ages_f <- length(unique(dm_sim$age))

n_individuals_f <- rep(0, uniq_ages_f)
real_moves_f <- rep(0, uniq_ages_f)
est_moves_mean_f <- rep(0, uniq_ages_f)
est_moves_hpdi_f <- matrix(0, 2, uniq_ages_f)
est_moves_pi_25_f <- matrix(0, 2, uniq_ages_f)
est_moves_pi_50_f <- matrix(0, 2, uniq_ages_f)
est_moves_pi_75_f <- matrix(0, 2, uniq_ages_f)
age_indices_f <- seq(from = 1, to = uniq_ages_f, by = 1)

pb <- txtProgressBar(min(age_indices_f), max(age_indices_f), style=3)

for(age_i in age_indices_f) {
  
  est_moves_sum <- rep(0, 2000)
  
  dm_individuals_at_age <- dm_sim_f[which(dm_sim_f$age == age_i - 1),]
  
  if (nrow(dm_individuals_at_age) == 0) next # should never happen on whole dataset
  
  # loop through all persons at given age
  for(ind in 1:nrow(dm_individuals_at_age)) {
    
    # Simply counting the number of individuals per age group
    # to check that it looks reasonable
    n_individuals_f[age_i] <- n_individuals_f[age_i] + 1
    
    # get real moves for individual
    n_moves <- dm_individuals_at_age[ind, "n_moves"]
    # add real moves for individual to real_moves[age_i]
    real_moves_f[age_i] <- real_moves_f[age_i] + n_moves
    
    # calculate N lambdas
    person_id <- dm_individuals_at_age[ind, "person_id"]
    lambdas <- exp(post_pois$mu + rnorm(2000,0,post_pois$sd_id) + post_pois$beta[, 1, age_i])
    
    # calculate N est_moves for individual
    est_moves <- sapply(lambdas, function(i) rpois(1, i) )
    
    # Sum est_moves across individuals, keeping the sample dimension
    est_moves_sum <- est_moves_sum + est_moves # est_moves is for one person, here we sum all people in age group
    
    # add mean of est moves
    est_moves_mean_f[age_i] <- est_moves_mean_f[age_i] + mean(est_moves)
  }
  
  # calculate HPDI est_moves for individual
  hpdi <- HPDI(est_moves_sum, prob=0.89) 
  PI_25 <- PI(est_moves_sum, prob=0.25)
  PI_50 <- PI(est_moves_sum, prob=0.50)
  PI_75 <- PI(est_moves_sum, prob=0.75)
  
  # add HPDI est_moves for individual to est_moves_hpdi
  est_moves_hpdi_f[1, age_i] <- est_moves_hpdi_f[1, age_i] + hpdi[1]
  est_moves_hpdi_f[2, age_i] <- est_moves_hpdi_f[2, age_i] + hpdi[2]
  
  # add PI est moves for individual to est_moves_pi
  est_moves_pi_25_f[1, age_i] <- est_moves_pi_25_f[1, age_i] + PI_25[1]
  est_moves_pi_25_f[2, age_i] <- est_moves_pi_25_f[2, age_i] + PI_25[2]
  
  est_moves_pi_50_f[1, age_i] <- est_moves_pi_50_f[1, age_i] + PI_50[1]
  est_moves_pi_50_f[2, age_i] <- est_moves_pi_50_f[2, age_i] + PI_50[2]
  
  est_moves_pi_75_f[1, age_i] <- est_moves_pi_75_f[1, age_i] + PI_75[1]
  est_moves_pi_75_f[2, age_i] <- est_moves_pi_75_f[2, age_i] + PI_75[2]
  
  setTxtProgressBar(pb, age_i)
}
close(pb)

# post-stratification males only

dm_sim_m <- dm_sim[which(dm_sim$gender == "2"),]

uniq_ages_m <- length(unique(dm_sim$age))

n_individuals_m <- rep(0, uniq_ages_m)
real_moves_m <- rep(0, uniq_ages_m)
est_moves_mean_m <- rep(0, uniq_ages_m)
est_moves_hpdi_m <- matrix(0, 2, uniq_ages_m)
est_moves_pi_25_m <- matrix(0, 2, uniq_ages_m)
est_moves_pi_50_m <- matrix(0, 2, uniq_ages_m)
est_moves_pi_75_m <- matrix(0, 2, uniq_ages_m)
age_indices_m <- seq(from = 1, to = uniq_ages_m, by = 1)

pb <- txtProgressBar(min(age_indices_m), max(age_indices_m), style=3)

for(age_i in age_indices_m) {
  
  est_moves_sum <- rep(0, 2000)
  
  dm_individuals_at_age <- dm_sim_m[which(dm_sim_m$age == age_i - 1),]
  
  if (nrow(dm_individuals_at_age) == 0) next # should never happen on whole dataset
  
  # loop through all persons at given age
  for(ind in 1:nrow(dm_individuals_at_age)) {
    
    # Simply counting the number of individuals per age group
    # to check that it looks reasonable
    n_individuals_m[age_i] <- n_individuals_m[age_i] + 1
    
    # get real moves for individual
    n_moves <- dm_individuals_at_age[ind, "n_moves"]
    # add real moves for individual to real_moves[age_i]
    real_moves_m[age_i] <- real_moves_m[age_i] + n_moves
    
    # calculate N lambdas
    person_id <- dm_individuals_at_age[ind, "person_id"]
    lambdas <- exp(post_pois$mu + rnorm(2000,0,post_pois$sd_id) + post_pois$beta[, 2, age_i])
    
    # calculate N est_moves for individual
    est_moves <- sapply(lambdas, function(i) rpois(1, i) )
    
    # Sum est_moves across individuals, keeping the sample dimension
    est_moves_sum <- est_moves_sum + est_moves # est_moves is for one person, here we sum all people in age group
    
    # add mean of est moves
    est_moves_mean_m[age_i] <- est_moves_mean_m[age_i] + mean(est_moves)
  }
  
  # calculate HPDI est_moves for individual
  hpdi <- HPDI(est_moves_sum, prob=0.89) 
  PI_25 <- PI(est_moves_sum, prob=0.25)
  PI_50 <- PI(est_moves_sum, prob=0.50)
  PI_75 <- PI(est_moves_sum, prob=0.75)
  
  # add HPDI est_moves for individual to est_moves_hpdi
  est_moves_hpdi_m[1, age_i] <- est_moves_hpdi_m[1, age_i] + hpdi[1]
  est_moves_hpdi_m[2, age_i] <- est_moves_hpdi_m[2, age_i] + hpdi[2]
  
  # add PI est moves for individual to est_moves_pi
  est_moves_pi_25_m[1, age_i] <- est_moves_pi_25_m[1, age_i] + PI_25[1]
  est_moves_pi_25_m[2, age_i] <- est_moves_pi_25_m[2, age_i] + PI_25[2]
  
  est_moves_pi_50_m[1, age_i] <- est_moves_pi_50_m[1, age_i] + PI_50[1]
  est_moves_pi_50_m[2, age_i] <- est_moves_pi_50_m[2, age_i] + PI_50[2]
  
  est_moves_pi_75_m[1, age_i] <- est_moves_pi_75_m[1, age_i] + PI_75[1]
  est_moves_pi_75_m[2, age_i] <- est_moves_pi_75_m[2, age_i] + PI_75[2]
  
  setTxtProgressBar(pb, age_i)
}
close(pb)


png("Figures/moves_ages_PI_gen.png",  type="cairo", res = 300, height = 30, width = 20, units = "cm")

par(mfrow = c(2,1))
blank()
# split by gender
plot(y = real_moves_f,
     x = age_indices_f-1,
     ylim = c(0,5000), 
     xlim = c(0,max(age_indices_f-1)),
     main = "Total moves per age, split by gender", 
     xlab = "Age", 
     ylab = "Total number of moves", 
     bty = "n", 
     col = cols[1],
     pch = 19,
     xaxt = "n", 
     yaxt = "n", 
     font.main = 1,
     type = "n")
axis(1, at = seq(0, 100, by = 10))
axis(2, at = seq(0, 5000, by = 1000))

# intervals for females
shade(est_moves_pi_75_f, age_indices_f - 1, col = cols_int[1])
shade(est_moves_pi_50_f, age_indices_f - 1, col = cols_int[1])
shade(est_moves_pi_25_f, age_indices_f - 1, col = cols_int[1])

# intervals for males
shade(est_moves_pi_75_m, age_indices_m - 1, col = cols_int[2])
shade(est_moves_pi_50_m, age_indices_m - 1, col = cols_int[2])
shade(est_moves_pi_25_m, age_indices_m - 1, col = cols_int[2])

# lines for female subset
lines(x = age_indices_f-1, y = real_moves_f, pch = 1, col = cols[1])

# lines for male subset
lines(x = age_indices_m-1, y = real_moves_m, pch = 1, col = cols[2])

# TODO: add legend
legend(x = 60, y = 4000, legend = c("Female", "Male"), lty = 1, col = c(cols[1], cols[2]), cex = 2, bty = "n", lwd = 5)

mtext('Plot A', side=3, line=0.5, at=85, font = 2, cex = 1.5)


# contrast between genders
post_strat_con <- est_moves_mean_f - est_moves_mean_m
post_strat_con_int <- est_moves_pi_25_f - est_moves_pi_25_m

plot(y = 0:(length(post_strat_con)-1),
     x = 0:(length(post_strat_con)-1), 
     ylim = c(-250, 1250),  
     xlim = c(0,(length(post_strat_con)-1)),
     xlab = "Age", 
     ylab = "Gender difference in total moves per age", 
     main = "Sex-based contrast for total mobility over ages", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     yaxt = "n",
     type = "n")
axis(1, at = seq(0,100, by = 10))
axis(2, at = seq(-250, 1250, by = 250))

shade( post_strat_con_int, 0:(length(post_strat_con_int[1,])-1), col = cols_int_neut)
lines(x = 0:(length(post_strat_con)-1), y = post_strat_con, col = col_neut, lty = 1, lwd = 2 )

# line for zero
lines(x = 0:(length(post_strat_con)-1), y = rep(0, (length(post_strat_con))), col = "grey", lty = 2)

mtext('Plot B', side=3, line=0.5, at=85, font = 2, cex = 1.5)

dev.off()
#-----------------------------------------------------------------------------------------------------------------------
### plotting beta estimate, implied number of moves per year from model estimates

# reformatting dm_sim to obtain averages from data
dm_sim %>%
  group_by(age) %>%
  summarize(
    total_moves = sum(n_moves),
    average_moves = mean(n_moves)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> mean_moves

# processing posterior
exp_mu_f <- sapply( 1:uniq_ages , function(i) exp( post_pois$mu + post_pois$beta[, 1,i] + rnorm(2000,0,post_pois$sd_id) ) )
exp_mu_m <- sapply( 1:uniq_ages , function(i) exp( post_pois$mu + post_pois$beta[, 2,i] + rnorm(2000,0,post_pois$sd_id) ) )

beta_mod_int_f <- apply(exp_mu_f, 2, HPDI)
beta_mod_int_m <- apply(exp_mu_m, 2, HPDI)

beta_mod_f <- apply(exp_mu_f, 2, mean)
beta_mod_m <- apply(exp_mu_m, 2, mean)

# contrasts
beta_con <- beta_mod_f - beta_mod_m
beta_int_con <- beta_mod_int_f - beta_mod_int_m

exp_mu_con <- exp_mu_f - exp_mu_m

# plotting
png("Figures/beta_estimates_gen.png", type = "cairo", res = 300, height = 30, width = 20, units = "cm")

par(mfrow = c(2,1))

# beta plot
blank()
plot(y = 0:(length(beta_mod_f)-1),
     x = 0:(length(beta_mod_f)-1), 
     ylim = c(0, 0.6),  
     xlim = c(0,(length(beta_mod_f)-1)),
     xlab = "Age", 
     ylab = "Estimated number of moves per year", 
     main = "Age-based estimated moves per year for each sex", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     type = "n")
axis(1, at = seq(0,90, by = 10))
axis(2, at = seq(0, max(beta_mod_int_f), by = 0.1))

#shade( apply(exp_mu_f, 2, PI, 0.75 ) , 0:(length(beta_mod_f)-1), col = cols_int[1])
shade( apply(exp_mu_f, 2, PI, 0.5 ) , 0:(length(beta_mod_f)-1), col = cols_int[1])
#shade( apply(exp_mu_f, 2, PI, 0.25 ) , 0:(length(beta_mod_f)-1), col = cols_int[1])

#shade( apply(exp_mu_m, 2, PI, 0.75 ) , 0:(length(beta_mod_m)-1), col = cols_int[2])
shade( apply(exp_mu_m, 2, PI, 0.5 ) , 0:(length(beta_mod_m)-1), col = cols_int[2])
#shade( apply(exp_mu_m, 2, PI, 0.25 ) , 0:(length(beta_mod_m)-1), col = cols_int[2])

points(y = mean_moves$average_moves, x = 0:(length(mean_moves$average_moves)-1), col = "black", pch = 19, cex = 0.5)
lines(x = 0:(length(beta_mod_f)-1), y = beta_mod_f, col = cols[1], lty = 2, lwd = 3 )
lines(x = 0:(length(beta_mod_m)-1), y = beta_mod_m, col = cols[2], lty = 2, lwd = 3 )

legend(x = 60, y = 0.6, legend = c("female", "male"), lty = 1, col = c(cols[1], cols[2]), cex = 1.5, bty = "n", lwd = 3)
legend(x = 60, y = 0.5, legend = c("data from sample"), pch = 19, col = "black", cex = 1.5, bty = "n")

mtext('Plot A', side=3, line=0.5, at=85, font = 2, cex = 1.5)

# contrast plot
blank()
plot(y = 0:(length(beta_mod_f)-1),
     x = 0:(length(beta_mod_f)-1), 
     ylim = c(-0.3, 0.3),  
     xlim = c(0,(length(beta_mod_f)-1)),
     xlab = "Age", 
     ylab = "Difference in estimated number of moves per year", 
     main = "Sex-based contrast between mobility over ages", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     yaxt = "n",
     type = "n")
axis(1, at = seq(0,90, by = 10))
axis(2, at = seq(-0.3, 0.3, by = 0.1), labels = c("-0.3", "-0.2", "-0.1", "0","0.1", "0.2", "0.3" ))

shade( apply(exp_mu_con, 2, PI, 0.5 ) , 0:(length(beta_mod_m)-1), col = cols_int_neut)
lines(x = 0:(length(beta_con)-1), y = beta_con, col = col_neut, lty = 1, lwd = 3 )

# line for zero
lines(x = 0:(length(beta_con)-1), y = rep(0, (length(beta_con))), col = "grey", lty = 2)

mtext('Plot B', side=3, line=0.5, at=85, font = 2, cex = 1.5)

dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting cohort model
#cohorts <- seq(from = 1850, to = 1922 , by = 1) # note: rewrites cohorts from 2_fit_modes.R

years_length <- length(cohorts) + 20     #ncol(m_coh_samples[[length(cohorts)]][["beta"]][,,1][,1])
coh_label <- seq(from = 1850, to = 1850 + years_length - 1 , by = 1)

# constructing matrix for plotting
coh_mat_f <- matrix(data = NA, nrow = 100, ncol = years_length+10)
coh_mat_m <- matrix(data = NA, nrow = 100, ncol = years_length+10)


for(x in 1:length(cohorts)){
  
  cohort_ages <- dim(m_coh_samples[[x]][["beta"]])
  cohort_ages <- cohort_ages[3]
  
  exp_mu_coh_f <- sapply( 1:cohort_ages , function(i) exp( m_coh_samples[[x]][["mu"]] + m_coh_samples[[x]][["beta"]][,,i][,1] + rnorm(2000,0,m_coh_samples[[x]][["sd_id"]]) ) )
  exp_mu_coh_m <- sapply( 1:cohort_ages , function(i) exp( m_coh_samples[[x]][["mu"]] + m_coh_samples[[x]][["beta"]][,,i][,2] + rnorm(2000,0,m_coh_samples[[x]][["sd_id"]]) ) )
  
  m_f <- apply(exp_mu_coh_f, 2, mean)
  m_m <- apply(exp_mu_coh_m, 2, mean)
  
  for(y in 1:length(m_f)){
    if (x+(y-1) < years_length) {
      # matrix construction for years of time on x axis
      coh_mat_f[y, x+(y-1)] <-  m_f[y]
    }
  }
  
  for(y in 1:length(m_m)){
    if (x+(y-1) < years_length) {
      # matrix construction for years of time on x axis
      coh_mat_m[y, x+(y-1)] <-  m_m[y]
    }
  }
  
}

coh_mat_f <- t(coh_mat_f)
coh_mat_m <- t(coh_mat_m)

# defining colors
grey <- rgb(0,0,0,alpha=0.2)


cols <- met.brewer(name = "Tam", 8)[c(5,8)]  #"VanGogh2"
col_neut <- met.brewer(name = "Tam", 8)[3]
cols_int <- adjustcolor(cols, 0.6)
cols_int_neut <- adjustcolor(col_neut, 0.4)


ramp_f <- colorRampPalette(c("mistyrose3", cols[1]))
cols_f <- ramp_f(1000)

ramp_m <- colorRampPalette(c("thistle3", cols[2]))
cols_m <- ramp_m(1000)



# note that observation ends at 1945 in HSN sample

# cohort heatmap
png("Figures/heatmap_gen.png", type = "cairo", res = 300, height = 40, width = 20, units = "cm")

# TODO: might also need to readjust figure margins
par(mfrow = c(2,1))

# females cohort plot
image(x = 0:(years_length - 1 + 10) , # + 10 for plotting
      y = 0:100, 
      z = coh_mat_f, 
      col = cols_f, 
      ylim = c(0,100), 
      xlab = "Birth year", 
      ylab = "Age", 
      bty = "n",
      xaxt = "n" )
axis(1, at= seq(0 ,years_length + 10, by = 10), labels= seq(1850, 1950, by = 10))
abline(h = c(20, 40, 60, 80), col = grey )
abline(v = seq(10, 100, by = 10), col = grey)
title("Mobility over age and time  \nFemales", font.main = 1)

mtext('Plot A', side=3, line=0.5, at=95, font = 2, cex = 1.5 )

#legend
points(x = rep(2, times = 1000), y = seq(from = 85, to = 95, length.out = 1000), col = cols_f)
text(x = 2, y = 82, "low mobility", adj = 0, cex = 0.8)     
text(x = 2, y = 98, "high mobility", adj = 0, cex = 0.8)  

# males 
image(x = 0:(years_length - 1 + 10) , # + 10 for plotting
      y = 0:100, 
      z = coh_mat_m, 
      col = cols_m, 
      ylim = c(0,100), 
      xlab = "Birth year", 
      ylab = "Age", 
      bty = "n",
      xaxt = "n" )
axis(1, at= seq(0 ,years_length + 10, by = 10), labels= seq(1850, 1950, by = 10))
abline(h = c(20, 40, 60, 80), col = grey )
abline(v = seq(10, 100, by = 10), col = grey)
title("Males", font.main = 1)

mtext('Plot B', side=3, line=0.5, at=95, font = 2, cex = 1.5 )

points(x = rep(2, times = 1000), y = seq(from = 85, to = 95, length.out = 1000), col = cols_m)
text(x = 2, y = 82, "low mobility", adj = 0, cex = 0.8)     
text(x = 2, y = 98, "high mobility", adj = 0, cex = 0.8)  

dev.off()

### Supplementary plots
#-------------------------------------------------------------------------------
# pois model estimates
png("Figures/model_estimates_gen.png", type = "cairo", res = 300, height = 6, width = 8, units = "cm")
plot(precis(stanfit_pois, pars = c("mu", "rho", "eta"), depth = 2), ylab = "Parameter", main ="Parameter estimates")
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# Rhat values against effective samples
png("Figures/rhat_neff_gen.png",type = "cairo", res = 300, height = 15, width = 20, units = "cm")
dashboard(stanfit_pois)
dev.off()

#------------------------------------------------------------------------------
# individual effects
exp_mu_a <- sapply( 1:36595 , function(i) exp( post_pois$a[,i]  ) )
a_mod <- apply(exp_mu_a, 2, mean)
a_mod_int <- apply(exp_mu_a, 2, PI, 0.5)

results <- matrix(0, 3, length(unique(dm_sim$person_id)))
results[1,] <- a_mod
results[2,] <- a_mod_int[1,]
results[3,] <- a_mod_int[2,]

# order by mean
results <- results[,order(results[1,])]

png("Figures/alphas.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = 0:(length(a_mod)-1),
     x = 0:(length(a_mod)-1), 
     ylim = c(0, max(a_mod_int)),  
     xlim = c(0,(length(a_mod)-1)),
     xlab = "RP index", 
     ylab = "Alpha estimate", 
     main = "Alpha estimate for each RP", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     yaxt = "n",
     type = "n")
axis(1, at = seq(0,40000, by = 10000))
axis(2, at = seq(0, 25, by = 5))

shade(results[c(2,3),] , 0:(length(a_mod)-1), col = vir_int[15])
lines(results[1,], col = vir[5])
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting individual trajectories - how RPs acquire moves over their observation periods
# create person tables and cumulative sum tables
dm_sim %>%
  group_by(person_id) %>%
  summarize(
    moves = sum(n_moves),
    gender = first(gender)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> dm_sim_person # person table with total number of lifetime moves

# cumulative counts of moves
dm_sim %>%
  group_by(person_id) %>%
  mutate(
    moves_cumsum = cumsum(n_moves),
  ) %>%
  ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE) -> dm_sim_cumsum

# split by gender
dm_sim_f <- dm_sim[which(dm_sim$gender == "1"),]
dm_sim_m <- dm_sim[which(dm_sim$gender == "2"),]

cols <- met.brewer(name = "Tam", 8)[c(5,8)]  
col_neut <- met.brewer(name = "Tam", 8)[3]
cols_int <- adjustcolor(cols, 0.6)
cols_int_neut <- adjustcolor(col_neut, 0.4)

ramp_f <- colorRampPalette(c("mistyrose2", cols[1]))
cols_f <- ramp_f(130)
plot_cols_f <- adjustcolor(cols_f, 0.7)

ramp_m <- colorRampPalette(c("thistle", cols[2]))
cols_m <- ramp_m(130)
plot_cols_m <- adjustcolor(cols_m, 0.7)



for(i in 1:nrow(dm_sim_person)){
  moves <- dm_sim_person$moves[i]
  if(dm_sim_person$gender[i] == "1"){
    dm_sim_person$col_real[i] <- plot_cols_f[moves + 1]
  } else{
    dm_sim_person$col_real[i] <- plot_cols_m[moves + 1]
  }
}



# Plot move accumulation
png("Figures/trajectories_nmoves_gen.png", type = "cairo" ,res = 300, height = 30, width = 15, units = "cm")

# Plotting for real trajectories
par(mfrow = c(2,1))
blank()

plot(dm_sim_cumsum$moves_cumsum ~ dm_sim$age ,
     type = "n", 
     xlab  = "Age", 
     ylab = "Number of moves", 
     main = "Cumulative totals of moves \nfor female RPs",
     font.main = 1,
     bty = "n",
     ylim = c(0,150),
     xlim = c(0,100))

# for females
for (i in 1:length(unique(dm_sim_f$person_id))){
  rp_id <- unique(dm_sim_f$person_id)[i]
  temp <- dm_sim_f[which(dm_sim_f$person_id == rp_id), ]
  temp_m <- dm_sim_cumsum[which(dm_sim_cumsum$person_id == rp_id), ]
  
  lines(temp$age, temp_m$moves_cumsum, lwd = 1, col = dm_sim_person[which(dm_sim_person$person_id == rp_id),"col_real"])
}

mtext('Plot A', side=3, line=0.5, at=90, font = 2, cex = 1.5)

#legend
points(x = rep(100, times = 130), y = seq(from = 100, to = 150, length.out = 130), col = plot_cols_f)
text(x = 80, y = 100, "few moves", adj = 0, cex = 0.9)     # update y coordinate
text(x = 80, y = 150, "many moves", adj = 0, cex = 0.9)    # update y coordinate   


plot(dm_sim_cumsum$moves_cumsum ~ dm_sim$age ,
     type = "n", 
     xlab  = "Age", 
     ylab = "Number of moves", 
     main = "Cumulative totals of moves \nfor male RPs",
     font.main = 1,
     bty = "n",
     ylim = c(0,150),
     xlim = c(0,100))

# for males
for (i in 1:length(unique(dm_sim_m$person_id))){
  rp_id <- unique(dm_sim_m$person_id)[i]
  temp <- dm_sim_m[which(dm_sim_m$person_id == rp_id), ]
  temp_m <- dm_sim_cumsum[which(dm_sim_cumsum$person_id == rp_id), ]
  
  lines(temp$age, temp_m$moves_cumsum, lwd = 1, col = dm_sim_person[which(dm_sim_person$person_id == rp_id),"col_real"])
}

mtext('Plot B', side=3, line=0.5, at=90, font = 2, cex = 1.5)

#legend
points(x = rep(100, times = 130), y = seq(from = 100, to = 150, length.out = 130), col = plot_cols_m)
text(x = 80, y = 100, "few moves", adj = 0, cex = 0.9)     # update y coordinate
text(x = 80, y = 150, "many moves", adj = 0, cex = 0.9)    # update y coordinate   

dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### sample age structure 

# png("Figures/age_representation.png", res = 300, height = 15, width = 20, units = "cm")
# plot(table(d$age),
#      main = "Representation of ages", 
#      xlab = "Age", 
#      ylab = "Number of observations", 
#      xlim = c(0,100),
#      xaxt = "n",
#      ylim = c(0,40000), 
#      bty = "n", 
#      col = "grey41",
#      font.main = 1)
# axis(1, at = seq(0, 100, by = 10))
# dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### sample cohort representation 

# d %>%
#   group_by(person_id) %>%
#   summarise(b_y = first(b_y)) %>%
#   as.data.frame(stringsAsFactors = FALSE) -> cohorts
# 
# png("Figures/birth_year_rep.png", res = 300, height = 15, width = 20, units = "cm")
# plot(table(cohorts$b_y), 
#      main = "Birth year representation in the HSN",
#      font.main = 1,
#      ylim = c(0,max(table(cohorts$b_y))),
#      xlim = c(1830,1940),
#      xaxt = "n",
#      xlab = "Year",
#      ylab = "Number of RPs",
#      col = "grey41",
#      bty = "n")
# axis(1, at = seq(1840, 1930, 10), labels = seq(1840, 1930, 10))
# dev.off()

#------------------------------------------------------------------------------
# gamma pois model estimates

png("Figures/negbin_estimates_gen.png", type = "cairo", res = 300, height = 6, width = 8, units = "cm")
plot(precis(stanfit_negbin, pars = c("mu", "rho", "eta", "phi"), depth = 2), ylab = "Parameter", main = "Parameter estimates") 
dev.off()


#------------------------------------------------------------------------------
# gamma pois betas

# processing posterior
exp_mu_f <- sapply( 1:uniq_ages , function(i) exp( post_negbin$mu + post_negbin$beta[, 1,i] + rnorm(2000,0,post_negbin$sd_id) ) )
exp_mu_m <- sapply( 1:uniq_ages , function(i) exp( post_negbin$mu + post_negbin$beta[, 2,i] + rnorm(2000,0,post_negbin$sd_id) ) )

beta_mod_int_f <- apply(exp_mu_f, 2, HPDI)
beta_mod_int_m <- apply(exp_mu_m, 2, HPDI)

beta_mod_f <- apply(exp_mu_f, 2, mean)
beta_mod_m <- apply(exp_mu_m, 2, mean)


# plotting
png("Figures/beta_estimates_gen_negbin.png", type = "cairo", res = 300, height = 15, width = 20, units = "cm")

# beta plot
blank()
plot(y = 0:(length(beta_mod_f)-1),
     x = 0:(length(beta_mod_f)-1), 
     ylim = c(0, 0.6),  
     xlim = c(0,(length(beta_mod_f)-1)),
     xlab = "Age", 
     ylab = "Estimated number of moves per year", 
     main = "Age-based estimated moves per year for each sex", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     type = "n")
axis(1, at = seq(0,90, by = 10))
axis(2, at = seq(0, max(beta_mod_int_f), by = 0.1))

#shade( apply(exp_mu_f, 2, PI, 0.75 ) , 0:(length(beta_mod_f)-1), col = cols_int[1])
shade( apply(exp_mu_f, 2, PI, 0.5 ) , 0:(length(beta_mod_f)-1), col = cols_int[1])
#shade( apply(exp_mu_f, 2, PI, 0.25 ) , 0:(length(beta_mod_f)-1), col = cols_int[1])

#shade( apply(exp_mu_m, 2, PI, 0.75 ) , 0:(length(beta_mod_m)-1), col = cols_int[2])
shade( apply(exp_mu_m, 2, PI, 0.5 ) , 0:(length(beta_mod_m)-1), col = cols_int[2])
#shade( apply(exp_mu_m, 2, PI, 0.25 ) , 0:(length(beta_mod_m)-1), col = cols_int[2])

points(y = mean_moves$average_moves, x = 0:(length(mean_moves$average_moves)-1), col = "black", pch = 19, cex = 0.5)
lines(x = 0:(length(beta_mod_f)-1), y = beta_mod_f, col = cols[1], lty = 2, lwd = 3 )
lines(x = 0:(length(beta_mod_m)-1), y = beta_mod_m, col = cols[2], lty = 2, lwd = 3 )

legend(x = 60, y = 0.6, legend = c("female", "male"), lty = 1, col = c(cols[1], cols[2]), cex = 1.5, bty = "n", lwd = 3)
legend(x = 60, y = 0.5, legend = c("data from sample"), pch = 19, col = "black", cex = 1.5, bty = "n")

dev.off()




