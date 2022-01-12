Dutch-historical-mobility
============

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Analysis script for "The complex life course of mobility: Quantitative description of 300,000 residential moves in 1850-1950 Netherlands", by Natalia Fedorova, Richard McElreath, Bret Beheim
DOI: 

# Requirements

- R (3.3.6 or greater) https://cran.r-project.org/
- tidyverse package
- rethinking package (v1.59 or greater), http://xcelab.net/rm/software/
- CmdStanR interface https://mc-stan.org/users/interfaces/cmdstan
- foreign package
- viridis package (for figures)
- testthat package (for reproducibility diagnostics)
- parallel package
- posterior package


# Instructions:

In R, set the working directory to that containing this readme file. For example, on a Mac or Linux machine, you might type into the command prompt

```
  setwd('~/Desktop/Dutch-historical-mobility')
```

if the folder containing the project is named 'Dutch-historical-mobility' and on your Desktop. You can tell if you are in the right place by typing in `dir()` and seeing this readme.md file.

The analysis takes input from the Historical Sample of the Netherlands (https://iisg.amsterdam/en/hsn) through the `1_prep_hsn_data.R` file. However, the default run assumes you need to simulate structurally comparable data.

Run sequence:
- run 2_fit_models.R, this simulates data and runs the poisson regression, beta-gamma regression, and cohort regressions 
- run 3_make_plots.R, this produces the plots present in the paper and supplementary

If you are working with the HSN sample, skip the simulation and uncomment relevant lines

Run sequence:
- move the HSN data to a ./Data_files folder within the repository
- run 1_prep_hsn_data.R, this does all the necessary data cleaning and creates the datasets for analysis
- run 2_fit_models.R, this runs the poisson regression, beta-gamma regression, and cohort regressions 
- run 3_make_plots.R, this produces the plots present in the paper and supplementary

The project is maintained by Natalia Fedorova in a Github repository at https://github.com/Naty-fedorova/Dutch-historical-mobility and licensed under Creative Commons [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). See LICENSE.md for details.
