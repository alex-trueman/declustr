library(tidyverse)

samples <- read_csv('data-raw/samples.csv')
mask <- read_csv('data-raw/mask.csv')

devtools::use_data(samples, overwrite = TRUE, compress = 'xz')
devtools::use_data(mask, overwrite = TRUE, compress = 'xz')
