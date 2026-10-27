# Load project functions and both base datasets with the smoke configuration.
for (f in list.files(test_path("..", "..", "R"), full.names = TRUE)) source(f)
library(data.table)
cfg_test <- config::get(file = test_path("..", "..", "config", "config.yml"), config = "smoke")
spec_m <- cfg_test$datasets$maringe
spec_s <- cfg_test$datasets$synpuf
spec_d <- cfg_test$datasets$maringe_surv12
base_test <- prepare_base_data(test_path("..", "..", spec_m$file), spec_m, "maringe")
base_syn <- prepare_base_data(test_path("..", "..", spec_s$file), spec_s, "synpuf")
