#####libraries needed:
is_inst <- function(pkg) {
  nzchar(system.file(package = pkg))
}

if (!is_inst("pacman")) {install.packages("pacman")}
library("pacman")

p_load(tidyverse,pbapply,optimx,purrr,matrixStats,pracma,magrittr)
