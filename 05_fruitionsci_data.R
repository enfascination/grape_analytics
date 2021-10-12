library(tidyverse)

dataPath <- "../data/forrestel/FruitionSciences/"
temp = list.files(dataPath, pattern="*.csv")
tbls <- list()

if (TRUE) {
  for (i in 1:length(temp)) {
    tbls[[temp[i]]] <- read_csv(paste0(dataPath, temp[i]))
  }
}

col_names <- c()
for (i in 1:length(tbls)) {
  col_names <- c(names(tbls[[i]]), col_names)
}
col_names <- unique(col_names)
print(col_names)
