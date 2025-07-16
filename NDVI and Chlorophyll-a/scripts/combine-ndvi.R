library(readr)
library(dplyr)
library(purrr)
options(readr.show_col_types = FALSE)
options(readr.show_progress = FALSE)

ff <- list.files(file.path("data", "ndvi"), full.names = TRUE)

# load and merge NDVI files -----------
d <- lapply(
  ff,
  \(f) {
    fw <- strsplit(f, "/")[[1]]
    fw <- fw[length(fw)]
    fw <- gsub("[.]csv", "", fw)
    f |>
      read_csv() |>
      mutate(
        FW_name = fw,
        ndvi = mean,
        .keep = "none"
      )
  }
)
d <- d[sapply(d, nrow) > 0] # Constable lake + Baltic sea
d <- bind_rows(d)

all_fw <- file.path("data", "Master dataset_FuSED.csv") |>
  read_csv() |> 
  filter(
    study_ID != "Intertidal rockpools",
    ecosystem.type != "Marine"
  ) |>
  pull(FW_name)

missing <- setdiff(all_fw, d |> pull(FW_name))

d <- bind_rows(
  d,
  tibble(FW_name = missing, ndvi = NA)
)

d |> write_csv(file.path("data", "ndvi.csv"))
