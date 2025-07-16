library(readr)
library(dplyr)
library(purrr)
options(readr.show_col_types = FALSE)
options(readr.show_progress = FALSE)

d <- read_csv(file.path("data", "Master dataset_FuSED.csv")) |>
  select("FW_name", "ecosystem.type", "study_ID", "lat") |>
  left_join(
    read_csv(file.path("data", "ndvi.csv")),
    by = "FW_name"
  ) |>
  full_join(
    read_csv(file.path("data", "chl_a.csv")),
    by = "FW_name"
  ) |>
  mutate(
    metric = ifelse(is.na(ndvi), "Chl-a", "NDVI"),
    avg = ifelse(is.na(ndvi), chl_a, ndvi)
  ) |>
  select(-"ndvi", -"chl_a")

message(
  " - Range NDVI: (",
  d |> 
    filter(metric == "NDVI") |> 
    pull(avg) |> 
    range(na.rm = TRUE) |> 
    round(3) |> 
    paste(collapse = " , ") |> 
    paste0(")")
)

d |> 
  filter(metric == "NDVI") |> 
  pull(avg) |> 
  summary() |> 
  round(3)

message("")

message(
  " - Range Chlorophyll-a: (",
  d |> 
    filter(metric == "Chl-a") |> 
    pull(avg) |> 
    range(na.rm = TRUE) |> 
    round(3) |> 
    paste(collapse = " , ") |> 
    paste0(")")
)

d |> 
  filter(metric == "Chl-a") |> 
  pull(avg) |> 
  summary() |> 
  round(3)

d |> write_csv(file.path("data", "proxy-npp.csv"))
