# BEF-in-Food-Webs
Code and data for to reproducing all results in Barnes *et al.* "Food Web Complexity Underlies the Relationship Between Biodiversity and Ecosystem Functioning". 

## Data
A) Datasets with node-level information used to calculate energy fluxes and food web metrics. Each folder contains (i) one adjacency matrix (a resource by consumer matrix of trophic interactions) and (ii) one taxon attributes file (containing information such as body mass and biomass of individual taxa) per food web, as well as (iii) a single 'meta' file containing food web-level information for all food webs contained within each folder. Data are contained in the following folders:
  1. Biodiversity Exploratories
  2. Brazilian streams_Saito
  3. IcelandicStreams
  4. Intertidalrockpools
  5. Lakes
  6. Rivers_Mario
  7. Russian soils
  8. UKstreams
  9. baltic_sea

B) Derived datasets with food web-level information used to carry out statistical analyses and produce results in Barnes *et al.* (2026). Datasets are divided by ecosystem type, accompanied by a metadata file with information on each data column, as follows:
  1. Barnes et al Metadata.xlsx
  2. meta.Lakes.csv
  3. meta.Marine_correction.csv
  4. meta.Soils_correction.csv
  5. meta.Streams_correction.csv

## Code
Custom R code to calculate energy fluxes and food web metrics to reproduce results presented in the forthcoming author correction of Barnes *et al.*. Original results presented in Barnes *et al.* (2026) can be reproduced using the code in the folder "Original code for Barnes et al 2026".

The following R scripts call datasets contained within the node-level data folder structure (A) to calculate taxa richness, energy fluxes (using the 'fluxweb' R package), and food web metrics found in the derived datasets (B): '1. FuSED Lakes_correction.r', '2. FuSED Marine_correction.r', '3. FuSED Soils_correction.r', and '4. FuSED Streams_correction.r'. 

Each of these scripts additionally call functions from the file 'Food_web_functions.r' created by Benoit Gauzens (https://github.com/gauzens), which contains custom functions to calculate various food web properties.

The folder 'NDVI and Chlorophyll-a' provides code and data for estimation of NPP via remotely sensed NDVI & Chlorophyll-a. We used the Normalized Difference Vegetation Index (NDVI) and the concentration of chlorophyll-a as proxies for the net primary productivity (NPP) for land and marine environments, respectively.

The folder 'sensitivities' provides code used to evaluate how uncertainty in metabolic rates and assimilation efficiencies influences energy‑flux estimates and the resulting structural equation models. By repeatedly resampling published regression coefficients, it generates 1000 alternative model outputs and corresponding distributions of p‑values and R² values.








