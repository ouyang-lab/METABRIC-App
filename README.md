# Getting Started:
## Installation
### Dependencies
This software uses the Shiny framework. We recommend using R with RStudio. 
The following R package dependencies are required for running this app (`Rscript ./install_dependencies.R` or open in RStudio and select the "Run" button):
 - `shiny`
 - `readr`
 - `MVA`
 - `randomForest`
 - `rfUtilities`
 - `glmnet`
 - `cluster`
 - `fgsea`
 - `factoextra`
 - `tidyverse`
 - `survminer`
 - `survival`
 - `pheatmap`
 - `ggfortify`
 - `plotly`
 - `viridis`
 - `enrichR` : `remotes::install_github("ycl6/enrichR@bugfix_2024")`
 - `reshape2`

## Running the METABRIC App
Open the `app.R` file in RStudio, select the "Run App" button, and the app will load in a new window.

## App features: 
For each gene in the dataset:
 - Unsupervised random forest (URF) importance rank
 - Kaplan-Meier survival analysis across clusters
 - Cellularity anylsis across quantiles
 - HER2 status across quantiles
 - ER status across quantiles
 - Cancer types across quantiles
 - Gene set enrichment (Reactome 2022)
 - Gene expression correlations


