#!/usr/bin/env Rscript

if (!requireNamespace("shiny", quietly = TRUE)) {
	install.packages("shiny")
}
if (!requireNamespace("factoextra", quietly = TRUE)) {
	install.packages("factoextra")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
	install.packages("tidyverse")
}
if (!requireNamespace("readr", quietly = TRUE)) {
	install.packages("readr")
}
if (!requireNamespace("MVA", quietly = TRUE)) {
	install.packages("MVA")
}
if (!requireNamespace("randomForest", quietly = TRUE)) {
	install.packages("randomForest")
}
if (!requireNamespace("rfUtilities", quietly = TRUE)) {
	install.packages("rfUtilities")
}
if (!requireNamespace("glmnet", quietly = TRUE)) {
	install.packages("glmnet")
}
if (!requireNamespace("cluster", quietly = TRUE)) {
	install.packages("cluster")
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
	install.packages("fgsea")
}
if (!requireNamespace("survminer", quietly = TRUE)) {
	install.packages("survminer")
}
if (!requireNamespace("survival", quietly = TRUE)) {
	install.packages("survival")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
	install.packages("pheatmap")
}
if (!requireNamespace("ggfortify", quietly = TRUE)) {
	install.packages("ggfortify")
}
if (!requireNamespace("plotly", quietly = TRUE)) {
	install.packages("plotly")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
	install.packages("viridis")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
	install.packages("reshape2")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
	install.packages("remotes")
}
if (!requireNamespace("enrichR", quietly = TRUE)) {
	remotes::install_github("ycl6/enrichR@bugfix_2024")
}


