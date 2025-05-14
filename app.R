#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(MVA)
library(randomForest)
library(rfUtilities)
library(glmnet)
library(cluster)
library(fgsea)
library(factoextra)
library(tidyverse)
library(survminer)
library(survival)
library(pheatmap)
library(ggfortify)
library(plotly)
library(viridis)
library(enrichR)
library(reshape2)

setEnrichrSite("Enrichr")
getOption("enrichR.live")
meta_full <- read_csv("filt_meta_full.csv", show_col_types = FALSE)
var_imp <- read_csv("var_imp.csv", show_col_types = FALSE)
var_imp <- var_imp %>% arrange(desc(MeanDecreaseGini)) %>%
  mutate(rank = rank(-MeanDecreaseGini))
gene_names <- var_imp %>% filter(str_detect(var_name, "quat") == TRUE) %>%
  dplyr::pull(.,labels2)
knn_vars2 <- meta_full %>% mutate(across(c(FOXA1:TCN1),
            ~ ntile(.,4), .names = "{.col}_quat")) %>% na.omit() %>%
   mutate(cancer_type = 
               case_when(
                 HISTOLOGICAL_SUBTYPE %in% c("Medullary", "Mucinous",
          "Other", "Tubular/ cribriform") ~ "Other",
                 TRUE ~ HISTOLOGICAL_SUBTYPE))

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("METABRIC Analysis App"),

    # Sidebar with a slider input for number of bins 
    fluidRow(
      conditionalPanel(condition="input.tabselected==1",
        column(2,
            selectInput("gene", "Gene:",
                        choices = gene_names))),
      conditionalPanel(condition="input.tabselected==3",
                       column(2,
                              selectInput("gene1", "Gene:",
                                          choices = gene_names,
                                          selected = "ESR1"),
                              selectInput("gene2", "Gene:",
                                          choices = gene_names,
                                          selected = "GATA3"),
                              selectInput("gene3", "Gene:",
                                          choices = gene_names,
                                          selected = "CA12"),
                              selectInput("gene4", "Gene:",
                                          choices = gene_names,
                                          selected = "TBC1D9"),
                              selectInput("gene5", "Gene:",
                                          choices = gene_names,
                                          selected = "PSAT1"),
                              selectInput("gene6", "Gene:",
                                          choices = gene_names,
                                          selected = "AGR3"),
                              selectInput("gene7", "Gene:",
                                          choices = gene_names,
                                          selected = "IL6ST"))),
      conditionalPanel(condition="input.tabselected==2",
                       column(2,
              numericInput("obs", "Number of Genes:", 30, min = 5,
                           max = 500),
              selectInput('data_base',"Choose a Database",
                          c("Reactome 2022" = "Reactome_2022",
                            "Kegg 2021" = "KEGG_2021_Human",
                    "WikiPathways 2023" = "WikiPathway_2023_Human",
  "GO Biological Process 2023" = "GO_Biological_Process_2023",
  "GO Molecular Function 2023" = "GO_Molecular_Function_2023",
  "GO Cellular Component 2023" = "GO_Cellular_Component_2023")))),
      mainPanel(
        tabsetPanel(type = "tabs", id = "tabselected",
          tabPanel("Unsupervised Random Forest Rank",
                   value = 1, textOutput("unsup")),
          tabPanel("KM Plot",value = 1,plotOutput("test")),
          tabPanel("Cellularity Plot",value = 1, plotOutput("cell")),
          tabPanel("HER2 Status Plot",value = 1, plotOutput("her2")),
          tabPanel("ER Status Plot",value = 1, plotOutput("er")),
          tabPanel("Cancer Type Plot",value = 1, plotOutput("type")),
        tabPanel("Gene Enrichment Plot", value = 2, plotlyOutput("enrich")),
  tabPanel("Gene Correlation Plot", value = 3, plotlyOutput("corr_plot"))
        )
      )

    )
)

# Define server logic required to draw a histogram
server <- function(input,session, output) {
  default_val <- 30
    output$test <- renderPlot({
      new_name <- paste0(input$gene,"_quat")
      lr_real <- survdiff(Surv(cens_time, cens_stat)
                       ~get(new_name),data = knn_vars2)
      lr_real_pval <-signif(pchisq(lr_real$chisq, length(lr_real$n)-1, lower.tail = FALSE),digits=5)
      #print(lr_real_pval)
      rf_km <- survfit(Surv(cens_time, cens_stat)
                       ~get(new_name),data = knn_vars2)
      ct <- broom::tidy(rf_km)
      ct %>% ggplot(aes(x = time, y = estimate, color = strata)) +
        geom_step(linewidth = 1.2) +
        geom_point(data = ct[ct$n.censor >= 1,],
                   aes(x = time, y = estimate, color = strata),
                   shape = "|", size = 3, stroke = 3) +
        annotate("text", x= 15, y = 0.55,
                 label = "Log-Rank Test", size = 5) +
        annotate("text", x= 15, y = 0.5,
                 #label = paste("P-Value","=", as.character(lr_real_pval)),
                 label = paste("p-value","=",as.character(lr_real_pval)),
                 size = 5) +
        labs(x = 'Time in Months', y = 'Survival Probability',
             color = "Expression Level:",
             title = paste(input$gene, "Expression Levels KM Plot")) +
        scale_y_continuous(limits = c(0.5,1)) +
        scale_x_continuous(breaks = c(0,30,60,90,120)) +
        scale_color_discrete(labels = c("Low",
                                        "Medium-Low","Medium-High","High")) +
        theme_classic(base_size = 15) +
        theme(legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.position = 'top',
              plot.title = element_text(hjust = 0.5))
    })
    output$unsup <- renderText({
      new_name_unsup <- paste0(input$gene,"_quat")
      new_dat <- var_imp %>% filter(var_name == new_name_unsup)
      paste(input$gene,"is the",
      paste0(new_dat$rank,"th"),"most important variable
            in the Unsupervised Random Forest out of 506 variables")
    })
    output$cell <- renderPlot({
      new_name <- paste0(input$gene,"_quat")
      chi <- chisq.test(knn_vars2[[new_name]], knn_vars2[["CELLULARITY"]])
      knn_vars2 %>% group_by(CELLULARITY, get(new_name)) %>%
        summarize(total = n()) %>% rename(., var = `get(new_name)`) %>%
        ggplot(aes(x = CELLULARITY)) +
        geom_bar(aes(y =total,
                     fill = as.factor(var)),
                 stat = "identity", position = "dodge") +
        geom_text(aes(y = total, label = paste("N:", total),
                      group = as.factor(var)),
    position = position_dodge(width = .9), vjust = -0.5, size = 5) +
        annotate("text", x = 2, y = 175,
              label = paste("Chi Square p-value:\n", round(chi$p.value,4))) +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "Cellularity", y = "Total Amount of Patients",
             fill = "Quantile", title = input$gene) +
        theme_classic2(base_size = 18) +
        theme(plot.title = element_text(hjust = 0.5))
    })
    output$her2 <- renderPlot({
      new_name <- paste0(input$gene,"_quat")
      chi <- chisq.test(knn_vars2[[new_name]], knn_vars2[["her2_real"]])
      knn_vars2 %>% group_by(her2_real, get(new_name)) %>%
        summarize(total = n()) %>% rename(., var = `get(new_name)`) %>%
        ggplot(aes(x = her2_real)) +
        geom_bar(aes(y =total,
                     fill = as.factor(var)),
                 stat = "identity", position = "dodge") +
        geom_text(aes(y = total, label = paste("N:", total),
                      group = as.factor(var)),
        position = position_dodge(width = .9), vjust = -0.5, size = 5) +
        annotate("text", x = 2.3, y = 275,
       label = paste("Chi Square P-Value:\n", round(  chi$p.value,4))) +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "HER2 Status", y = "Total Amount of Patients",
             fill = "Quantile", title = input$gene) +
        theme_classic2(base_size = 18) +
        theme(plot.title = element_text(hjust = 0.5))
    })
    output$er <- renderPlot({
      new_name <- paste0(input$gene,"_quat")
      chi <- chisq.test(knn_vars2[[new_name]], knn_vars2[["ER_IHC"]])
      knn_vars2 %>% group_by(ER_IHC, get(new_name)) %>%
        summarize(total = n()) %>% rename(., var = `get(new_name)`) %>%
        ggplot(aes(x = ER_IHC)) +
        geom_bar(aes(y =total,
                     fill = as.factor(var)),
                 stat = "identity", position = "dodge") +
        geom_text(aes(y = total, label = paste("N:", total),
                      group = as.factor(var)),
        position = position_dodge(width = .9), vjust = -0.5, size = 5) +
        annotate("text", x = 0.7, y = 375,
        label = paste("Chi Square P-Value:\n", round(chi$p.value,4))) +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "ER Status", y = "Total Amount of Patients",
             fill = "Quantile", title = input$gene) +
        theme_classic2(base_size = 18) +
        theme(plot.title = element_text(hjust = 0.5))
    })
    output$type <- renderPlot({
      new_name <- paste0(input$gene,"_quat")
      chi <- chisq.test(knn_vars2[[new_name]], knn_vars2[["cancer_type"]])
      knn_vars2 %>% group_by(cancer_type, get(new_name)) %>%
        summarize(total = n()) %>% rename(., var = `get(new_name)`) %>%
        ggplot(aes(x = cancer_type)) +
        geom_bar(aes(y =total,
                     fill = as.factor(var)),
                 stat = "identity", position = "dodge") +
        geom_text(aes(y = total, label = paste("N:", total),
                group = as.factor(var)),
          position = position_dodge(width = .9), vjust = -0.5, size = 5) +
        annotate("text", x = 2.5, y = 274,
        label = paste("Chi Square P-Value:\n", round(chi$p.value,4))) +
        scale_fill_viridis(discrete = TRUE) +
        labs(x = "ER Status", y = "Total Amount of Patients",
             fill = "Quantile", title = input$gene) +
        theme_classic2(base_size = 18) +
        theme(plot.title = element_text(hjust = 0.5))
    })
    output$enrich <- renderPlotly({
      top_genes_rf <- var_imp %>% dplyr::arrange(desc(MeanDecreaseGini)) %>%
        filter(str_detect(var_name, "quat")) %>%
        slice(1:input$obs) %>% pull(., var_name)
      top_genes_rf <- stringr::word(top_genes_rf, 1, sep="_")
      
      result <- enrichr(top_genes_rf, input$data_base)
      dat <- as.data.frame(result[[input$data_base]])
      dat$total <- gsub("/.*$","",dat$Overlap)
      dat$term_real <- sub("\\s*\\(.*", "", dat$Term)
      dat$gene_labs <- gsub(";", " ", dat$Genes, fixed=TRUE)
  ggplotly(dat %>% mutate(total = as.numeric(total)) %>%
             arrange(Adjusted.P.value) %>% slice(1:10) %>%
        ggplot() + geom_bar(aes(x=fct_reorder(term_real, Adjusted.P.value,
                                              .desc = TRUE),
                                y=total,fill=desc(Adjusted.P.value),
                                text = paste("Pathway:",term_real,
                                             "<br>Genes:",gene_labs,
                                             "<br> Adjusted P-Value:",
                                             round(Adjusted.P.value,4))),
                            stat="identity") + 
  scale_x_discrete(labels = function(x) str_wrap(substr(x, 1, 35),25)) +
    scale_y_continuous(name ="Total Number of Genes") +
        scale_fill_viridis_c(option = "D") +
        labs(x = 'Pathway', fill = 'Adjusted p-value') +
        theme_classic() +
        theme (axis.text.y = element_text(size=8.5),
               axis.text.x = element_text(size=12)) + coord_flip(),
  tooltip = "text")
      
    })
    output$corr_plot <- renderPlotly({
all_genes <- c(input$gene1,input$gene2,input$gene3,input$gene4,input$gene5,
                  input$gene6,input$gene7)
      top_genes_rf <- var_imp %>% dplyr::arrange(desc(MeanDecreaseGini)) %>%
        filter(str_detect(var_name, "quat")) %>%
        filter(labels2 %in% all_genes) %>% pull(., var_name)
      top_genes_rf <- stringr::word(top_genes_rf, 1, sep="_")
      
      alt_heatmap_vars <- knn_vars2 %>%
        dplyr::select(all_of(top_genes_rf)) 
      corr_var <- cor(alt_heatmap_vars)
      #p_val <- ggcorrplot::cor_pmat(alt_heatmap_vars)
      melt_var <- melt(corr_var)
      p <- melt_var %>% ggplot(aes(x = Var1, y = Var2, fill = value,
                                   text = paste("Gene 1:", Var1,
                                                "<br>Gene 2:", Var2,
                                                "<br>Corr:",
                                                round(value,4)))) +
        geom_tile() +
  scale_fill_gradient2(low = "#fcfdbf", high = "#000004", mid = "#b73779", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_classic2(base_size = 10) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(x = 'Gene 1', y = "Gene 2")
      ggplotly(p,tooltip = "text")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
