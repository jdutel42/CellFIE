###################################################################################"
# ---- Initialization ---- #

#############
# Libraries #
#############

library(shiny)
library(ggplot2)
library(SingleCellExperiment)
library(shinycssloaders)
library(dplyr)
library(shinyalert)
library(shinyFeedback)
library(shinyjs)


###########
# Options #
###########

options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB
useShinyjs()

###################################################################################"
# ---- User Interface ---- #

ui <- fluidPage(
  
  
  ###################
  # ---- Title ---- #
  ###################
  
  titlePanel("Single-cellExperiment FeaturePlot Explorer"),
  
  
  #####################
  # ---- Sidebar ---- #
  #####################
  
  sidebarLayout(
    sidebarPanel(
      
      ########################
      # ---- Input file ---- #
      ########################
      
      fileInput(
        "sce_rds",
        "1. Upload data (SCE object .rds)",
        accept = ".rds"
      ),
      
      ###################
      # ---- Assay ---- #
      ###################
      
      uiOutput("assay_ui"),
      
      #######################
      # ---- Embedding ---- #
      #######################
      
      uiOutput("embedding_ui"),
      
      # Change this part to dynamically get the available reductions from the SCE object !!!!!
      
      # selectInput(
      #   "reduction",
      #   "3. Which embedding ?",
      #   # choices = c("UMAP_uwot", "TSNE", "PCA"),
      #   choices = c("PCA", "TSNE", "UMAP"),
      #   selected = "UMAP"
      # ),
      
      #######################
      # ---- Features ---- #
      #######################
      
      # Change this part to dynamically get the available reductions from the SCE object !!!!!
      
      uiOutput("feature_ui"),
      
      ###################
      # ---- Genes ---- #
      ###################
      
      shinycssloaders::withSpinner(
        uiOutput("gene_ui")
      ),
      
      
      
      #########################
      # ---- Plot Tuning ---- #
      #########################
      
      sliderInput(
        "pt_size",
        "Point size",
        min = 0.1,
        max = 2,
        value = 0.6,
        step = 0.1
      ),
      
      
      ##################################
      # ---- Plot Tuning Metadata ---- #
      ##################################
      
      # shinycssloaders::withSpinner(
      #   uiOutput("metadata_ui")
      # ),
      # 
      # shinycssloaders::withSpinner(
      #   uiOutput("metadata_val_ui")
      # ),
    
      
      
      
      # selectInput(
      #   "metadata",
      #   "Choose a metadata to display",
      #   nested_list_metadata
      # ),
      
      
      ###########################
      # ---- Download Plot ---- #
      ###########################
      
      # Button to download the PCA plot
      downloadButton(
        "download_plot", 
        "Download Plot"
      )
      
      
      
      
      
    ),
    
    ################################
    # ---- Plot visualization ---- #
    ################################
    
    mainPanel(
      shinycssloaders::withSpinner(
        plotOutput(
          "featureplot",
          height = "750px"
        )
      )
    )
  )
)
