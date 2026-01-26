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
# BiocManager::install("shinyjs")
# BiocManager::install("shinyFeedback")
# BiocManager::install("shinyalert")
library(shinyalert)
library(shinyFeedback)
library(shinyjs)
# BiocManager::install("qs")
library(qs2)


###########
# Options #
###########

options(shiny.maxRequestSize = 5 * 1024^3)  # 5 Go
useShinyjs()

###################################################################################"
# ---- User Interface ---- #

ui <- fluidPage(
  
  ###################
  # ---- Title ---- #
  ###################
  
  titlePanel("CellFIE :Cells & Features Interactive Explorer"),
  
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
        "1. Upload SingleCellExperiment RDS file",
        accept = c(".rds", ".RDS")
      ),
      
      # uiOutput("qs_ui"),
      
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
        "download_plot_png", 
        "Download PNG"
      ),
      
      downloadButton(
        "download_plot_pdf", 
        "Download PDF"
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








#########################################################################################

# How to access the Shiny interface:
#
# 1. On the machine hosting the Shiny application, run the following command in R:
#
#    runApp(
#      appDir = "~/Documents/Project/Heimdall",
#      host = "0.0.0.0",
#      port = 3838,
#      launch.browser = FALSE
#    )
#
# 2. From any other machine on the same network, open a web browser and navigate to:
#
#    http://10.31.208.117:3838
#
#    (Note: the IP address corresponds to the current host machine and may change.)

#########################################################################################