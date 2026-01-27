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
library(bslib)


###########
# Options #
###########

options(shiny.maxRequestSize = 5 * 1024^3)  # 5 Go
useShinyjs()

###################################################################################"
# ---- User Interface ---- #

ui <- page_fillable(
  
  ###################
  # ---- Title ---- #
  ###################
  
  titlePanel("CellFIE : Cells & Features Interactive Explorer"),
  
  
  ####################
  # ---- NavTab ---- #
  ####################
  
  navset_card_tab(
    
    navset_card_tab(
      
      # =========================
      # ---- Home / Welcome ---- #
      # =========================
      
      nav_panel(
        "Home",
        
        layout_column_wrap(
          width = 1/2,
          
          # ---- Welcome & Overview ----
          card(
            card_header("Welcome to CellFIE"),
            p(
              "CellFIE (Cells & Features Interactive Explorer) is an interactive Shiny application ",
              "designed for the exploration and visualization of single-cell data ",
              "based on SingleCellExperiment objects."
            ),
            tags$ul(
              tags$li("Load and inspect SingleCellExperiment datasets"),
              tags$li("Explore embeddings and features interactively"),
              tags$li("Visualize gene expression and metadata"),
              tags$li("Export publication-ready figures")
            ),
            hr(),
            h5("Application workflow"),
            tags$ol(
              tags$li("Upload your SingleCellExperiment object"),
              tags$li("Select assay and dimensional reduction"),
              tags$li("Explore features, genes and cell populations"),
              tags$li("Download figures")
            )
          ),
          
          # ---- Data Loading ----
          card(
            card_header("Load your data"),
            
            ########################
            # ---- Input file ---- #
            ########################
            
            fileInput(
              "sce_rds",
              "1. Upload SingleCellExperiment RDS file",
              accept = c(".rds", ".RDS")
            ),
            
            # uiOutput("qs_ui"),  
            
            p(
              "The uploaded object must be a valid ",
              tags$code("SingleCellExperiment"),
              " containing assays, reduced dimensions, and metadata."
            ),
            
            # shinyFeedback::feedbackDanger(
            #   "sce_rds",
            #   show = FALSE,
            #   text = "Invalid or unsupported object."
            # ),
            
          # ---- Navigation / Index ----
          card(
            card_header("Explore the application"),
            
            p("Once your data is loaded, navigate to the following sections:"),
            
            tags$ul(
              tags$li(
                tags$b("FeaturePlot"), ": visualize gene or feature expression on embeddings"
              ),
              tags$li(
                tags$b("Cell Explorer"), ": explore cell-level metadata and annotations"
              )
            ),
            
            p(
              "Use the navigation tabs above to access each module."
            )
          )
        )
      ),
      
    #   # ======================
    #   # ---- FeaturePlot ---- #
    #   # ======================
    #   
    #   nav_panel(
    #     "FeaturePlot",
    #     "FeaturePlot UI here"
    #   ),
    #   
    #   # =========================
    #   # ---- Cell Explorer ---- #
    #   # =========================
    #   
    #   nav_panel(
    #     "Cell Explorer",
    #     "Cell Explorer UI here"
    #   )
    ),
    
    #############################
    # ---- FeaturePlot Tab ---- #
    #############################
    
    nav_panel("FeaturePlot", 
    
      #####################
      # ---- Sidebar ---- #
      #####################
      
      sidebarLayout(
        sidebarPanel(
          
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
    ),
    
    nav_panel("VolcanoPlot",
              "Work in progress..."
    ),
    
    nav_panel("About",
              ################
              # ---- About ---- #
              ################
              
              "Some information about the application.",
              # includeMarkdown("www/About.md")
    )
              
  ),
  

  

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