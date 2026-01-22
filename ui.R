library(shiny)
library(ggplot2)
library(SingleCellExperiment)
library(shinycssloaders)

options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

# =======================
# UI
# =======================
ui <- fluidPage(
    
    titlePanel("Single-cellExperiment FeaturePlot Explorer"),
    
    sidebarLayout(
        sidebarPanel(
            fileInput(
                "sce_rds",
                "Upload SCE object (.rds)",
                accept = ".rds"
            ),
            
            shinycssloaders::withSpinner(uiOutput("assay_ui")),
            
            shinycssloaders::withSpinner(uiOutput("gene_ui")),
            
            selectInput(
                "reduction",
                "Embedding",
                # choices = c("UMAP_uwot", "TSNE", "PCA"),
                choices = c("PCA", "TSNE", "UMAP"),
                selected = "UMAP"
            ),
            
            sliderInput(
                "pt_size",
                "Point size",
                min = 0.1,
                max = 2,
                value = 0.6,
                step = 0.1
            ),
            
            # selectInput(
            #   "metadata",
            #   "Choose a metadata to display",
            #   nested_list_metadata
            # ),

            # Button to download the PCA plot
            downloadButton("download_plot", "Download Plot")
        ),
        
        mainPanel(
          shinycssloaders::withSpinner(
            plotOutput("featureplot", height = "750px")
          )
        )
    )
)
