server <- function(input, output, session) {
  
  # ---- Load SCE object ----
  sce_obj <- reactive({
    req(input$sce_rds)
    sce <- readRDS(input$sce_rds$datapath)
    
    validate(
      need(
        inherits(sce, "SingleCellExperiment"),
        "The uploaded file is not a SingleCellExperiment object"
      )
    )
    
    sce
  })
  
  # ---- Assay selection ----
  output$assay_ui <- renderUI({
    req(sce_obj())
    selectInput(
      "assay",
      "Assay",
      choices = assayNames(sce_obj()),
      selected = assayNames(sce_obj())[1]
    )
  })
  
  
  # ---- Gene selection ----
  output$gene_ui <- renderUI({
    req(sce_obj(), input$assay)
    
    genes <- rownames(assay(sce_obj(), input$assay))
    
    selectizeInput(
      "genes",
      "Gene(s)",
      choices = genes,
      multiple = TRUE,
      options = list(
        maxItems = 4,
        placeholder = "Type gene name(s)..."
      )
    )
  })
  
  
  
  # ---- FeaturePlot ----
  ## Calculate
  featureplot_obj <- reactive({
    sce <- sce_obj()
    
    validate(
      need(
        input$reduction %in% reducedDimNames(sce),
        paste(
          "Reduction",
          input$reduction,
          "not found in object. Ask bioinfo to compute it."
        )
      )
    )
    
    # print("Ok1")
    
    # Coordinates
    emb <- as.data.frame(reducedDim(sce, input$reduction))
    colnames(emb)[1:2] <- c("Dim1", "Dim2")
    
    # print("Ok2")
    
    # Expression
    expr <- assay(sce, input$assay)[input$genes, , drop = FALSE]
    
    # print("Ok3")
    
    # Convert sparse matrix -> dense matrix
    expr_mat <- as.matrix(expr)
    
    # print("Ok3.1")
    
    # Combine
    df <- cbind(
      as.data.frame(emb),
      as.data.frame(t(expr_mat))
    )
    
    # print("Ok4")
    
    df_long <- tidyr::pivot_longer(
      df,
      cols = all_of(input$genes),
      names_to = "gene",
      values_to = "Expression"
    )
    
    # print("Ok5")
    
    ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = Dim1, y = Dim2, color = Expression)
    ) +
      ggplot2::geom_point(size = input$pt_size) +
      ggplot2::scale_color_gradient(low = "lightgrey", high = "blue") +
      ggplot2::facet_wrap(~ gene) +
      ggplot2::theme_classic(base_size = 16, base_line_size = 0.5, base_rect_size = 0.5, ink = "black")
      # ggplot2::theme_bw()
  })
  
  
  
  
  ## Ploting featureplot
  output$featureplot <- renderPlot({
    req(sce_obj(), input$genes, input$assay, input$reduction)
  
    
    
    featureplot_obj()
    
  })
  
  # ---- Download plot ----
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("FeaturePlot_SCE_", Sys.Date(), "_", format(Sys.time(), "%X"), ".png", sep = "") # To improve to dynamically change the name (eg gene name or whatever)
    },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 150)
      print(featureplot_obj())
      dev.off()
    }
  )
  
  ######################
  # ---- Metadata ---- #
  ######################
  
  # Create a nested list to store all possible metavalues in each col
  nested_list_metadata <- reactive({
    
    # Store names of metadata present in the SCE object
    colnames_metadata_list <- colnames(colData(sce_obj()))

    # Create the list of list with names according colnames
    setNames(
      lapply(colnames_metadata_list, function(col) {
        unique(colData(sce_obj())[[col]])
      }),
      colnames_metadata_list
    )
  })
  
  
  
}
