server <- function(input, output, session) {

  ##########################
  # ---- Set Work Dir ---- #
  ##########################
  
  setwd(dir = "~/Documents/Project/Heimdall/")
  
  #############################
  # ---- Load SCE object ---- #
  #############################
  
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

  observeEvent(sce_obj(), {
    showNotification(
      "Data successfully loaded!",
      type = "message",
      duration = 1
    )
  })
  
  
  
  
  
  
  
  # output$qs_ui <- renderUI({
  #   fileInput(
  #     "qs_file",
  #     "1. Indicate path to data (.qs file)",
  #     accept = c(".qs")
  #   )
  # })
  # 
  # observeEvent(input$qs_file, {
  #   file_path <- input$qs_file$datapath
  #   # Get data from qs file path
  #   data <- qs_read(file_path)
  # })
  
  
  
  # observeEvent(sce_obj(), {
  #   shinyalert(
  #     title = "Information",
  #     text  = "Data successfully loaded!",
  #     type  = "info",
  #     timer = 1000,
  #     showConfirmButton = FALSE
  #   )
  # })
  
  
  # observeEvent(input$sce_rds, {
  #   withProgress(message = "Loading data...", {
  #     sce_obj()
  #   })
  # })
  
  # runjs("toastr.success('Data successfully loaded!')")


  #############################
  # ---- Assay selection ---- #
  #############################
  
  output$assay_ui <- renderUI({
    
    req(input$sce_rds)
    
    selectInput(
      "assay",
      "2. Which assay ?",
      # selected = assayNames(sce_obj())[1],
      choices = assayNames(sce_obj())
    )
  })
  
  # observeEvent(, {
  #   showNotification("Assay successfully loaded!", type = "message")
  # })
  
  
  #################################
  # ---- Embedding selection ---- #
  #################################
  
  output$embedding_ui <- renderUI({
    
    req(input$assay)
    
    selectInput(
      "embedding",
      "3. Which embedding ?",
      # selected = reducedDimNames(sce_obj())[1],
      choices = reducedDimNames(sce_obj())
    )
  })
  
  
  

  ###############################
  # ---- Feature selection ---- #
  ###############################
  
  output$feature_ui <- renderUI({
    
    req(input$embedding)
    
    selectInput(
      "feature",
      "4. Which features ?",
      choices = c("Expression", colnames(colData(sce_obj()))),
      # selected = "Expression"
    )
  })


  ############################
  # ---- Gene selection ---- #
  ############################

  output$gene_ui <- renderUI({

    req(input$feature == "Expression")
    
    genes <- rownames(assay(sce_obj(), input$assay))
    
    selectizeInput(
      "genes",
      "Which Gene(s) ?",
      choices = genes,
      multiple = TRUE,
      options = list(
        maxItems = 4,
        placeholder = "Type gene name(s)..."
      )
    )
  })
  

  
  
  
  ######################
  # ---- Metadata ---- #
  ######################
  
  # # Create a nested list to store all possible metavalues in each col
  output$metadata_ui <- renderUI({
    req(sce_obj())
    
    # Store names of metadata present in the SCE object
    colnames_metadata_list <- colnames(colData(sce_obj()))
    # If col have more than 20 unique value are present, remove them from the choices
    # colnames_metadata_list <- colnames_metadata_list[sapply(colnames_metadata_list, function(col) {
    #   length(unique(colData(sce_obj())[[col]])) <= 30
    # })]
    # colnames_metadata_list <- setdiff(colnames_metadata_list, c("nCounts", "nFeatures"))
    
    selectInput(
      "metadata_col",
      "Choose data to color by",
      choices = colnames_metadata_list
    )
  })
  
  # output$metadata_val_ui <- renderUI({
  #   req(sce_obj(), input$metadata_col)
  #   
  #   vals <- unique(colData(sce_obj())[[input$metadata_col]])
  #   
  #   selectInput(
  #     "metadata_val",
  #     "Choose metadata value",
  #     choices = vals
  #   )
  # })
  
    # selectInput(
    #   "metadata_val",
    #   "Choose a metadata to display",
    #   choices = setNames(
    #     lapply(colnames_metadata_list, function(col) {
    #       unique(colData(sce_obj())[[col]])
    #     }),
    #     colnames_metadata_list
    #   ),
    # )



  
  #########################
  # ---- FeaturePlot ---- #
  #########################
  
    #################################
    # ---- Gather informations ---- #
    #################################
  
  # Create a reactive data.frame for ggplot
  
  df_plot <- reactive({
    sce <- sce_obj()
    
    df_emb <- as.data.frame(reducedDim(sce, input$embedding))
    colnames(df_emb)[1:2] <- c("Dim1", "Dim2")
    df_emb <- df_emb[, 1:2]
    
    df_meta <- as.data.frame(colData(sce))
    
    if (input$feature == "Expression") {
      
      expr_mat <- assay(sce, input$assay)[input$genes, , drop = FALSE] |>
        as.matrix() |>
        t() |>
        as.data.frame()
      
      df <- cbind(df_emb, df_meta, expr_mat)
      
      tidyr::pivot_longer(
        df,
        cols = all_of(input$genes),
        names_to = "gene",
        values_to = "Expression"
      )
      
    } else {
      cbind(df_emb, df_meta)
    }
  })
  


  
  ## Calculate
  featureplot_obj <- reactive({
    df <- df_plot()
    
    if (input$feature == "Expression") {
      ggplot2::ggplot(
        df,
        ggplot2::aes(Dim1, Dim2, color = Expression)
      ) +
        ggplot2::geom_point(size = input$pt_size) +
        ggplot2::scale_color_gradient(low = "lightgrey", high = "blue") +
        ggplot2::facet_wrap(~ gene) +
        ggplot2::theme_classic(base_size = 16)
      
    } else {
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          Dim1,
          Dim2,
          color = .data[[input$feature]]
        )
      ) +
        ggplot2::geom_point(size = input$pt_size) +
        ggplot2::scale_color_brewer(palette = "Set2") +
        ggplot2::theme_classic(base_size = 16)
    }
  })
  
  
  
  
  # featureplot_obj <- reactive({
  #   sce <- sce_obj()
  #   
    
    
    #   ggplot2::ggplot(
    #     df_long,
    #     ggplot2::aes(x = Dim1, y = Dim2, color = Expression)
    #   ) +
    #     ggplot2::geom_point(size = input$pt_size) +
    #     ggplot2::scale_color_gradient(low = "lightgrey", high = "blue") +
    #     ggplot2::facet_wrap(~ gene) +
    #     ggplot2::theme_classic(base_size = 16, base_line_size = 0.5, base_rect_size = 0.5, ink = "black")
    #     # ggplot2::theme_bw()
    # })
    
    ######################
    # ---- Plotting ---- #
    ######################
    
    featureplot_obj <- reactive({
      
      df <- df_plot()
      
      if (input$feature == "Expression") {
        ggplot2::ggplot(
          df,
          ggplot2::aes(Dim1, Dim2, color = Expression)
        ) +
          ggplot2::geom_point(size = input$pt_size) +
          ggplot2::scale_color_gradient(low = "lightgrey", high = "blue") +
          ggplot2::facet_wrap(~ gene) +
          ggplot2::theme_classic(base_size = 16)
        
      } else {
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            Dim1,
            Dim2,
            color = .data[[input$feature]]
          )
        ) +
          ggplot2::geom_point(size = input$pt_size) +
          ggplot2::scale_color_brewer(palette = "Set2") +
          ggplot2::theme_classic(base_size = 16)
      }
    })
    
    
    
    # if (input$feature == "Expression") {
    #   p <- p +
    #     ggplot2::aes(color = Expression) +
    #     ggplot2::scale_color_gradient(low = "lightgrey",
    #                                   high = "blue")
    # } else {
    #   p <- p +
    #     ggplot2::aes(color = .data[[input$feature]]) +
    #     ggplot2::scale_color_brewer(palette = "Set2")
    # }
    
    # p
    
    # if (is.numeric(df_long[[input$feature]])) {
    #   p <- p +
    #     ggplot2::scale_color_gradient(low = "lightgrey",
    #                                   high = "blue")
    #   
    # } else {
    #   p <- p +
    #     ggplot2::scale_color_brewer(palette = "Set2")

      # ## Optionnel : adapter l’échelle selon le type
      # if (is.numeric(df_long[[input$metadata_col]])) {
      #   p <- p + ggplot2::scale_color_viridis_c()
      # } else {
      #   p <- p + ggplot2::scale_color_brewer(palette = "Set2")
      # }
    # }
    # if (is.null(input$metadata_col) || input$metadata_col == "") {
    #   
    #   ## Par défaut : expression
    #   p <- p +
    #     ggplot2::aes(color = Expression) +
    #     ggplot2::scale_color_gradient(
    #       low = "lightgrey",
    #       high = "blue"
    #     )
    #   
    # } else {
    #   
    #   ## Metadata sélectionnée
    #   p <- p +
    #     ggplot2::aes(color = .data[[input$metadata_col]])
    #   
    #   # ## Optionnel : adapter l’échelle selon le type
    #   # if (is.numeric(df_long[[input$metadata_col]])) {
    #   #   p <- p + ggplot2::scale_color_viridis_c()
    #   # } else {
    #   #   p <- p + ggplot2::scale_color_brewer(palette = "Set2")
    #   # }
    # }
    # 
    # p
  # })
  
  ## Ploting featureplot
  output$featureplot <- renderPlot({
    
    req(sce_obj(), input$genes, input$assay, input$embedding)
    
    validate(
      need(
        input$embedding %in% reducedDimNames(sce_obj()),
        paste(
          "embedding",
          input$embedding,
          "not found in object. Ask bioinfo to compute it."
        )
      )
    )
    
    featureplot_obj()
    
  })
  

  ###########################
  # ---- Download plot ---- #
  ###########################
  
  output$download_plot_png <- downloadHandler(
    filename = function() {
      paste("FeaturePlot_SCE_", Sys.Date(), "_", format(Sys.time(), "%X"), ".png", sep = "") # To improve to dynamically change the name (eg gene name or whatever)
    },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 150)
      print(featureplot_obj())
      dev.off()
    }
  )
  
  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste("FeaturePlot_SCE_", Sys.Date(), "_", format(Sys.time(), "%X"), ".pdf", sep = "") # To improve to dynamically change the name (eg gene name or whatever)
    },
    content = function(file) {
      pdf(file, width = 1200, height = 800)
      print(featureplot_obj())
      dev.off()
    }
  )
  
}
