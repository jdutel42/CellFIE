####################################################################################
# ---- CellFIE server.R — Version 2.0 ---- #
####################################################################################

server <- function(input, output, session) {

  ##############################################################################
  # ---- Helpers ---- #
  ##############################################################################

  # Retourne le thème ggplot selon le choix utilisateur
  get_gg_theme <- function(theme_id, base_size = 14) {
    switch(theme_id,
      "classic" = ggplot2::theme_classic(base_size = base_size),
      "minimal"  = ggplot2::theme_minimal(base_size = base_size),
      "bw"       = ggplot2::theme_bw(base_size = base_size),
      "void"     = ggplot2::theme_void(base_size = base_size),
      "light"    = ggplot2::theme_light(base_size = base_size),
      ggplot2::theme_classic(base_size = base_size)
    )
  }

  # Sauvegarde générique PNG / PDF / SVG
  make_download_handler <- function(plot_reactive, prefix_input, w_input, h_input) {
    list(
      png = downloadHandler(
        filename = function() paste0(prefix_input, "_", Sys.Date(), ".png"),
        content  = function(file) {
          w <- if (!is.null(w_input())) w_input() else 1200
          h <- if (!is.null(h_input())) h_input() else 900
          png(file, width = w, height = h, res = 150)
          print(plot_reactive())
          dev.off()
        }
      ),
      pdf = downloadHandler(
        filename = function() paste0(prefix_input, "_", Sys.Date(), ".pdf"),
        content  = function(file) {
          w <- if (!is.null(w_input())) w_input() / 72 else 14
          h <- if (!is.null(h_input())) h_input() / 72 else 10
          pdf(file, width = w, height = h)
          print(plot_reactive())
          dev.off()
        }
      ),
      svg = downloadHandler(
        filename = function() paste0(prefix_input, "_", Sys.Date(), ".svg"),
        content  = function(file) {
          w <- if (!is.null(w_input())) w_input() / 96 else 12
          h <- if (!is.null(h_input())) h_input() / 96 else 9
          svg(file, width = w, height = h)
          print(plot_reactive())
          dev.off()
        }
      )
    )
  }

  ##############################################################################
  # ---- Chargement des données ---- #
  ##############################################################################

  output$qs_ui <- renderUI({
    fileInput(
      "data_input",
      "Sélectionnez votre fichier",
      accept = c(".rds", ".RDS", ".qs", ".QS"),
      buttonLabel = "Parcourir...",
      placeholder = "Aucun fichier sélectionné"
    )
  })

  sce_obj <- reactive({
    req(input$data_input)

    ext <- tolower(tools::file_ext(input$data_input$name))

    withProgress(message = "Chargement des données...", value = 0.3, {
      obj <- tryCatch({
        if (ext == "qs") {
          qs2::qread(input$data_input$datapath)
        } else if (ext == "rds") {
          readRDS(input$data_input$datapath)
        } else {
          showNotification("Format non supporté (.rds ou .qs uniquement)", type = "error", duration = 5)
          return(NULL)
        }
      }, error = function(e) {
        showNotification(paste("Erreur de lecture :", e$message), type = "error", duration = 8)
        return(NULL)
      })

      setProgress(0.7, message = "Conversion de l'objet...")

      if (inherits(obj, "Seurat")) {
        message("Objet Seurat détecté → conversion en SCE")
        sce <- Seurat::as.SingleCellExperiment(obj)
      } else if (inherits(obj, "SingleCellExperiment")) {
        sce <- obj
      } else {
        showNotification("L'objet doit être un Seurat ou SingleCellExperiment", type = "error")
        return(NULL)
      }

      setProgress(1, message = "Prêt !")
      sce
    })
  })

  observeEvent(sce_obj(), {
    showNotification(
      paste0("✅ Données chargées : ", ncol(sce_obj()), " cellules, ", nrow(sce_obj()), " gènes"),
      type = "message", duration = 3
    )
  })

  # ---- Badge statut dans la navbar ----
  output$data_status_badge <- renderUI({
    if (is.null(sce_obj())) {
      tags$span("⚠️ Aucune donnée", class = "badge bg-warning text-dark")
    } else {
      tags$span(
        paste0("✅ ", ncol(sce_obj()), " cellules"),
        class = "badge bg-success"
      )
    }
  })

  # ---- Résumé data (Home) ----
  output$data_summary_ui <- renderUI({
    req(sce_obj())
    sce <- sce_obj()
    tagList(
      hr(),
      tags$b("✅ Objet chargé avec succès"),
      br(), br(),
      fluidRow(
        column(6, div(class = "stat-box",
          div(class = "stat-val", format(ncol(sce), big.mark = " ")),
          div(class = "stat-lab", "Cellules")
        )),
        column(6, div(class = "stat-box",
          div(class = "stat-val", format(nrow(sce), big.mark = " ")),
          div(class = "stat-lab", "Gènes")
        ))
      ),
      br(),
      fluidRow(
        column(6, div(class = "stat-box",
          div(class = "stat-val", length(assayNames(sce))),
          div(class = "stat-lab", "Assays")
        )),
        column(6, div(class = "stat-box",
          div(class = "stat-val", length(reducedDimNames(sce))),
          div(class = "stat-lab", "Réductions")
        ))
      )
    )
  })

  output$stats_boxes_ui <- renderUI({
    req(sce_obj())
    sce <- sce_obj()
    tagList(
      fluidRow(
        column(3, div(class = "stat-box",
          div(class = "stat-val", format(ncol(sce), big.mark = " ")),
          div(class = "stat-lab", "Cellules")
        )),
        column(3, div(class = "stat-box",
          div(class = "stat-val", format(nrow(sce), big.mark = " ")),
          div(class = "stat-lab", "Gènes")
        )),
        column(3, div(class = "stat-box",
          div(class = "stat-val", length(assayNames(sce))),
          div(class = "stat-lab", "Assays")
        )),
        column(3, div(class = "stat-box",
          div(class = "stat-val", length(reducedDimNames(sce))),
          div(class = "stat-lab", "Réductions dim.")
        ))
      )
    )
  })

  ##############################################################################
  # ---- Sélecteurs partagés ---- #
  ##############################################################################

  output$assay_ui <- renderUI({
    req(sce_obj())
    selectInput("assay", "Assay", choices = assayNames(sce_obj()), selected = assayNames(sce_obj())[1])
  })

  output$embedding_ui <- renderUI({
    req(sce_obj())
    rdnames <- reducedDimNames(sce_obj())
    # Priorité UMAP > TSNE > PCA
    default <- if ("UMAP" %in% rdnames) "UMAP" else if ("TSNE" %in% rdnames) "TSNE" else rdnames[1]
    selectInput("embedding", "Réduction dimensionnelle", choices = rdnames, selected = default)
  })

  output$feature_ui <- renderUI({
    req(sce_obj())
    selectInput("feature", "Afficher",
                choices = c("Expression génique" = "Expression", colnames(colData(sce_obj()))),
                selected = "Expression")
  })

  output$gene_ui <- renderUI({
    req(input$feature == "Expression", sce_obj(), input$assay)
    genes <- rownames(assay(sce_obj(), input$assay))
    selectizeInput("genes", "Gène(s) à afficher",
                   choices = genes, multiple = TRUE,
                   options = list(maxItems = 6, placeholder = "Tapez un nom de gène..."))
  })

  ##############################################################################
  # ---- FeaturePlot ---- #
  ##############################################################################

  output$fp_card_header <- renderUI({
    req(sce_obj())
    if (!is.null(input$genes) && length(input$genes) > 0) {
      paste("Expression :", paste(input$genes, collapse = ", "), "—", input$embedding)
    } else if (!is.null(input$feature) && input$feature != "Expression") {
      paste(input$feature, "—", input$embedding)
    } else {
      "FeaturePlot"
    }
  })

  df_featureplot <- reactive({
    req(sce_obj(), input$embedding, input$feature)
    sce <- sce_obj()

    df_emb <- as.data.frame(reducedDim(sce, input$embedding))[, 1:2]
    colnames(df_emb) <- c("Dim1", "Dim2")
    df_meta <- as.data.frame(colData(sce))

    if (input$feature == "Expression") {
      req(input$genes, length(input$genes) > 0)
      expr_mat <- assay(sce, input$assay)[input$genes, , drop = FALSE] |>
        as.matrix() |> t() |> as.data.frame()
      df <- cbind(df_emb, df_meta, expr_mat)
      tidyr::pivot_longer(df, cols = all_of(input$genes),
                          names_to = "gene", values_to = "Expression")
    } else {
      cbind(df_emb, df_meta)
    }
  })

  featureplot_obj <- reactive({
    req(df_featureplot())
    df  <- df_featureplot()
    th  <- get_gg_theme(req(input$fp_theme), req(input$fp_base_size))
    ttl <- if (!is.null(input$fp_title) && nchar(input$fp_title) > 0) input$fp_title else NULL
    show_legend <- isTRUE(input$fp_legend)

    if (input$feature == "Expression") {
      p <- ggplot2::ggplot(df, ggplot2::aes(Dim1, Dim2, color = Expression)) +
        ggplot2::geom_point(size = input$fp_pt_size, alpha = input$fp_pt_alpha) +
        ggplot2::scale_color_gradient(low = input$fp_col_low, high = input$fp_col_high) +
        ggplot2::facet_wrap(~ gene) +
        ggplot2::labs(title = ttl,
                      x = paste0(input$embedding, "_1"),
                      y = paste0(input$embedding, "_2")) +
        th
    } else {
      feat_col <- input$feature
      is_num   <- is.numeric(df[[feat_col]])
      p <- ggplot2::ggplot(df, ggplot2::aes(Dim1, Dim2, color = .data[[feat_col]])) +
        ggplot2::geom_point(size = input$fp_pt_size, alpha = input$fp_pt_alpha) +
        ggplot2::labs(title = ttl, color = feat_col,
                      x = paste0(input$embedding, "_1"),
                      y = paste0(input$embedding, "_2")) +
        th

      if (is_num) {
        p <- p + ggplot2::scale_color_gradient(low = input$fp_col_low, high = input$fp_col_high)
      } else {
        p <- p + ggplot2::scale_color_brewer(palette = input$fp_palette_cat)
      }
    }

    if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
    p
  })

  output$featureplot <- renderPlot({
    req(sce_obj())
    if (input$feature == "Expression") {
      req(input$genes, length(input$genes) > 0)
    }
    validate(
      need(input$embedding %in% reducedDimNames(sce_obj()),
           paste("La réduction", input$embedding, "n'est pas disponible dans cet objet."))
    )
    featureplot_obj()
  })

  # Downloads FeaturePlot
  handlers_fp <- reactive(make_download_handler(
    featureplot_obj, "FeaturePlot",
    reactive(input$fp_dl_width), reactive(input$fp_dl_height)
  ))
  output$fp_dl_png <- downloadHandler(
    filename = function() paste0("FeaturePlot_", Sys.Date(), ".png"),
    content  = function(file) { png(file, width = input$fp_dl_width, height = input$fp_dl_height, res = 150); print(featureplot_obj()); dev.off() }
  )
  output$fp_dl_pdf <- downloadHandler(
    filename = function() paste0("FeaturePlot_", Sys.Date(), ".pdf"),
    content  = function(file) { pdf(file, width = input$fp_dl_width / 72, height = input$fp_dl_height / 72); print(featureplot_obj()); dev.off() }
  )
  output$fp_dl_svg <- downloadHandler(
    filename = function() paste0("FeaturePlot_", Sys.Date(), ".svg"),
    content  = function(file) { svg(file, width = input$fp_dl_width / 96, height = input$fp_dl_height / 96); print(featureplot_obj()); dev.off() }
  )

  ##############################################################################
  # ---- ViolinPlot ---- #
  ##############################################################################

  output$vln_assay_ui <- renderUI({
    req(sce_obj())
    selectInput("vln_assay", "Assay", choices = assayNames(sce_obj()))
  })

  output$vln_gene_ui <- renderUI({
    req(sce_obj(), input$vln_assay)
    genes <- rownames(assay(sce_obj(), input$vln_assay))
    selectizeInput("vln_genes", "Gène(s)",
                   choices = genes, multiple = TRUE,
                   options = list(maxItems = 8, placeholder = "Tapez un nom de gène..."))
  })

  output$vln_group_ui <- renderUI({
    req(sce_obj())
    cols <- colnames(colData(sce_obj()))
    cat_cols <- cols[sapply(cols, function(c) {
      v <- colData(sce_obj())[[c]]
      is.character(v) || is.factor(v) || (is.numeric(v) && length(unique(v)) <= 30)
    })]
    selectInput("vln_group", "Grouper par", choices = cat_cols)
  })

  violinplot_obj <- reactive({
    req(sce_obj(), input$vln_genes, input$vln_group, length(input$vln_genes) > 0)
    sce <- sce_obj()
    th  <- get_gg_theme(req(input$vln_theme), req(input$vln_base_size))
    ttl <- if (!is.null(input$vln_title) && nchar(input$vln_title) > 0) input$vln_title else NULL

    expr_mat <- assay(sce, input$vln_assay)[input$vln_genes, , drop = FALSE] |>
      as.matrix() |> t() |> as.data.frame()
    df <- cbind(expr_mat, group = as.factor(colData(sce)[[input$vln_group]]))
    df_long <- tidyr::pivot_longer(df, cols = all_of(input$vln_genes),
                                   names_to = "gene", values_to = "Expression")

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = Expression, fill = group)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.85) +
      ggplot2::scale_fill_brewer(palette = input$vln_palette) +
      ggplot2::facet_wrap(~ gene, scales = "free_y") +
      ggplot2::labs(title = ttl, x = input$vln_group, fill = input$vln_group) +
      th +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (isTRUE(input$vln_boxplot)) {
      p <- p + ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.6, fill = "white")
    }

    if (isTRUE(input$vln_show_points)) {
      p <- p + ggplot2::geom_jitter(size = input$vln_pt_size, alpha = input$vln_pt_alpha,
                                     width = 0.15, shape = 16)
    }

    if (isTRUE(input$vln_flip)) p <- p + ggplot2::coord_flip()
    if (!isTRUE(input$vln_legend)) p <- p + ggplot2::theme(legend.position = "none")
    p
  })

  output$violinplot <- renderPlot({
    req(sce_obj(), input$vln_genes, length(input$vln_genes) > 0)
    violinplot_obj()
  })

  output$vln_dl_png <- downloadHandler(
    filename = function() paste0("ViolinPlot_", Sys.Date(), ".png"),
    content  = function(file) { png(file, width = input$vln_dl_width, height = input$vln_dl_height, res = 150); print(violinplot_obj()); dev.off() }
  )
  output$vln_dl_pdf <- downloadHandler(
    filename = function() paste0("ViolinPlot_", Sys.Date(), ".pdf"),
    content  = function(file) { pdf(file, width = input$vln_dl_width / 72, height = input$vln_dl_height / 72); print(violinplot_obj()); dev.off() }
  )
  output$vln_dl_svg <- downloadHandler(
    filename = function() paste0("ViolinPlot_", Sys.Date(), ".svg"),
    content  = function(file) { svg(file, width = input$vln_dl_width / 96, height = input$vln_dl_height / 96); print(violinplot_obj()); dev.off() }
  )

  ##############################################################################
  # ---- DotPlot ---- #
  ##############################################################################

  output$dot_assay_ui <- renderUI({
    req(sce_obj())
    selectInput("dot_assay", "Assay", choices = assayNames(sce_obj()))
  })

  output$dot_genes_ui <- renderUI({
    req(sce_obj(), input$dot_assay)
    genes <- rownames(assay(sce_obj(), input$dot_assay))
    selectizeInput("dot_genes", "Gène(s)",
                   choices = genes, multiple = TRUE,
                   options = list(maxItems = 30, placeholder = "Tapez un nom de gène..."))
  })

  output$dot_group_ui <- renderUI({
    req(sce_obj())
    cols <- colnames(colData(sce_obj()))
    cat_cols <- cols[sapply(cols, function(c) {
      v <- colData(sce_obj())[[c]]
      is.character(v) || is.factor(v) || (is.numeric(v) && length(unique(v)) <= 30)
    })]
    selectInput("dot_group", "Grouper par", choices = cat_cols)
  })

  dotplot_obj <- reactive({
    req(sce_obj(), input$dot_genes, input$dot_group, length(input$dot_genes) > 0)
    sce   <- sce_obj()
    th    <- get_gg_theme(req(input$dot_theme), req(input$dot_base_size))
    ttl   <- if (!is.null(input$dot_title) && nchar(input$dot_title) > 0) input$dot_title else NULL
    genes <- input$dot_genes

    expr_mat <- assay(sce, input$dot_assay)[genes, , drop = FALSE] |> as.matrix()
    group_vec <- as.character(colData(sce)[[input$dot_group]])
    groups    <- unique(group_vec)

    # Calcul expression moyenne et % exprimé
    dot_data <- do.call(rbind, lapply(genes, function(g) {
      do.call(rbind, lapply(groups, function(gr) {
        cells <- which(group_vec == gr)
        vals  <- expr_mat[g, cells]
        data.frame(
          gene        = g,
          group       = gr,
          avg_exp     = mean(vals),
          pct_exp     = 100 * mean(vals > 0),
          stringsAsFactors = FALSE
        )
      }))
    }))

    # Clusterisation optionnelle
    if (isTRUE(input$dot_cluster) && length(genes) > 2) {
      mat_avg <- reshape2::acast(dot_data, gene ~ group, value.var = "avg_exp")
      gene_order <- rownames(mat_avg)[hclust(dist(mat_avg))$order]
      dot_data$gene <- factor(dot_data$gene, levels = gene_order)
    } else {
      dot_data$gene <- factor(dot_data$gene, levels = rev(genes))
    }

    p <- ggplot2::ggplot(dot_data,
           ggplot2::aes(x = group, y = gene,
                        color = avg_exp, size = pct_exp)) +
      ggplot2::geom_point() +
      ggplot2::scale_size(range = c(0, input$dot_scale_size), name = "% exprimé") +
      ggplot2::labs(title = ttl, x = input$dot_group, y = "Gène", color = "Exp. moy.") +
      th +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (input$dot_palette == "viridis") {
      p <- p + ggplot2::scale_color_viridis_c(option = "viridis")
    } else if (input$dot_palette == "plasma") {
      p <- p + ggplot2::scale_color_viridis_c(option = "plasma")
    } else {
      p <- p + ggplot2::scale_color_distiller(palette = input$dot_palette, direction = 1)
    }

    if (isTRUE(input$dot_flip)) p <- p + ggplot2::coord_flip()
    if (!isTRUE(input$dot_legend)) p <- p + ggplot2::theme(legend.position = "none")
    p
  })

  output$dotplot <- renderPlot({
    req(sce_obj(), input$dot_genes, length(input$dot_genes) > 0)
    dotplot_obj()
  })

  output$dot_dl_png <- downloadHandler(
    filename = function() paste0("DotPlot_", Sys.Date(), ".png"),
    content  = function(file) { png(file, width = input$dot_dl_width, height = input$dot_dl_height, res = 150); print(dotplot_obj()); dev.off() }
  )
  output$dot_dl_pdf <- downloadHandler(
    filename = function() paste0("DotPlot_", Sys.Date(), ".pdf"),
    content  = function(file) { pdf(file, width = input$dot_dl_width / 72, height = input$dot_dl_height / 72); print(dotplot_obj()); dev.off() }
  )
  output$dot_dl_svg <- downloadHandler(
    filename = function() paste0("DotPlot_", Sys.Date(), ".svg"),
    content  = function(file) { svg(file, width = input$dot_dl_width / 96, height = input$dot_dl_height / 96); print(dotplot_obj()); dev.off() }
  )

  ##############################################################################
  # ---- Cell Explorer ---- #
  ##############################################################################

  output$ce_embedding_ui <- renderUI({
    req(sce_obj())
    rdnames <- reducedDimNames(sce_obj())
    default <- if ("UMAP" %in% rdnames) "UMAP" else rdnames[1]
    selectInput("ce_embedding", "Réduction", choices = rdnames, selected = default)
  })

  output$ce_color_ui <- renderUI({
    req(sce_obj())
    cols <- colnames(colData(sce_obj()))
    selectInput("ce_color", "Colorer par", choices = cols, selected = cols[1])
  })

  cell_explorer_plot_obj <- reactive({
    req(sce_obj(), input$ce_embedding, input$ce_color)
    sce <- sce_obj()
    th  <- get_gg_theme(req(input$ce_theme), req(input$ce_base_size))
    ttl <- if (!is.null(input$ce_title) && nchar(input$ce_title) > 0) input$ce_title else input$ce_color

    df_emb  <- as.data.frame(reducedDim(sce, input$ce_embedding))[, 1:2]
    colnames(df_emb) <- c("Dim1", "Dim2")
    df_meta <- as.data.frame(colData(sce))
    df      <- cbind(df_emb, df_meta)

    color_col <- input$ce_color
    is_num    <- is.numeric(df[[color_col]])

    p <- ggplot2::ggplot(df, ggplot2::aes(Dim1, Dim2, color = .data[[color_col]])) +
      ggplot2::geom_point(size = input$ce_pt_size, alpha = input$ce_pt_alpha) +
      ggplot2::labs(title = ttl, color = color_col,
                    x = paste0(input$ce_embedding, "_1"),
                    y = paste0(input$ce_embedding, "_2")) +
      th

    if (is_num) {
      p <- p + ggplot2::scale_color_viridis_c()
    } else {
      if (input$ce_palette == "Polychrome") {
        n_col <- length(unique(df[[color_col]]))
        pal   <- scales::hue_pal()(n_col)
        p     <- p + ggplot2::scale_color_manual(values = pal)
      } else {
        p <- p + ggplot2::scale_color_brewer(palette = input$ce_palette)
      }
    }

    # Labels centroids
    if (isTRUE(input$ce_label) && !is_num) {
      centroids <- df |>
        dplyr::group_by(.data[[color_col]]) |>
        dplyr::summarise(Dim1 = median(Dim1), Dim2 = median(Dim2), .groups = "drop")
      p <- p + ggplot2::geom_label(data = centroids,
                                    ggplot2::aes(label = .data[[color_col]]),
                                    size = input$ce_label_size,
                                    color = "black", fill = "white", alpha = 0.8,
                                    label.size = 0.2)
    }

    if (!isTRUE(input$ce_legend)) p <- p + ggplot2::theme(legend.position = "none")
    p
  })

  output$cell_explorer_plot <- renderPlot({
    req(sce_obj(), input$ce_embedding, input$ce_color)
    cell_explorer_plot_obj()
  })

  # Composition cellulaire
  output$cell_composition_plot <- renderPlot({
    req(sce_obj(), input$ce_color)
    sce   <- sce_obj()
    color_col <- input$ce_color
    vals  <- as.factor(colData(sce)[[color_col]])
    df    <- as.data.frame(table(vals))
    colnames(df) <- c("Groupe", "N")
    df    <- df |> dplyr::arrange(dplyr::desc(N)) |>
             dplyr::mutate(Groupe = factor(Groupe, levels = Groupe))

    ggplot2::ggplot(df, ggplot2::aes(x = Groupe, y = N, fill = Groupe)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.3, size = 3) +
      ggplot2::labs(title = "Cellules par groupe", x = NULL, y = "Nombre de cellules") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  })

  output$cell_table <- DT::renderDataTable({
    req(sce_obj(), input$ce_color)
    sce <- sce_obj()
    vals <- as.factor(colData(sce)[[input$ce_color]])
    df   <- as.data.frame(table(Groupe = vals))
    df   <- df |> dplyr::arrange(dplyr::desc(Freq)) |>
            dplyr::mutate(Pourcentage = round(100 * Freq / sum(Freq), 1))
    DT::datatable(df, options = list(pageLength = 10, dom = "tp"), rownames = FALSE)
  })

  output$ce_dl_png <- downloadHandler(
    filename = function() paste0("CellExplorer_", Sys.Date(), ".png"),
    content  = function(file) { png(file, width = input$ce_dl_width, height = input$ce_dl_height, res = 150); print(cell_explorer_plot_obj()); dev.off() }
  )
  output$ce_dl_pdf <- downloadHandler(
    filename = function() paste0("CellExplorer_", Sys.Date(), ".pdf"),
    content  = function(file) { pdf(file, width = input$ce_dl_width / 72, height = input$ce_dl_height / 72); print(cell_explorer_plot_obj()); dev.off() }
  )
  output$ce_dl_svg <- downloadHandler(
    filename = function() paste0("CellExplorer_", Sys.Date(), ".svg"),
    content  = function(file) { svg(file, width = input$ce_dl_width / 96, height = input$ce_dl_height / 96); print(cell_explorer_plot_obj()); dev.off() }
  )

  ##############################################################################
  # ---- QC ---- #
  ##############################################################################

  output$qc_cols_ui <- renderUI({
    req(sce_obj())
    cols <- colnames(colData(sce_obj()))
    num_cols <- cols[sapply(cols, function(c) is.numeric(colData(sce_obj())[[c]]))]
    if (length(num_cols) == 0) return(p("Aucune métrique numérique détectée."))
    checkboxGroupInput("qc_cols", "Métriques à visualiser",
                       choices = num_cols,
                       selected = num_cols[seq_len(min(4, length(num_cols)))])
  })

  output$qc_scatter_x_ui <- renderUI({
    req(sce_obj(), input$qc_scatter)
    cols <- colnames(colData(sce_obj()))
    num_cols <- cols[sapply(cols, function(c) is.numeric(colData(sce_obj())[[c]]))]
    selectInput("qc_scatter_x", "Axe X", choices = num_cols, selected = num_cols[1])
  })

  output$qc_scatter_y_ui <- renderUI({
    req(sce_obj(), input$qc_scatter)
    cols <- colnames(colData(sce_obj()))
    num_cols <- cols[sapply(cols, function(c) is.numeric(colData(sce_obj())[[c]]))]
    selectInput("qc_scatter_y", "Axe Y", choices = num_cols,
                selected = if (length(num_cols) > 1) num_cols[2] else num_cols[1])
  })

  observeEvent(sce_obj(), {
    cols     <- colnames(colData(sce_obj()))
    cat_cols <- cols[sapply(cols, function(c) {
      v <- colData(sce_obj())[[c]]
      is.character(v) || is.factor(v)
    })]
    updateSelectInput(session, "qc_color_by",
                      choices = c("Aucun" = "none", cat_cols))
  })

  qc_violin_obj <- reactive({
    req(sce_obj(), input$qc_cols, length(input$qc_cols) > 0)
    sce <- sce_obj()
    df  <- as.data.frame(colData(sce))[, input$qc_cols, drop = FALSE]
    df_long <- tidyr::pivot_longer(df, cols = everything(), names_to = "Metrique", values_to = "Valeur")

    ggplot2::ggplot(df_long, ggplot2::aes(x = Metrique, y = Valeur, fill = Metrique)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.8) +
      ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.7) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::facet_wrap(~ Metrique, scales = "free") +
      ggplot2::labs(title = "Distribution des métriques QC", x = NULL) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_blank())
  })

  output$qc_violin_plot <- renderPlot({
    req(input$qc_violin, length(input$qc_cols) > 0)
    qc_violin_obj()
  })

  qc_scatter_obj <- reactive({
    req(sce_obj(), input$qc_scatter_x, input$qc_scatter_y)
    sce <- sce_obj()
    df  <- as.data.frame(colData(sce))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[input$qc_scatter_x]],
                                           y = .data[[input$qc_scatter_y]]))

    if (!is.null(input$qc_color_by) && input$qc_color_by != "none") {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[input$qc_color_by]]),
                                    size = input$qc_pt_size, alpha = 0.6) +
               ggplot2::scale_color_brewer(palette = "Set2")
    } else {
      p <- p + ggplot2::geom_point(size = input$qc_pt_size, alpha = 0.5, color = "#3d6b8e")
    }

    p + ggplot2::geom_smooth(method = "lm", color = "#c0392b", se = FALSE, linewidth = 0.7) +
      ggplot2::labs(title = paste(input$qc_scatter_x, "vs", input$qc_scatter_y),
                    x = input$qc_scatter_x, y = input$qc_scatter_y) +
      ggplot2::theme_classic(base_size = 13)
  })

  output$qc_scatter_plot <- renderPlot({
    req(input$qc_scatter, input$qc_scatter_x, input$qc_scatter_y)
    qc_scatter_obj()
  })

  output$qc_dl_png <- downloadHandler(
    filename = function() paste0("QC_", Sys.Date(), ".png"),
    content  = function(file) {
      p <- if (isTRUE(input$qc_violin)) qc_violin_obj() else qc_scatter_obj()
      png(file, width = input$qc_dl_width, height = input$qc_dl_height, res = 150)
      print(p); dev.off()
    }
  )
  output$qc_dl_pdf <- downloadHandler(
    filename = function() paste0("QC_", Sys.Date(), ".pdf"),
    content  = function(file) {
      p <- if (isTRUE(input$qc_violin)) qc_violin_obj() else qc_scatter_obj()
      pdf(file, width = input$qc_dl_width / 72, height = input$qc_dl_height / 72)
      print(p); dev.off()
    }
  )
  output$qc_dl_svg <- downloadHandler(
    filename = function() paste0("QC_", Sys.Date(), ".svg"),
    content  = function(file) {
      p <- if (isTRUE(input$qc_violin)) qc_violin_obj() else qc_scatter_obj()
      svg(file, width = input$qc_dl_width / 96, height = input$qc_dl_height / 96)
      print(p); dev.off()
    }
  )

  ##############################################################################
  # ---- Tableau métadonnées ---- #
  ##############################################################################

  output$table_meta_ui <- renderUI({
    req(sce_obj())
    cols <- colnames(colData(sce_obj()))
    checkboxGroupInput("table_cols", "Colonnes à afficher",
                       choices = cols, selected = cols[seq_len(min(6, length(cols)))])
  })

  output$metadata_table <- DT::renderDataTable({
    req(sce_obj(), input$table_cols)
    df <- as.data.frame(colData(sce_obj()))[, input$table_cols, drop = FALSE]
    DT::datatable(df,
                  filter = "top",
                  extensions = "Buttons",
                  options = list(
                    pageLength = 20,
                    scrollX = TRUE,
                    dom = "Bfrtip",
                    buttons = c("copy", "csv", "excel")
                  ),
                  rownames = TRUE)
  })

  output$table_dl_csv <- downloadHandler(
    filename = function() paste0("metadata_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(input$table_cols)
      df <- as.data.frame(colData(sce_obj()))[, input$table_cols, drop = FALSE]
      write.csv(df, file, row.names = TRUE)
    }
  )

  output$table_dl_xlsx <- downloadHandler(
    filename = function() paste0("metadata_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      req(input$table_cols)
      df <- as.data.frame(colData(sce_obj()))[, input$table_cols, drop = FALSE]
      openxlsx::write.xlsx(df, file, rowNames = TRUE)
    }
  )
}
