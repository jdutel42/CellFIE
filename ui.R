###################################################################################
# ---- CellFIE : Cells & Features Interactive Explorer ---- #
# ui.R — Version 2.0
###################################################################################

library(shiny)
library(ggplot2)
library(SingleCellExperiment)
library(shinycssloaders)
library(dplyr)
library(shinyalert)
library(shinyFeedback)
library(shinyjs)
library(qs2)
library(bslib)

options(shiny.maxRequestSize = 5 * 1024^3)  # 5 Go

# ---- Thème bslib ----
cellfie_theme <- bs_theme(
  version = 5,
  bg = "#f8f9fb",
  fg = "#1a1a2e",
  primary = "#3d6b8e",
  secondary = "#5a9e8f",
  success = "#4caf80",
  warning = "#e09b3d",
  danger = "#c0392b",
  base_font = font_google("Source Sans 3"),
  heading_font = font_google("IBM Plex Sans"),
  code_font = font_google("JetBrains Mono"),
  bootswatch = NULL
)

# ---- Helpers UI ----
download_panel <- function(id_prefix, label = "Figure") {
  tagList(
    hr(),
    tags$b(paste("💾 Télécharger", label)),
    br(), br(),
    fluidRow(
      column(4, downloadButton(paste0(id_prefix, "_png"), "PNG",
                               class = "btn-outline-primary btn-sm w-100")),
      column(4, downloadButton(paste0(id_prefix, "_pdf"), "PDF",
                               class = "btn-outline-danger btn-sm w-100")),
      column(4, downloadButton(paste0(id_prefix, "_svg"), "SVG",
                               class = "btn-outline-secondary btn-sm w-100"))
    ),
    br(),
    fluidRow(
      column(6, numericInput(paste0(id_prefix, "_width"), "Largeur (px)", value = 1200, min = 400, max = 5000, step = 100)),
      column(6, numericInput(paste0(id_prefix, "_height"), "Hauteur (px)", value = 900, min = 300, max = 5000, step = 100))
    )
  )
}

plot_theme_panel <- function(id_prefix) {
  tagList(
    hr(),
    tags$b("🎨 Apparence du plot"),
    br(), br(),
    selectInput(paste0(id_prefix, "_theme"), "Thème ggplot",
                choices = c("Classic" = "classic", "Minimal" = "minimal",
                            "BW" = "bw", "Void" = "void", "Light" = "light"),
                selected = "classic"),
    numericInput(paste0(id_prefix, "_base_size"), "Taille de police", value = 14, min = 8, max = 24, step = 1),
    textInput(paste0(id_prefix, "_title"), "Titre du plot", placeholder = "Laisser vide = auto"),
    checkboxInput(paste0(id_prefix, "_legend"), "Afficher la légende", value = TRUE)
  )
}

###################################################################################
# ---- UI ---- #
###################################################################################

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "https://cdn-icons-png.flaticon.com/512/2932/2932519.png",
             height = "30px", style = "margin-right:8px; vertical-align:middle;"),
    "CellFIE"
  ),
  theme = cellfie_theme,
  bg = "#1a1a2e",
  inverse = TRUE,

  # ---- HEAD extras ----
  header = tagList(
    useShinyjs(),
    shinyFeedback::useShinyFeedback(),
    tags$head(
    tags$style(HTML("
      .navbar-brand { font-family: 'IBM Plex Sans', sans-serif; font-weight:700; font-size:1.3rem; letter-spacing:1px; }
      .card { border:none; box-shadow: 0 2px 12px rgba(0,0,0,0.07); border-radius:10px; }
      .card-header { background: linear-gradient(90deg,#3d6b8e22,#5a9e8f11); font-weight:600; border-radius:10px 10px 0 0 !important; }
      .btn-outline-primary { border-color:#3d6b8e; color:#3d6b8e; }
      .btn-outline-primary:hover { background:#3d6b8e; color:white; }
      .sidebar { background:#f0f4f8 !important; }
      .badge-data { background:#3d6b8e; color:white; font-size:0.75rem; padding:3px 8px; border-radius:20px; }
      .stat-box { background:white; border-radius:8px; padding:12px 16px; box-shadow:0 1px 6px rgba(0,0,0,0.08); text-align:center; }
      .stat-box .stat-val { font-size:1.6rem; font-weight:700; color:#3d6b8e; }
      .stat-box .stat-lab { font-size:0.8rem; color:#666; }
      hr { border-color:#e0e0e0; }
      .shiny-notification { border-radius:8px; }
    "))
  )),

  #######################
  # ---- Home ---- #
  #######################

  nav_panel(
    title = "🏠 Accueil",
    icon = NULL,

    layout_columns(
      col_widths = c(5, 7),

      # ---- Card chargement ----
      card(
        card_header("📂 Charger vos données"),
        card_body(
          p("Formats acceptés : ", tags$code(".rds"), " (Seurat ou SCE) et ", tags$code(".qs")),
          uiOutput("qs_ui"),
          br(),
          shinycssloaders::withSpinner(
            uiOutput("data_summary_ui"),
            type = 4, color = "#3d6b8e"
          )
        )
      ),

      # ---- Card bienvenue + stats ----
      card(
        card_header("👋 Bienvenue dans CellFIE v2.0"),
        card_body(
          p(
            strong("CellFIE"), " est une application interactive pour l'exploration visuelle de données single-cell.",
            " Elle supporte les objets ", tags$code("SingleCellExperiment"), " et ", tags$code("Seurat"), "."
          ),
          br(),
          uiOutput("stats_boxes_ui"),
          br(),
          hr(),
          h6("📋 Workflow recommandé"),
          tags$ol(
            tags$li("Charger votre objet SCE / Seurat (.rds ou .qs)"),
            tags$li("Sélectionner l'assay et la réduction dimensionnelle"),
            tags$li("Explorer les onglets : FeaturePlot, ViolinPlot, DotPlot, Cell Explorer"),
            tags$li("Personnaliser et télécharger vos figures")
          )
        )
      )
    )
  ),

  ##############################
  # ---- FeaturePlot ---- #
  ##############################

  nav_panel(
    title = "🔬 FeaturePlot",

    layout_sidebar(
      sidebar = sidebar(
        title = "Paramètres",
        width = 290,

        h6("📌 Données"),
        uiOutput("assay_ui"),
        uiOutput("embedding_ui"),
        uiOutput("feature_ui"),
        shinycssloaders::withSpinner(uiOutput("gene_ui"), type = 4, color = "#3d6b8e"),

        hr(),
        h6("⚙️ Points"),
        sliderInput("fp_pt_size", "Taille des points", min = 0.05, max = 3, value = 0.6, step = 0.05),
        sliderInput("fp_pt_alpha", "Transparence", min = 0.1, max = 1, value = 0.85, step = 0.05),
        checkboxInput("fp_raster", "Rasteriser (> 50k cellules)", value = FALSE),

        hr(),
        h6("🎨 Couleurs (Expression)"),
        fluidRow(
          column(6, colourpicker::colourInput("fp_col_low", "Min", value = "#d3d3d3")),
          column(6, colourpicker::colourInput("fp_col_high", "Max", value = "#08306b"))
        ),
        selectInput("fp_palette_cat", "Palette (métadonnées catégoriques)",
                    choices = c("Set2", "Set1", "Dark2", "Paired", "Tableau10" = "Tableau 10",
                                "Polychrome36" = "alphabet"),
                    selected = "Set2"),

        plot_theme_panel("fp"),
        download_panel("fp_dl", "FeaturePlot")
      ),

      # ---- Main ----
      card(
        card_header(uiOutput("fp_card_header")),
        card_body(
          shinycssloaders::withSpinner(
            plotOutput("featureplot", height = "680px"),
            type = 6, color = "#3d6b8e"
          )
        )
      )
    )
  ),

  ##############################
  # ---- ViolinPlot ---- #
  ##############################

  nav_panel(
    title = "🎻 ViolinPlot",

    layout_sidebar(
      sidebar = sidebar(
        title = "Paramètres",
        width = 290,

        h6("📌 Données"),
        uiOutput("vln_assay_ui"),
        uiOutput("vln_gene_ui"),
        uiOutput("vln_group_ui"),

        hr(),
        h6("⚙️ Options"),
        checkboxInput("vln_show_points", "Afficher les points", value = TRUE),
        sliderInput("vln_pt_size", "Taille des points", 0.1, 2, 0.3, 0.1),
        sliderInput("vln_pt_alpha", "Transparence points", 0.05, 1, 0.3, 0.05),
        checkboxInput("vln_boxplot", "Ajouter boxplot", value = TRUE),
        checkboxInput("vln_flip", "Inverser les axes", value = FALSE),
        selectInput("vln_palette", "Palette",
                    choices = c("Set2", "Set1", "Dark2", "Paired"),
                    selected = "Set2"),

        plot_theme_panel("vln"),
        download_panel("vln_dl", "ViolinPlot")
      ),

      card(
        card_header("Expression par groupe cellulaire"),
        card_body(
          shinycssloaders::withSpinner(
            plotOutput("violinplot", height = "650px"),
            type = 6, color = "#3d6b8e"
          )
        )
      )
    )
  ),

  ##############################
  # ---- DotPlot ---- #
  ##############################

  nav_panel(
    title = "🔵 DotPlot",

    layout_sidebar(
      sidebar = sidebar(
        title = "Paramètres",
        width = 290,

        h6("📌 Données"),
        uiOutput("dot_assay_ui"),
        uiOutput("dot_genes_ui"),
        uiOutput("dot_group_ui"),

        hr(),
        h6("⚙️ Options"),
        sliderInput("dot_scale_size", "Taille max des points", 2, 15, 8, 0.5),
        checkboxInput("dot_flip", "Inverser les axes", value = TRUE),
        checkboxInput("dot_cluster", "Clusteriser gènes", value = FALSE),
        selectInput("dot_palette", "Palette couleur",
                    choices = c("Blues", "RdYlBu", "Viridis" = "viridis",
                                "Plasma" = "plasma", "YlOrRd"),
                    selected = "Blues"),

        plot_theme_panel("dot"),
        download_panel("dot_dl", "DotPlot")
      ),

      card(
        card_header("Expression moyenne et proportion d'expression par groupe"),
        card_body(
          shinycssloaders::withSpinner(
            plotOutput("dotplot", height = "650px"),
            type = 6, color = "#3d6b8e"
          )
        )
      )
    )
  ),

  ##############################
  # ---- Cell Explorer ---- #
  ##############################

  nav_panel(
    title = "🧬 Cell Explorer",

    layout_sidebar(
      sidebar = sidebar(
        title = "Paramètres",
        width = 290,

        h6("📌 Affichage"),
        uiOutput("ce_embedding_ui"),
        uiOutput("ce_color_ui"),

        hr(),
        h6("⚙️ Options"),
        sliderInput("ce_pt_size", "Taille des points", 0.05, 3, 0.6, 0.05),
        sliderInput("ce_pt_alpha", "Transparence", 0.1, 1, 0.85, 0.05),
        checkboxInput("ce_label", "Afficher les labels", value = TRUE),
        sliderInput("ce_label_size", "Taille des labels", 2, 8, 4, 0.5),
        selectInput("ce_palette", "Palette",
                    choices = c("Set2", "Set1", "Dark2", "Paired", "Polychrome"),
                    selected = "Set2"),
        checkboxInput("ce_legend", "Afficher la légende", value = TRUE),

        plot_theme_panel("ce"),
        download_panel("ce_dl", "CellExplorer")
      ),

      layout_columns(
        col_widths = c(8, 4),

        card(
          card_header("Embedding — Métadonnées"),
          card_body(
            shinycssloaders::withSpinner(
              plotOutput("cell_explorer_plot", height = "620px"),
              type = 6, color = "#3d6b8e"
            )
          )
        ),

        card(
          card_header("📊 Composition cellulaire"),
          card_body(
            shinycssloaders::withSpinner(
              plotOutput("cell_composition_plot", height = "300px"),
              type = 6, color = "#3d6b8e"
            ),
            hr(),
            shinycssloaders::withSpinner(
              DT::dataTableOutput("cell_table"),
              type = 4, color = "#3d6b8e"
            )
          )
        )
      )
    )
  ),

  ##############################
  # ---- QC ---- #
  ##############################

  nav_panel(
    title = "🧪 QC",

    layout_sidebar(
      sidebar = sidebar(
        title = "Métriques QC",
        width = 280,

        uiOutput("qc_cols_ui"),
        hr(),
        h6("⚙️ Options"),
        checkboxInput("qc_violin", "ViolinPlot", value = TRUE),
        checkboxInput("qc_scatter", "ScatterPlot", value = TRUE),
        uiOutput("qc_scatter_x_ui"),
        uiOutput("qc_scatter_y_ui"),
        selectInput("qc_color_by", "Colorer par :", choices = c("Aucun" = "none")),
        sliderInput("qc_pt_size", "Taille des points", 0.1, 2, 0.5, 0.1),

        download_panel("qc_dl", "QC")
      ),

      layout_columns(
        col_widths = c(6, 6),

        card(
          card_header("Distribution des métriques QC"),
          shinycssloaders::withSpinner(
            plotOutput("qc_violin_plot", height = "400px"),
            type = 6, color = "#3d6b8e"
          )
        ),

        card(
          card_header("Corrélation entre métriques"),
          shinycssloaders::withSpinner(
            plotOutput("qc_scatter_plot", height = "400px"),
            type = 6, color = "#3d6b8e"
          )
        )
      )
    )
  ),

  ##############################
  # ---- Statistics Table ---- #
  ##############################

  nav_panel(
    title = "📋 Tableau",

    layout_sidebar(
      sidebar = sidebar(
        title = "Filtres",
        width = 260,

        uiOutput("table_meta_ui"),
        uiOutput("table_filter_ui"),
        hr(),
        downloadButton("table_dl_csv", "📥 Télécharger CSV", class = "btn-outline-primary w-100"),
        br(), br(),
        downloadButton("table_dl_xlsx", "📥 Télécharger Excel", class = "btn-outline-success w-100")
      ),

      card(
        card_header("Métadonnées cellulaires"),
        card_body(
          shinycssloaders::withSpinner(
            DT::dataTableOutput("metadata_table"),
            type = 4, color = "#3d6b8e"
          )
        )
      )
    )
  ),

  ##############################
  # ---- About ---- #
  ##############################

  nav_panel(
    title = "ℹ️ À propos",

    layout_columns(
      col_widths = c(6, 6),

      card(
        card_header("CellFIE — Cells & Features Interactive Explorer"),
        card_body(
          h5("Version 2.0"),
          p("Application Shiny pour l'exploration interactive de données single-cell."),
          hr(),
          h6("📦 Fonctionnalités"),
          tags$ul(
            tags$li("Chargement SCE / Seurat (.rds, .qs)"),
            tags$li("FeaturePlot : expression génique sur embedding"),
            tags$li("ViolinPlot : distribution par groupe"),
            tags$li("DotPlot : expression moyenne × proportion"),
            tags$li("Cell Explorer : métadonnées sur embedding"),
            tags$li("QC : métriques qualité"),
            tags$li("Export PNG / PDF / SVG haute résolution")
          ),
          hr(),
          h6("🔧 Dépendances R"),
          tags$code("shiny, bslib, ggplot2, SingleCellExperiment, qs2, DT, colourpicker, scales, tidyr, dplyr")
        )
      ),

      card(
        card_header("📘 Guide d'utilisation"),
        card_body(
          tags$ol(
            tags$li(strong("Charger"), " votre objet SCE ou Seurat via l'onglet Accueil"),
            tags$li(strong("FeaturePlot"), " : visualiser l'expression d'un ou plusieurs gènes sur l'embedding"),
            tags$li(strong("ViolinPlot"), " : comparer la distribution d'expression entre populations"),
            tags$li(strong("DotPlot"), " : vue d'ensemble de plusieurs gènes × groupes"),
            tags$li(strong("Cell Explorer"), " : explorer les annotations cellulaires"),
            tags$li(strong("QC"), " : vérifier les métriques qualité (nCounts, nFeatures, pct_mito...)"),
            tags$li(strong("Tableau"), " : explorer et exporter les métadonnées")
          ),
          hr(),
          p("💡 Tous les plots sont personnalisables (couleurs, thème, taille) et téléchargeables en PNG, PDF ou SVG.")
        )
      )
    )
  ),

  footer = tags$div(
    style = "position:fixed; bottom:10px; right:16px; z-index:9999;",
    uiOutput("data_status_badge")
  )
)
