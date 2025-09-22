# app.R

library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(tidyverse)

build_enrichment_network <- function(data, min_overlap = 0.2) {
  # Parse genes per term
  data$gene_list <- strsplit(data$geneID, "/")
  term_ids <- data$Description
  
  edges <- expand.grid(from = seq_along(term_ids), to = seq_along(term_ids)) %>%
    filter(from < to) %>%
    mutate(
      overlap = mapply(function(i, j) {
        a <- data$gene_list[[i]]
        b <- data$gene_list[[j]]
        length(intersect(a, b)) / length(union(a, b))
      }, from, to)
    ) %>%
    filter(overlap >= min_overlap)
  
  # Create igraph object
  g <- graph_from_data_frame(
    d = data.frame(
      from = term_ids[edges$from],
      to = term_ids[edges$to],
      weight = edges$overlap
    ),
    vertices = data.frame(
      name = term_ids,
      p.adjust = data$p.adjust,
      Count = data$Count
    ),
    directed = FALSE
  )
  
  return(g)
}

plot_enrichment_network <- function(g, title = "Enrichment Network") {
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = weight), color = "grey70") +
    geom_node_point(aes(size = Count, color = p.adjust)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_viridis_c() +
    theme_void() +
    ggtitle(title)
}



options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB limit

ui <- fluidPage(
  navbarPage(
    "ORA Results Visualization App",
    
    # --- Upload Tab ---
    tabPanel("Input data",
             sidebarLayout(
               sidebarPanel(
                 fileInput("input_data_file", "Upload .rds File", accept = ".rds")
               ),
               mainPanel(
                 verbatimTextOutput("file_info"),
                 tableOutput("data_preview")
               )
             )
    ),
    
    # --- Barplot Tab ---
    tabPanel("Barplot", 
             sidebarLayout(
               sidebarPanel(
                 textInput("plot_title", "Plot title:", "GO Pathway Enrichment"),
                 sliderInput("title_size", "Title size:", 1, 30, 20, 1),
                 selectInput("Gene_Ontology", "Gene Ontology:",
                             choices = c("ALL", "CC", "BP", "MF"), selected = "ALL"),
                 selectInput("Color", "Choose colors",
                             choices = c("Standard", "Monochromatic", "Custom"), selected = "Standard"),
                 conditionalPanel(
                   condition = "input.Color == 'Custom'",
                   textInput("color1", "First color", "#993399"),
                   textInput("color2", "Second color", "#FFFF00")
                 ),
                 selectInput("x_content", "X axis:", c("GeneRatio", "Count"), selected = "GeneRatio"),
                 sliderInput("n_terms", "Number of terms:", 1, 50, 10),
                 sliderInput("font_size", "Font size:", 1, 30, 10),
                 selectInput("extension", "File extension:", c("pdf", "png")),
                 downloadButton("download_barplot", "Download Barplot")
               ),
               mainPanel(plotOutput("bar_plot"))
             )
    ),
    
    # --- Dotplot Tab ---
    tabPanel("Dotplot", 
             sidebarLayout(
               sidebarPanel(
                 textInput("plot_title", "Plot title:", "GO Pathway Enrichment"),
                 sliderInput("title_size", "Title size:", 1, 30, 20),
                 selectInput("Gene_Ontology", "Gene Ontology:",
                             choices = c("ALL", "CC", "BP", "MF"), selected = "ALL"),
                 selectInput("Color", "Choose colors",
                             choices = c("Standard", "Monochromatic", "Custom"), selected = "Standard"),
                 conditionalPanel(
                   condition = "input.Color == 'Custom'",
                   textInput("color1", "First color", "#993399"),
                   textInput("color2", "Second color", "#FFFF00")
                 ),
                 sliderInput("n_terms", "Number of terms:", 1, 50, 10),
                 sliderInput("font_size", "Font size:", 1, 30, 10),
                 selectInput("extension", "File extension:", c("pdf", "png")),
                 downloadButton("download_dotplot", "Download Dotplot")
               ),
               mainPanel(plotOutput("dot_plot"))
             )
    ),
    
    tabPanel("Network Plot",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("min_edge", "Min % gene overlap", min = 0, max = 1, value = 0.2, step = 0.05),
                 sliderInput("n_terms", "Number of terms", min = 2, max = 50, value = 10),
                 textInput("plot_title", "Title", value = "Enrichment Network")
               ),
               mainPanel(
                 plotOutput("net_plot")
               )
             )
    )
    
  )
)

server <- function(input, output, session) {
  
  input_data <- reactive({
    req(input$input_data_file)
    readRDS(input$input_data_file$datapath)
  })
  
  output$file_info <- renderPrint({
    req(input$input_data_file)
    input$input_data_file
  })
  
  output$data_preview <- renderTable({
    head(input_data())
  })
  
  filtered_data <- reactive({
    data <- input_data()
    
    # Convert GeneRatio from "x/y" to numeric
    if ("GeneRatio" %in% names(data) && is.character(data$GeneRatio)) {
      data$GeneRatio <- sapply(strsplit(data$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    }
    
    if (input$Gene_Ontology == "ALL") return(data)
    data %>% filter(ONTOLOGY == input$Gene_Ontology)
  })
  
  
  color_scale <- reactive({
    if (input$Color == "Monochromatic") {
      scale_fill_gradient(low = "#08306B", high = "#CCFFFF")
    } else if (input$Color == "Custom") {
      scale_fill_gradient(low = input$color1, high = input$color2)
    } else {
      scale_fill_continuous()  # default
    }
  })
  
  output$bar_plot <- renderPlot({
    df <- filtered_data() %>%
      arrange(p.adjust) %>%
      head(input$n_terms)
    
    ggplot(df, aes(
      x = .data[[input$x_content]],
      y = reorder(Description, .data[[input$x_content]]),
      fill = p.adjust
    )) +
      geom_bar(stat = "identity") +
      labs(
        title = input$plot_title,
        x = input$x_content,
        y = "Term"
      ) +
      color_scale() +
      theme_minimal(base_size = input$font_size) +
      theme(plot.title = element_text(size = input$title_size, face = "bold"))
  })
  
  output$dot_plot <- renderPlot({
    df <- filtered_data() %>%
      arrange(desc(Count)) %>%
      head(input$n_terms)
    
    ggplot(df, aes(
      x = GeneRatio,
      y = reorder(Description, GeneRatio),
      size = Count,
      color = p.adjust
    )) +
      geom_point() +
      labs(
        title = input$plot_title,
        x = "GeneRatio",
        y = "Term"
      ) +
      color_scale() +
      theme_minimal(base_size = input$font_size) +
      theme(plot.title = element_text(size = input$title_size, face = "bold"))
  })
  
  output$download_barplot <- downloadHandler(
    filename = function() {
      paste0("barplot_ORA_", Sys.Date(), ".", input$extension)
    },
    content = function(file) {
      device <- if (input$extension == "pdf") pdf else png
      device(file)
      print(output$bar_plot())
      dev.off()
    }
  )
  
  output$download_dotplot <- downloadHandler(
    filename = function() {
      paste0("dotplot_ORA_", Sys.Date(), ".", input$extension)
    },
    content = function(file) {
      device <- if (input$extension == "pdf") pdf else png
      device(file)
      print(output$dot_plot())
      dev.off()
    }
  )
  
  netplot_data <- reactive({
    df <- filtered_data() %>%
      arrange(p.adjust) %>%
      head(input$n_terms)
    
    build_enrichment_network(df, min_overlap = input$min_edge)
  })
  
  output$net_plot <- renderPlot({
    g <- netplot_data()
    plot_enrichment_network(g, title = input$plot_title)
  })
  
}

shinyApp(ui = ui, server = server)

