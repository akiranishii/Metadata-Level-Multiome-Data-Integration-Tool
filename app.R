# Load additional required package
if (!requireNamespace("purrr", quietly = TRUE)) {
  install.packages("purrr")
  install.packages("DT")
  install.packages("plotly")
}

library(purrr)
library(stringr)
library(DT)
library(tidyverse)
library(plotly)
library(shiny)

# Update UI
ui <- fluidPage(
  titlePanel("RNA-seq, Proteomics, and Lipidomics Data Explorer"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition = "input.tabs === 'Metadata'",
                       selectInput("experiment_type",
                                   "Experiment type:",
                                   c("All", "Bulk RNA-seq", "Proteomics", "Lipidomics"),
                                   selected = "All"),
                       selectInput("species",
                                   "Species:",
                                   c("All", "Mouse", "Human"),
                                   selected = "All"),
                       textInput("search_id", "Search by ID:", ""),
                       actionButton("filter_data", "Filter"),
                       checkboxInput("show_exp_design", "Show experimental design", value = FALSE)
      ),
      conditionalPanel(condition = "input.tabs === 'RNA-seq'",
                       selectInput("rnaseq_dataset", "Select RNA-seq dataset:", NULL),
                       textInput("rnaseq_search_symbol", "Search by Symbol:", ""),
                       numericInput("rnaseq_log2FC_threshold", "min log2FoldChange Threshold:", 0),
                       numericInput("rnaseq_padj_threshold", "max padj Threshold:", 0.05, min = 0, max = 1, step = 0.01)
      ),
      conditionalPanel(condition = "input.tabs === 'Proteomics'",
                       selectInput("proteomics_dataset", "Select Proteomics dataset:", NULL),
                       textInput("proteomics_search_symbol", "Search by Symbol:", ""),
                       numericInput("proteomics_log2FC_threshold", "min log2FoldChange Threshold:", 0),
                       numericInput("proteomics_padj_threshold", "max padj Threshold:", 0.05, min = 0, max = 1, step = 0.01)
      ),
      auth0::logoutButton()
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Metadata", 
                           DTOutput("metadata_table"),
                           downloadButton("download_metadata", "Download Metadata")),
                  tabPanel("RNA-seq",
                           DTOutput("rnaseq_table"),
                           div(style = "height: 20px;"),
                           downloadButton("download_rnaseq", "Download RNA-seq Data"),
                           checkboxInput("show_plot", "MA plot (will not work for 'All datasets')", value = FALSE),
                           uiOutput("gene_selectize_group"),
                           plotlyOutput("rnaseq_plot")),
                  tabPanel("Proteomics",
                           DTOutput("proteomics_table"),
                           div(style = "height: 20px;"),
                           downloadButton("download_proteomics", "Download Proteomics Data"),
                           checkboxInput("show_proteomics_plot", "Volcano plot (will not work for 'All datasets')", value = FALSE),
                           uiOutput("protein_selectize_group"),
                           plotlyOutput("proteomics_plot"))
      )
    )
  )
)

# Update server logic
server <- function(input, output, session) {
  # Read metadata file
  metadata <- reactiveFileReader(1000, NULL, "metadata.csv", read.csv)
  
  # Filter metadata
  filtered_metadata <- eventReactive(input$filter_data, {
    temp_metadata <- metadata()
    
    if (input$experiment_type != "All") {
      temp_metadata <- temp_metadata %>% filter(Experiment.type == input$experiment_type)
    }
    
    if (input$species != "All") {
      temp_metadata <- temp_metadata %>% filter(Species == input$species)
    }
    
    if (input$search_id != "") {
      temp_metadata <- temp_metadata %>% filter(grepl(tolower(input$search_id), tolower(ID), fixed = TRUE))
    }
    
    
    temp_metadata
  }, ignoreNULL = FALSE)
  
  # Render metadata table
  output$metadata_table <- renderDT({
    filtered <- filtered_metadata()
    
    if (!input$show_exp_design) {
      filtered <- filtered %>% select(-Experimental.design)
    }
    
    datatable(filtered,
              options = list(lengthMenu = c(5, 10, 20, 50), pageLength = 10))
  })
  
  # Read RNA-seq datasets
  rna_seq_data <- reactive({
    rna_seq_files <- list.files("rna-seq", pattern = "*.csv", full.names = TRUE)
    rna_seq_list <- lapply(rna_seq_files, function(x) {
      rnaseq_dataset <- read.csv(x)
      rnaseq_dataset$ID <- gsub(".csv", "", basename(x))
      rnaseq_dataset
    })
    rna_seq_list <- setNames(rna_seq_list, gsub(".csv", "", basename(rna_seq_files)))
    rna_seq_list
  })
  
  # Filter RNA-seq data
  filtered_rnaseq_data <- reactive({
    if (input$rnaseq_dataset == "All datasets") {
      rna_seq_combined <- bind_rows(rna_seq_data())
    } else {
      rna_seq_combined <- rna_seq_data()[[input$rnaseq_dataset]]
    }
    if (input$rnaseq_search_symbol != "") {
      rna_seq_combined <- rna_seq_combined %>% filter(grepl(tolower(input$rnaseq_search_symbol), tolower(Symbol), fixed = TRUE))
    }
    
    rna_seq_combined <- rna_seq_combined %>% filter(abs(log2FoldChange) >= input$rnaseq_log2FC_threshold, padj <= input$rnaseq_padj_threshold)
    
    rna_seq_combined
  })
  
  # Update dataset selection input
  observe({
    updateSelectInput(session, "rnaseq_dataset",
                      choices = c("All datasets", gsub(".csv", "", basename(list.files("rna-seq", pattern = "*.csv")))),
                      selected = "All datasets")
  })
  
  # Render RNA-seq table
  output$rnaseq_table <- renderDT({
    datatable(filtered_rnaseq_data(),
              options = list(lengthMenu = c(5, 10, 20, 50), pageLength = 10))
  })
  
  output$gene_selectize_group <- renderUI({
    if (!is.null(filtered_rnaseq_data())) {
      selectizeInput("selected_genes", "Select genes to label on the plot (only filtered genes appear here):",
                     choices = unique(filtered_rnaseq_data()$Symbol),
                     multiple = TRUE,
                     options = list(placeholder = "Search and select genes"))
    } else {
      return(NULL)
    }
  })
  
  
  # return plot data
  plot_data <- reactive({
    if (input$rnaseq_dataset == "All datasets") {
      return(NULL)
    } else {
      data <- rna_seq_data()[[input$rnaseq_dataset]]
      data <- data %>% mutate(log2bm = log2(baseMean + 1))
      plot_data <- data %>% select(Symbol, log2FoldChange, padj, log2bm)
      
      plot_data["group"] <- "NotSignificant(p>0.05)"
      plot_data[which(plot_data['padj'] < 0.05 & plot_data['log2FoldChange'] < 0 ),"group"] <- "Down(p<0.05 & log2FC<0)"
      plot_data[which(plot_data['padj'] < 0.05 & plot_data['log2FoldChange'] > 0 ),"group"] <- "Up(p<0.05 & log2FC>0)"
      
      plot_data["selected"] <- "Not Selected"
      plot_data[which(grepl(tolower(input$rnaseq_search_symbol), tolower(plot_data[['Symbol']]), fixed = TRUE)),"selected"] <- "Selected"
      
      return(plot_data)
    }
  })
  
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  
  x <- list(
    title = "log2 mean expression",
    titlefont = f
  )
  y <- list(
    title = "log2 fold change",
    titlefont = f
  )
  
  output$rnaseq_plot <- renderPlotly({
    if (input$show_plot) {
      if (is.null(plot_data())) {
        return(NULL)
      } else {
        selected_genes <- input$selected_genes
        selected_data <- plot_data() %>% filter(Symbol %in% selected_genes)
        unselected_data <- plot_data() %>% filter(!(Symbol %in% selected_genes))
        
        plot_ly() %>%
          add_trace(data = unselected_data,
                    x = ~log2bm, y = ~log2FoldChange,
                    text = ~paste("Gene name: ", Symbol, '<br>log2FC:', log2FoldChange,
                                  '<br>log2meanexpression:', log2bm, '<br>P:', padj),
                    mode = "markers",
                    color = ~group,
                    marker = list(symbol = "circle")) %>%
          add_trace(data = selected_data,
                    x = ~log2bm, y = ~log2FoldChange,
                    text = ~paste("Gene name: ", Symbol, '<br>log2FC:', log2FoldChange,
                                  '<br>log2meanexpression:', log2bm, '<br>P:', padj),
                    mode = "markers",
                    # color = ~group,
                    marker = list(symbol = "x", color = "red")) %>%
          layout(title = paste(input$rnaseq_dataset, "MA_plot"), xaxis = x, yaxis = y)
      }
    } else {
      return(NULL)
    }
  })
  
  
  # Read proteomics datasets
  proteomics_data <- reactive({
    proteomics_files <- list.files("proteomics", pattern = "*.csv", full.names = TRUE)
    proteomics_list <- lapply(proteomics_files, function(x) {
      proteomics_dataset <- read.csv(x)
      proteomics_dataset$ID <- gsub(".csv", "", basename(x))
      proteomics_dataset
    })
    proteomics_list <- setNames(proteomics_list, gsub(".csv", "", basename(proteomics_files)))
    proteomics_list
  })
  
  # Filter proteomics data
  filtered_proteomics_data <- reactive({
    if (input$proteomics_dataset == "All datasets") {
      proteomics_combined <- bind_rows(proteomics_data())
    } else {
      proteomics_combined <- proteomics_data()[[input$proteomics_dataset]]
    }
    if (input$proteomics_search_symbol != "") {
      proteomics_combined <- proteomics_combined %>% filter(grepl(tolower(input$proteomics_search_symbol), tolower(Symbol), fixed = TRUE))
    }
    proteomics_combined <- proteomics_combined %>% filter(abs(log2FoldChange) >= input$proteomics_log2FC_threshold, padj <= input$proteomics_padj_threshold)
    proteomics_combined
  })
  
  # Update dataset selection input
  observe({
    updateSelectInput(session, "proteomics_dataset",
                      choices = c("All datasets", gsub(".csv", "", basename(list.files("proteomics", pattern = "*.csv")))),
                      selected = "All datasets")
  })
  
  # Render proteomics table
  output$proteomics_table <- renderDT({
    datatable(filtered_proteomics_data(),
              options = list(lengthMenu = c(5, 10, 20, 50), pageLength = 10))
  })
  
  output$protein_selectize_group <- renderUI({
    if (!is.null(filtered_proteomics_data())) {
      selectizeInput("selected_proteins", "Select proteins to label on the plot (only filtered proteins appear here):",
                     choices = unique(filtered_proteomics_data()$Symbol),
                     multiple = TRUE,
                     options = list(placeholder = "Search and select proteins"))
    } else {
      return(NULL)
    }
  })
  
  
  # return plot data
  proteomics_plot_data <- reactive({
    if (input$proteomics_dataset == "All datasets") {
      return(NULL)
    } else {
      data <- proteomics_data()[[input$proteomics_dataset]]
      pseudo_count <- 1e-10
      data <- data %>% mutate(log10_padj = -1*log10(padj + pseudo_count))
      plot_data <- data %>% select(Symbol, log2FoldChange, padj, log10_padj)
      
      plot_data["group"] <- "NotSignificant(p>0.05)"
      plot_data[which(plot_data['padj'] < 0.05 & plot_data['log2FoldChange'] < 0 ),"group"] <- "Down(p<0.05 & log2FC<0)"
      plot_data[which(plot_data['padj'] < 0.05 & plot_data['log2FoldChange'] > 0 ),"group"] <- "Up(p<0.05 & log2FC>0)"
      
      plot_data["selected"] <- "Not Selected"
      plot_data[which(grepl(tolower(input$proteomics_search_symbol), tolower(plot_data[['Symbol']]), fixed = TRUE)),"selected"] <- "Selected"
      
      return(plot_data)
    }
  })
  
  x_p <- list(
    title = "log2 Fold Change",
    titlefont = f
  )
  y_p <- list(
    title = "-log10 padj",
    titlefont = f
  )
  
  output$proteomics_plot <- renderPlotly({
    if (input$show_proteomics_plot) {
      if (is.null(proteomics_plot_data())) {
        return(NULL)
      } else {
        selected_proteins <- input$selected_proteins
        selected_data <- proteomics_plot_data() %>% filter(Symbol %in% selected_proteins)
        unselected_data <- proteomics_plot_data() %>% filter(!(Symbol %in% selected_proteins))
        
        plot_ly() %>%
          add_trace(data = unselected_data,
                    x = ~log2FoldChange, y = ~log10_padj,
                    text = ~paste("Gene name: ", Symbol, '<br>log2FC:', log2FoldChange, '<br>P:', padj),
                    mode = "markers",
                    color = ~group,
                    marker = list(symbol = "circle")) %>%
          add_trace(data = selected_data,
                    x = ~log2FoldChange, y = ~log10_padj,
                    text = ~paste("Gene name: ", Symbol, '<br>log2FC:', log2FoldChange, '<br>P:', padj),
                    mode = "markers",
                    # color = ~group,
                    marker = list(symbol = "x", color = "red")) %>%
          layout(title = paste(input$proteomics_dataset, "MA_plot"), xaxis = x_p, yaxis = y_p)
      }
    } else {
      return(NULL)
    }
  })
  
  
  ####### Add option to download tables
  
  # Download filtered metadata
  output$download_metadata <- downloadHandler(
    filename = function() {
      paste("filtered_metadata", Sys.Date(), ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(filtered_metadata(), file, row.names = FALSE)
    }
  )
  
  # Download filtered RNA-seq data
  output$download_rnaseq <- downloadHandler(
    filename = function() {
      paste("filtered_rnaseq", Sys.Date(), ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(filtered_rnaseq_data(), file, row.names = FALSE)
    }
  )
  
  # Download filtered proteomics data
  output$download_proteomics <- downloadHandler(
    filename = function() {
      paste("filtered_proteomics", Sys.Date(), ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(filtered_proteomics_data(), file, row.names = FALSE)
    }
  )
  
}

# Run the application
auth0::shinyAppAuth0(ui, server)
# shinyApp(ui, server)
                               