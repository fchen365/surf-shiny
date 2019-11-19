# Load packages
library(shiny)
library(shinyBS)
library(shinyLP)
library(shinythemes)
options(repos = BiocManager::repositories())

source("src/utilities.R")
source("src/parse.R") 
source("src/drseq.R") 
source("src/associate.R")
source("src/daseq.R")

## 104 RBP names
targets <- read.delim("data/RBP.txt", stringsAsFactors = F)$RBP 

# Load surf results
# surf.results <- list(AQR = readRDS(paste0("data/AQR.results.rds")))
surf.results <- lapply(targets, function(target)
  readRDS(paste0("data/",target,".results.rds")))
names(surf.results) = targets

## gene_name
gene_name = read.table("data/gene_name.txt", header = F, 
                       col.names = c("gene_id","gene_name"), 
                       row.names = 1, stringsAsFactors = F)

# Define UI
ui <- navbarPage(
  "SURF", 
  theme = shinytheme("flatly"),
  tabPanel(
    "Home", icon = icon("home"),
    jumbotron(
      "Welcome to SURF!", 
      "The Statistical Utility for RBP Functions (SURF) is an integrative analysis framework to identify alternative splicing (AS), alternative transcription initiation (ATI), and alternative polyadenylation (APA) events regulated by individual RBPs and elucidate RNA-RBP interactions governing these events. This shiny app presents SURF results on 104 RBPs available from the ENCODE consortium.",
      button = FALSE),
    fluidRow(
      column(
        width = 6, 
        panel_div(
          class_type = "primary", 
          panel_title = "Background",
          content = "Post-transcriptional regulation by RNA binding proteins (RBPs) is a major contributor to protein diversity in mammalian genomes. RBPs interact with pre-mature messenger RNA transcripts and orchestrate formation of mature RNA transcripts through regulation of alternative splicing events such as exon skipping, 3’ or 5’ splicing, intron retention, alternative transcription initiation, and alternative polyadenylation. Recent advances in ultraviolet cross-linking immunoprecipitation followed by high throughput sequencing (CLIP-seq) resulted in large collections of RBP binding data coupled with transcriptome profiling by RNA-seq across multiple conditions (e.g. ENCODE).")),
      column(
        width = 5, 
        panel_div(
          class_type = "primary", 
          "How does it work?",
          content = "SURF leverages large-scale CLIP-seq and RNA-seq data with and without RNA interference screening and infers rules of RBPs in alternative transcriptional regulation (ATR), including AS, ATI, and APA. For an overview figure of the SURF framework, please go on the bottom of this Home page. The multi-moduled SURF first extents the versatile differential exon usage analysis method DEXSeq for detection of differential ATR events and associates these events to local RNA-RBP interactions as measured by CLIP-seq."
        )
      )
    ),
    fluidRow(
      column(
        width = 10, 
        panel_div(
          class_type = "primary", 
          panel_title = "What does it do?",
          content = fluidRow(
            column(
              width = 4, 
              panel_div(
                class_type = "success", 
                panel_title = "ATR Event",
                content = "This app allows you to search ATR events relevant to a specific RBP. You can search by gene names, ATR event types, etc."
              )
            ),
            column(
              width = 4, 
              panel_div(
                class_type = "success", 
                panel_title = "Volcano Plot",
                content = "Volcano plot provides you the global view of what ATR event types an RBP is likely to regulate."
              )
            ),
            column(
              width = 4, 
              panel_div(
                class_type = "success", 
                "FA Plot",
                content = "Functional association plot depicts the potential genomic location where RBP-RNA interaction takes place in ATR events."
              )
            )
          )
        )
      )
      
    ),  # end of fluidRow
    fluidRow(
      column(
        width = 5, 
        panel_div(
          class_type = "primary", 
          "Citation",
          HTML("Fan Chen and Sunduz Keles. “SURF: Integrative analysis of a compendium of RNA-seq and eCLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins.” ")
        )
      ),
      column(
        width = 3, 
        panel_div(
          class_type = "primary", 
          "Contact and License",
          HTML("Email us: <a href='mailto:fan.chen@wisc.edu'>Fan Chen</a>, <a href='mailto:keles@stat.wisc.edu'>Sunduz Keles</a><br><br>Copyright (c) 2019 Fan Chen, Sunduz Keles")
        )
      )
    ),
    fluidRow(
      column(
        width = 10, 
        textOutput(outputId = "pipeline.desc"),
        imageOutput("pipeline")
      )
    )
    # fluidPage(
    #   titlePanel("Statistical Utility for RBP Functions (SURF)"),
    #   sidebarLayout(
    #     # Sidebar with a slider input
    #     sidebarPanel(
    #       textOutput(outputId = "citation")
    #     ),
    # 
    #     # Show a plot of the generated distribution
    #     mainPanel(
    #       textOutput(outputId = "abstract"),
    #       textOutput(outputId = "pipeline.desc"),
    #       imageOutput("pipeline")
    #     )
    #   )
    # )
  ),
  tabPanel(
    "ATR Event", icon = icon("search"),
    fluidPage(
      titlePanel("Differential ATR Event"),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "atr.factor",
            label = strong("RNA-binding protein:"),
            choices = sort(targets),
            selected = "AQR"
          ),
          selectInput(
            inputId = "atr.event",
            label = strong("ATR event:"),
            choices = surf.events,
            selected = c("SE", "AFE"),
            multiple = T
          ),
          selectInput(
            inputId = "atr.group",
            label = strong("Differential group:"),
            choices = c("increase", "decrease", "no change"),
            selected = c("increase", "decrease"),
            multiple = T
          ),
          sliderInput("atr.fdr.cutoff", "Significance level:", 
                      min = 0, max = 0.05, value = 0.05)
        ),
        mainPanel(
          DT::dataTableOutput("surfData")
        )
      )
    )
  ),
  tabPanel(
    "Volcano Plot", icon = icon("chart-area"),
    fluidPage(
      titlePanel("Volcano Plot"),
      sidebarLayout(
        sidebarPanel(
          # # Input side bar
          selectInput(
            inputId = "volcano.factor",
            label = strong("RNA-binding protein:"),
            choices = sort(targets),
            selected = "AQR"
          ),
          sliderInput("volcano.fdr.cutoff", "DrSeq FDR cut-off:",
                      min = 0, max = 0.4, value = 0.01),
          
          sliderInput("volcano.lfc.cutoff", "|log2 fold change| cut-off:",
                      min = 0, max = 2.5, value = 1.0, round = F),
        ),
        # Output: Description, lineplot, and reference
        mainPanel(
          plotOutput(outputId = "volcano.plot", height = "400px"),
          textOutput(outputId = "volcano.desc")
        )
      )
    )
  ),
  tabPanel(
    "FA Plot", icon = icon("chart-line"),
    fluidPage(
      titlePanel("Functional Association (FA) Plot"),
      sidebarLayout(
        # Input side bar
        sidebarPanel(
          selectInput(
            inputId = "fa.factor",
            label = strong("RNA-binding protein:"),
            choices = sort(targets),
            selected = "AQR"
          ),
          
          selectInput(
            inputId = "fa.plot.event",
            label = strong("ATR event:"),
            choices = surf.events,
            selected = c("SE", "AFE"),
            multiple = T
          ),
          
          sliderInput("fa.trim", "Trim (quantile):",  
                      min = 0, max = 0.5, value = 0.025),
          
          sliderInput("fa.fdr.cutoff", "Significance level:", 
                      min = 0, max = 0.2, value = 0.05)
        ),
        # Output: Description, lineplot, and reference
        mainPanel(
          plotOutput(outputId = "fa.plot", height = "400px"), 
          textOutput(outputId = "fa.desc"),
          textOutput(outputId = "feature.desc"),
          imageOutput("feature_SE", height = "auto"),
          imageOutput("feature_RI", height = "auto"),
          imageOutput("feature_A3SS", height = "auto"),
          imageOutput("feature_A5SS", height = "auto"),
          imageOutput("feature_AFE", height = "auto"),
          imageOutput("feature_A5U", height = "auto"),
          imageOutput("feature_IAP", height = "auto"),
          imageOutput("feature_TAP", height = "auto"),
        )
      )
    )
  )
)

# Define server function
server <- function(input, output, session) {
  
  output$citation <- renderText("Fan Chen and Sunduz Keles. \n\"SURF: Integrative analysis of a compendium of RNA-seq and eCLIP-seq datasets highlights complex governing of alternative transcriptional regulation by RNA-binding proteins.\"")
  output$abstract <- renderText("Post-transcriptional regulation by RNA binding proteins (RBPs) is a major contributor to protein diversity in mammalian genomes. RBPs interact with pre-mature messenger RNA transcripts and orchestrate formation of mature RNA transcripts through regulation of alternative splicing events such as exon skipping, 3’ or 5’ splicing, intron retention, alternative transcription initiation, and alternative polyadenylation. Recent advances in ultraviolet cross-linking immunoprecipitation followed by high throughput sequencing (CLIP-seq) resulted in large collections of RBP binding data coupled with transcriptome profiling by RNA-seq across multiple conditions. We leveraged such a large collection of CLIP-seq and RNA-seq data with and without RNA interference screening from the ENCODE consortium and developed SURF, Statistical Utility for RBP Functions. SURF is an integrative analysis framework to identify alternative splicing (AS), alternative transcription initiation (ATI), and alternative polyadenylation (APA) events regulated by individual RBPs and elucidate RNA-RBP interactions governing these events. The multi-moduled SURF first extents the versatile differential exon usage analysis method DEXSeq for detection of differential alternative transcriptional regulation (ATR) of AS, ATI, and APA events and associates these events to local RNA-RBP interactions as measured by CLIP-seq. Large-scale application of SURF recovered known roles of a handful of RBPs while generating novel hypotheses for the others under a well-calibrated false discovery rate. Downstream analysis of RBP-RNA interaction regions in SURF-identified associations exhibited significant enrichment of somatic mutations in both the TCGA and ICGC datasets. Furthermore, a systematic comparison of SURF-identified transcript targets of RBPs across GTEx and TCGA compendia highlighted specific AS, ATI, and APA regulation roles for RBPs in adult acute myeloid leukemia.")
  output$pipeline.desc <- renderText("Here is an overview of SURF pipeline.")
  output$pipeline <- renderImage({
    width <- session$clientData$output_pipeline_width * .89
    list(src = "figure/1_pipeline.png",
         width = min(800, width),
         alt = "An overview of SURF pipeline.")
  }, deleteFile = F)
  
  
  ## diff ATR event 
  output$surfData <- DT::renderDataTable({
    sdat <- surf.results[[input$atr.factor]]@trainData
    sdat <- sdat[sdat$event_name %in% input$atr.event &
                   sdat$group %in% input$atr.group &
                   sdat$padj < input$atr.fdr.cutoff, ]
    gr <- range(sdat$genomicData) %>% unlist %>% data.frame
    data.frame(gr, sdat[,-c(5,6,17)]) %>% 
      mutate("adjusted p-value" = round(padj, digits = 7),
             "log2 fold change" = round(logFoldChange, digits = 5), 
             gene = gene_name[gene_id, "gene_name"], 
             event = event_name, 
             difference = group, 
             ) %>%
      dplyr::select(c(1:3,5,23,22,20,21,24)) %>% 
      DT::datatable(options = list(lengthMenu = c(50, 30, 50), pageLength = 50))
  })
  
  ## volcano
  output$volcano.plot <- renderPlot({
    volcano.plot(surf.results[[input$volcano.factor]]@trainData, 
                 lfc.cutoff = c(-1,1) * input$volcano.lfc.cutoff, 
                 fdr.cutoff = input$volcano.fdr.cutoff)
  })
  output$volcano.desc <- renderText({
    paste0("Volcano plot (-log10 transformed adjusted p-value versus log2 of fold change) of DrSeq results for ", input$volcano.factor, " shRNA-seq dataset (ENCODE), stratified by ATR event types.")
  })
  
  # FA plot
  output$feature.desc <- renderText({
    paste0("The location features for ", 
           paste(input$fa.plot.event, collapse = ", "), 
           " are illustrated below.")
  })
  output$feature_SE <- renderImage({
    if ("SE" %in% input$fa.plot.event)
    list(src = ifelse("SE" %in% input$fa.plot.event, 
                      "figure/feature_SE.png", ""),
         width = session$clientData$output_feature_SE_width * .6)
  }, deleteFile = F)
  output$feature_RI <- renderImage({
    list(src = ifelse("RI" %in% input$fa.plot.event, 
                      "figure/feature_RI.png", ""),
         width = session$clientData$output_feature_RI_width * .6)
  }, deleteFile = F)
  output$feature_A3SS <- renderImage({
    list(src = ifelse("A3SS" %in% input$fa.plot.event, 
                      "figure/feature_A3SS.png", ""),
         width = session$clientData$output_feature_A3SS_width * .6)
  }, deleteFile = F)
  output$feature_A5SS <- renderImage({
    list(src = ifelse("A5SS" %in% input$fa.plot.event, 
                      "figure/feature_A5SS.png", ""),
         width = session$clientData$output_feature_A5SS_width * .6)
  }, deleteFile = F)
  output$feature_AFE <- renderImage({
    list(src = ifelse("AFE" %in% input$fa.plot.event, 
                      "figure/feature_AFE.png", ""),
         width = session$clientData$output_feature_AFE_width * .6)
  }, deleteFile = F)
  output$feature_A5U <- renderImage({
    list(src = ifelse("A5U" %in% input$fa.plot.event, 
                      "figure/feature_A5U.png", ""),
         width = session$clientData$output_feature_A5U_width * .6)
  }, deleteFile = F)
  output$feature_IAP <- renderImage({
    list(src = ifelse("IAP" %in% input$fa.plot.event, 
                      "figure/feature_IAP.png", ""),
         width = session$clientData$output_feature_IAP_width * .6)
  }, deleteFile = F)
  output$feature_TAP <- renderImage({
    list(src = ifelse("TAP" %in% input$fa.plot.event, 
                      "figure/feature_TAP.png", ""),
         width = session$clientData$output_feature_TAP_width * .6)
  }, deleteFile = F)
  output$fa.plot <- renderPlot({
    fa.plot(surf.results[[input$fa.factor]], 
            plot.event = input$fa.plot.event, 
            trim = input$fa.trim, 
            fdr.cutoff = input$fa.fdr.cutoff)
  })
  output$fa.desc <- renderText({
    paste0("Functional association (FA) plot for ", input$fa.factor, " in ", 
           paste(input$fa.plot.event, collapse = ", "), ". ", 
           "The upper panel plots the feature signals for three differential ATR groups.", 
           "The lower panel depicts the -log10 transformed p-values for each association testing after multiplicity correction.",
           "The dashed lines indicate the significance level (of FDR).")
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)