rm(list = ls())
library(shiny)
library(shinydashboard)
library(readxl)
library(dplyr)
library(shinyjs)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(DT)
library(tidyverse)
library(enrichplot)
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("ReactomePA")
library(ReactomePA)
# devtools::install_github('davidgohel/gdtools')
# devtools::install_github("r-lib/svglite")
library(svglite)
library(limma)
library(preprocessCore)
library(gplots)
library(ggplot2)

library(RColorBrewer)
library(DT)

library(dplyr)
library(data.table)
library(DT)
library(markdown)
library(survival)
library(survminer)
library(survMisc)

source(file.path("helpers", "ggcoxzph.R"))

source(file.path("helpers", "nonpara_helpers.R"))
source(file.path("helpers", "semipara_helpers.R"))
source(file.path("helpers", "para_helpers.R"))

source(file.path("modules", "nonpara_tableplot_module.R"))
source(file.path("modules", "stratify_module.R"))
source(file.path("modules", "nonpara_module.R"))
source(file.path("modules", "semipara_module.R"))
source(file.path("modules", "para_module.R"))

button_style1 = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
button_style2 = "color: #fff; background-color: #ff9e42; border-color: #ff7f42"

options(shiny.maxRequestSize = 300 * 1024^2)

# 注意getwd() 和 该文件路径
# source("app/富集分析.R")


header <- dashboardHeader(title = "")

# sidebar ----
sidebar <- dashboardSidebar(
  tags$style(".shiny-input-container {margin-bottom: 0px}
                #file1_progress {margin-bottom: 0px}
                .checkbox {margin-top: 0px}"),
  # 动态菜单栏
  sidebarMenuOutput("side_menus")
)

# =================================================================================

# upload
upload_box <- box(title = "upload file", status = "primary",
                  solidHeader = TRUE, width = 12,
                  fluidRow(
                    column(6, fileInput("uploadFile",
                                        label = "Choose  upload file",
                                        multiple = FALSE),
                           # todo 设置上传文件的类型
                           # accept = c(".csv", ".tsv")
                    )
                  ),
                  # fluidRow(
                  #   column(6,
                  #          tags$strong("file information")
                  #   )
                  # ),
                  fluidRow(
                    column(6,
                           tableOutput("files")
                    )
                  )
)

# download
download_box <- box(title = "download file", status = "primary",
                    solidHeader = TRUE, width = 12,
                    fluidRow(
                      column(6,
                             uiOutput("chooseFile")
                      )),
                    fluidRow(
                      column(6,
                             downloadButton("downloadData", "Download")
                      ))
)

# ================================================================

# body ----
body <- dashboardBody(

  tabItems(
    tabItem(tabName = "home_page", fluidRow(
      column(11,
             includeMarkdown("README.md"),
             offset = 1)
    )),
    # tabItem(tabName = "prediction", fluidRow(
    #   prediction_box
    # )),
    tabItem(tabName = "download", fluidRow(
      download_box
    )),
    tabItem(tabName = "upload", fluidRow(
      upload_box
    )),
    # tabItem(tabName = "search", fluidRow(
    #   search_box
    # )),
    # tabItem(tabName = "help", fluidRow()),
    # tabItem(tabName = "Home",
    #         fluidRow(
    #           column(11,
    #                  includeMarkdown("README.md"),
    #                  offset = 1)
    #         )),

    tabItem(

      tabName = "Data",
      sidebarLayout(

        sidebarPanel(

          h4(strong("Upload Data", style = "color: steelblue")),

          fileInput("file1", HTML("Choose a CSV File: "),
                    accept = c("text/csv", ".csv", "text/comma-separated-values,text/plain"),
                    placeholder = "Example: ovarian {survival}"),

          column(12, align = "right", checkboxInput("header", "Header", TRUE)),

          radioButtons("separator", "Separator: ",
                       choices = c(',', ';', Tab = '\t'),
                       selected = ",", inline = TRUE),

          br(),

          selectInput("var_time", "Time Variable: ", ""),
          selectInput("var_event", "Event Variable: ", ""),

          width = 3
        ),

        mainPanel(
          DT::dataTableOutput("DataTable"))
      )
    ),

    tabItem(tabName = "NonPara", nonPara_UI("nonPara_Panel")),
    tabItem(tabName = "SemiPara", semiPara_UI("semiPara_Panel")),
    tabItem(tabName = "Para", Para_UI("Para_Panel")),

    # ====================================================================
    # 差异分析
    tabItem(tabName = "VariationAnalysisDataUpload", fluidRow(
      sidebarPanel(

        fileInput("file", "Upload Expression File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        # Input: Checkbox if file has header ----
        checkboxInput("header11", "Header", TRUE),

        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ",")


      ),

      #Output the uploaded file
      mainPanel(

        h4("Table Header"),
        tableOutput("table")
      )
    )),

    tabItem(tabName = "bbb", fluidRow(
      sidebarPanel(

        #Button to select the type of normalization method used
        radioButtons("norm", "Normalization Method",
                     choices = c("Raw data" = "0",
                                 "Log2" = "1",
                                 "Quantile" = "2",
                                 "1+Log2" = "3"),
                     selected = "0"),
        textInput("boxtitle", "Title:", "Boxplot"),
        textInput("xbox", "X-axis:", "Samples"),
        textInput("ybox", "Y-axis:", "Expression value"),
        downloadButton("normdata", "Download normalized data")
      ),

      #Output boxplot for the normalized method used
      mainPanel(
        plotOutput("boxplot", width = "100%", height = 600)
      )
    )),
    tabItem(tabName = "ccc", fluidRow(
      sidebarPanel(
        conditionalPanel(condition = "input.tabselected==1",
                         tags$h4("Differential Expression using Limma"),
                         tags$hr(),

                         #Upload Pheno file for class/time variable
                         fileInput("file11", "Upload Pheno File",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv")),

                         # Input: Checkbox if file has header ----
                         checkboxInput("header1", "Header", TRUE),

                         # Input: Select separator ----
                         radioButtons("sep1", "Separator",
                                      choices = c(Comma = ",",
                                                  Semicolon = ";",
                                                  Tab = "\t"),
                                      selected = ","),
                         #Title for the plot
                         textInput("voltitle", "Title:", "Volcano plot"),
                         #Download button for normalized data
                         downloadButton("voldata", "Download D.E. results")
        ),

        conditionalPanel(condition = "input.tabselected==2",
                         tags$h4("Heatmap"),
                         tags$hr(),
                         textInput("lfc", "LogFC cutoff", "2"),

                         textInput("heattitle1", "Title of Heatmap:", "Tumor vs Normal"),

                         selectInput("color11", "Choose color:",
                                     choices = c("RedYellowGreen", "GreenBlue", "RedYellowBlue", "RedBlue"))
        )
      ),

      #Output Volcano plot and Heatmap
      mainPanel(
        #tableOutput("contrast"),
        tabsetPanel(type = "tabs",
                    tabPanel("Volcano", value = 1, plotOutput("volcanoplot", width = "100%", height = 600)),
                    tabPanel("Heatmap", value = 2, plotOutput("volheatmap", width = "100%", height = 600)),
                    id = "tabselected"
        )
      )
    )),
    tabItem(tabName = "ddd", fluidRow(
      sidebarPanel(
        textInput("heattitle", "Title of Heatmap:", "Top 500 most variable genes across samples"),
        textInput("vargenes", "Top __ variable genes:", "500"),
        selectInput("color11", "Choose color:",
                    choices = c("RedYellowGreen", "GreenBlue", "RedYellowBlue", "RedBlue")),
        downloadButton("vardata", "Download top variance data")

      ),
      #Output heatmap
      mainPanel(
        plotOutput("heatmap", width = "100%", height = 600)
      )
    )),

    # 富集分析
    tabItem(tabName = "eee", fluidRow(
      sidebarPanel(
        fileInput(inputId = "genelist",
                  label = h3("Gene list input (csv format):"),
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        h5("tips: The first column is gene symbol, and the second column is logFC"),
        br(),
        br(),

        # radioButtons(inputId = "species",
        #              label = "Choose",
        #              choices = c("Human"= "hsa", "Mouse"="mmu")),

        selectInput(inputId = "species",
                    label = "Choose a species:",
                    choices = c("Human" = "hsa")),
        br(),
        h2("Parameters for analysis"),
        sliderInput(inputId = "pval",
                    label = "KEGG p value cutoff",
                    min = 0, max = 1,
                    value = 0.05),
        sliderInput(inputId = "qval",
                    label = "KEGG q value cutoff",
                    min = 0, max = 1,
                    value = 0.05),
        selectInput(inputId = "padjMethod",
                    label = "Choose a p adjustment method:",
                    choices = c("None" = "none", "Benjamini-Hochberg" = "BH")),
        actionButton("run", "Run!"), ## input won't be analyzed immediately
        br(),
        h2("Parameters for plotting"),
        selectInput(inputId = "x",
                    label = "Choose X axis:",
                    choices = c("GeneRatio" = "GeneRatio", "Count" = "Count")),
        selectInput(inputId = "color",
                    label = "Color by:",
                    choices = c("Adjusted P Value" = "p.adjust",
                                "P Value" = "pvalue",
                                "q Value" = "qvalue")),
        sliderInput(inputId = "show",
                    label = "Number of pathways to show:",
                    min = 0, max = 100,
                    value = 10),
        br(),
        selectInput(inputId = "format",
                    label = "Choose File format to save:",
                    choices = c("pdf" = "pdf", "png" = "png", "svg" = "svg"))
      ),
      mainPanel(
        DT::dataTableOutput("EADataTable")
      )

    )),
    tabItem(tabName = "fff",
            fluidRow(
              tabBox(
                title = tagList(icon("gear"), "Enrichment Analysis"),
                width = 12,

                tabPanel("DotPlot",
                         fluidRow(

                           column(12,
                                  # align = "center",
                                  tags$h1("KEGG"),
                                  plotOutput("keggdotplot"),
                                  downloadButton("keggdot_down", "Download the plot"),
                                  br(),
                                  tags$h1("Reactome"),
                                  plotOutput("reactomedotplot"),
                                  downloadButton("reactomedot_down", "Download the plot"),
                                  br()
                           )
                         )
                ),
                tabPanel("Table",
                         fluidRow(

                           column(12,
                                  # align = "center",
                                  tags$h1("KEGG"),
                                  dataTableOutput("kegg_dt"),
                                  downloadButton("keggdf", "Download the table"),
                                  br(),
                                  tags$h1("Reactome"),
                                  dataTableOutput("reactome_dt"),
                                  downloadButton("reactomedf", "Download the table"),
                                  br(),
                                  verbatimTextOutput("message_1"),
                                  br()
                           )
                         )
                ),
                tabPanel("BarPlot",
                         fluidRow(

                           column(12,
                                  # align = "center",
                                  tags$h1("KEGG"),
                                  plotOutput("keggbarplot"),
                                  downloadButton("keggbar_down", "Download the plot"),
                                  br(),
                                  tags$h1("Reactome"),
                                  plotOutput("reactomebarplot"),
                                  downloadButton("reactomebar_down", "Download the plot"),
                                  br()
                           )
                         )
                ),
                tabPanel("Network",
                         fluidRow(

                           column(12,
                                  # align = "center",
                                  tags$h1("KEGG"),
                                  plotOutput("keggcnetplot"),
                                  downloadButton("keggcnet_down", "Download the plot"),
                                  br(),
                                  tags$h1("Reactome"),
                                  plotOutput("reactomecnetplot"),
                                  downloadButton("reactomecnet_down", "Download the plot"),
                                  br()
                           )
                         )
                ),
              )
            )
    )
  )
)

# assemble UI
ui <- dashboardPage(header, sidebar, body, skin = "blue")

server <- function(input, output, session) {

  output$side_menus <- renderMenu({
    sidebarMenu(
      menuItem("Home", tabName = "home_page", icon = icon("home")),
      menuItem("VariationAnalysis", tabName = "VariationAnalysis", icon = icon("calendar"),
               menuSubItem("Data Upload",
                           tabName = "VariationAnalysisDataUpload",
                           icon = icon("angle-right")),
               menuSubItem("Normalization",
                           tabName = "bbb",
                           icon = icon("angle-right")),
               menuSubItem("DifferentialExpression",
                           tabName = "ccc",
                           icon = icon("angle-right")),
               menuSubItem("Variance",
                           tabName = "ddd",
                           icon = icon("angle-right"))
      ),
      menuItem("EnrichmentAnalysis", tabName = "EnrichmentAnalysis", icon = icon("list-alt"),
               menuSubItem("Data Upload",
                           tabName = "eee",
                           icon = icon("angle-right")),
               menuSubItem("Analysis",
                           tabName = "fff",
                           icon = icon("angle-right"))
      ),
      menuItem("SurvivalAnalysis", tabName = "SurvivalAnalysis", icon = icon("table"),
               # todo 生存分析home
               # menuSubItem("Home",
               #             tabName = "Home",
               #             icon = icon("home")),

               menuSubItem("Data Upload",
                           tabName = "Data",
                           icon = icon("angle-right")),
               menuSubItem("Non-Parametric Methods",
                           tabName = "NonPara",
                           icon = icon("angle-right")),

               menuSubItem("Semi-Parametric Methods",
                           tabName = "SemiPara",
                           icon = icon("angle-right")),

               menuSubItem("Parametric Methods",
                           tabName = "Para",
                           icon = icon("angle-right"))
      ),
      # menuItem("Prediction", tabName = "prediction", icon = icon("hourglass-start")),
      menuItem("Download", tabName = "download", icon = icon("folder-open")),
      menuItem("Upload", tabName = "upload", icon = icon("exchange-alt"))
      # menuItem("Search", tabName = "search", icon = icon("globe")),
      # menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  })
  # ==============================================================================================
  # 差异分析
  #Read Expresssion input from user
  df1 <- reactive({
    req(input$file)
    df1 <- read.csv(input$file$datapath,
                    header = input$header11,
                    sep = input$sep,
                    row.names = 1,
                    check.names = F)

    return(df1)
  })

  #Read Pheno input from user
  pheno <- reactive({
    req(input$file11)
    pheno <- read.csv(input$file11$datapath,
                      header = TRUE,
                      sep = input$sep1,

                      check.names = F)

    return(pheno)
  })


  #Ouput the head of Expression file for display

  # VAData <- reactive({
  #   inFile <- input$file
  #
  #   if (!is.null(inFile)) {
  #
  #     read.csv(inFile$datapath,
  #              header = input$header11,
  #              sep = input$sep1)
  #   }
  #   else {  read.csv(file.path("data", "VariationAnalysis/VariationAnalysis-Expressionfile-test.csv")) }
  # }
  # )

  output$table <- renderTable({
    req(input$file)
    df1 <- read.csv(input$file$datapath,
                    header = input$header11,
                    sep = input$sep,

                    check.names = F)

    return(head(df1))

    # DT::datatable(VAData(), options = list(scrollX = TRUE))
  })

  #Output example dataset to display on home page
  #Iris data from library(DT)
  output$mytable3 <- DT::renderDataTable({
    DT::datatable(iris[, -5], colnames = c('Genes', 'Sample1', 'Sample2', 'Sample3', 'Sample4'), options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })

  #select normalization method based on user selection
  norm <- reactive({
    if (input$norm == "0") {
      #raw data
      norm <- df1()
    }
    if (input$norm == "1") {
      #log2 transform
      norm <- log2(df1())
    }
    if (input$norm == "2") {
      #Quantile normalization
      norm <- normalize.quantiles(as.matrix(df1()))
      colnames(norm) <- colnames(df1())
    }
    if (input$norm == "3") {
      # 1+log2 transform
      norm <- log2(df1() + 1)
    }
    rownames(norm) <- rownames(df1())
    return(norm)
  })

  #Print the normalized data to a file named "normdata.csv"
  output$normdata <- downloadHandler(
    filename = function() {
      paste("normdata.csv", sep = "")
    },
    content = function(file) {
      write.csv(norm(), file, row.names = TRUE)
    }
  )

  #Function to output boxplot based on the normalization selected
  plotBox <- function() {
    boxplot(norm(), main = input$boxtitle, xlab = input$xbox, ylab = input$ybox) +
      geom_point()

  }

  #Output the boxplot generated above
  output$boxplot <- renderPlot({
    plotBox()
  })

  #Differential Expression using Limma package
  #Function to output the model fit for volcano plot
  volcano <- function() {
    group <- factor(pheno()[, 2])
    design <- model.matrix(~0 + group)
    ## Make the column names of the design matrix a bit nicer
    colnames(design) <- levels(group)
    fit <- lmFit(norm(), design)
    fit.cont <- contrasts.fit(fit, c(-1, 1))
    fit.cont <- eBayes(fit.cont)
    return(fit.cont)
  }

  #Output volcano plot from the result generated above
  output$volcanoplot <- renderPlot({
    volcanoplot(volcano(), coef = 1, main = input$voltitle)
  })

  #print the results from Differential expression analysis to a file named "Results.csv"
  output$voldata <- downloadHandler(
    filename = function() {
      paste("Results.csv", sep = "")
    },
    content = function(file) {
      write.csv(volcano(), file, quote = F, row.names = TRUE)
    }
  )

  #Take the LogFC filter from user and generate heatmap based on the D.E. genes
  output$volheatmap <- renderPlot({
    probeset.list1 <- topTable(volcano(), coef = 1, n = "Inf", lfc = input$lfc)
    genes <- row.names(probeset.list1)
    norm.data <- norm()[genes,]
    heatmap.2(as.matrix(norm.data), col = rev(color()(50)), trace = "none", main = input$heattitle1,
              scale = "row", Colv = "NA")
  })

  #Function to select number of genes with most variation
  heat <- reactive({
    logcounts <- norm()
    rownames(logcounts) <- rownames(df1())
    var_genes <- apply(logcounts, 1, var)

    # Get the gene names for the top 500 most variable genes
    select_var <- names(sort(var_genes, decreasing = TRUE))[1:input$vargenes]

    # Subset logcounts matrix
    highly_variable_lcpm <- logcounts[select_var,]
    return(highly_variable_lcpm)
  })

  #Variable to select different colors by user
  color <- reactive({
    if (input$color11 == "RedYellowGreen") {
      mypalette <- brewer.pal(11, "RdYlGn")
      morecols <- colorRampPalette(mypalette)
    }
    if (input$color11 == "GreenBlue") {
      morecols <- colorRampPalette(c("blue", "white", "green"))
    }
    if (input$color11 == "RedYellowBlue") {
      mypalette <- brewer.pal(11, "RdYlBu")
      morecols <- colorRampPalette(mypalette)
    }
    if (input$color11 == "RedBlue") {
      morecols <- colorRampPalette(c("blue", "white", "red"))
    }
    return(morecols)
  })

  #Generate heatmap for the most variable genes generated by function heat
  output$heatmap <- renderPlot({
    heatmap.2(as.matrix(heat()), col = rev(color()(50)), trace = "none", main = input$heattitle,
              scale = "row", Colv = "NA")

  })

  #Print the gene expression of most variable genes to a file named "variance.csv"
  output$vardata <- downloadHandler(
    filename = function() {
      paste("variance.csv", sep = "")
    },
    content = function(file) {
      write.csv(heat(), file, row.names = TRUE)
    }
  )


  # ==============================================================================================
  # 富集分析
  ## counters
  output$num <- renderText({
    if (!file.exists("data/counter.RData")) {
      counter <- 0
    } else
      load(file = "data/counter.RData")
    counter <- counter + 1
    save(counter, file = "data/counter.RData")
    paste0("You are the ", counter, " visitors.")
  })

  # todo  Ea
  EAData <- reactive({
    inFile <- input$genelist

    if (!is.null(inFile)) {

      read_csv(input$genelist$datapath, col_names = F)
    }
    else {
      read.csv(file.path("data", "EnrichmentAnalysis-test.csv"))
    }

  })

  df <- eventReactive(input$run, {
    shiny::validate(
      need(input$genelist, "Upload a file!")
    )
    df <- read_csv(input$genelist$datapath, col_names = F)

    # if (is.null(df)) {
    #   df <- EAData()
    # }

    return(df)
  })

  output$EADataTable <- DT::renderDataTable(
    DT::datatable(EAData(), options = list(scrollX = TRUE)))

  org <- eventReactive(input$run, {
    org <- ifelse(input$species == "hsa", "org.Hs.eg.db", "")
  })

  id <- eventReactive(input$run, {
    id <- bitr(toupper(df()$X1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org())
    return(id)
  })

  fc <- eventReactive(input$run, {
    if (ncol(df()) == 2) {
      fc <- df()$X2[match(id()$SYMBOL, df()$X1)]
      names(fc) <- id()$ENTREZID
    } else
      fc <- NULL
    return(fc)

  })

  kegg <- eventReactive(input$run, {
    kegg <- enrichKEGG(id()$ENTREZID, organism = input$species,
                       pvalueCutoff = input$pval,
                       pAdjustMethod = input$padjMethod,
                       minGSSize = 10,
                       maxGSSize = 500,
                       qvalueCutoff = input$qval) %>%
      setReadable(., OrgDb = org(), keyType = "ENTREZID")
    return(kegg)
  })

  reactome_sp <- eventReactive(input$run, {
    reactome_sp <- ifelse(input$species == "hsa", "human", "mouse")
  })

  reactome <- eventReactive(input$run, {
    reactome <- enrichPathway(id()$ENTREZID, organism = reactome_sp(),
                              pvalueCutoff = input$pval,
                              pAdjustMethod = input$padjMethod,
                              minGSSize = 10,
                              maxGSSize = 500,
                              qvalueCutoff = input$qval,
                              readable = T)
    return(reactome)
  })

  # output$kegg_dt <- renderDataTable(fc() %>% data.frame())
  output$kegg_dt <- renderDataTable({
    kegg() %>%
      data.frame() %>%
      datatable() %>%
      formatRound(columns = 5:7, digits = 3)
  })

  output$reactome_dt <- renderDataTable({
    reactome() %>%
      data.frame() %>%
      datatable() %>%
      formatRound(columns = 5:7, digits = 3)
  })

  output$message_1 <- renderText({
    paste0(nrow(id()), " out of ", nrow(df()), " genes are mapped!")
  })

  # 'length(x) = 3 > 1' in coercion to 'logical(1)'
  output$keggdotplot <- renderPlot({
    dotplot(kegg(), x = input$x, color = input$color,
            showCategory = input$show)
  })

  output$reactomedotplot <- renderPlot({
    dotplot(reactome(), x = input$x, color = input$color,
            showCategory = input$show)
  })

  output$keggbarplot <- renderPlot({
    barplot(kegg(), x = input$x, color = input$color,
            showCategory = input$show) + ylab(input$x)
  })

  output$reactomebarplot <- renderPlot({
    barplot(reactome(), x = input$x, color = input$color,
            showCategory = input$show) + ylab(input$x)
  })

  output$keggcnetplot <- renderPlot({
    cnetplot(kegg(), # x=input$x, color=input$color,
             foldChange = fc(),
             showCategory = input$show)
  })

  output$reactomecnetplot <- renderPlot({
    cnetplot(reactome(), # x=input$x, color=input$color,
             foldChange = fc(),
             showCategory = input$show)
  })


  ## from here, create download files
  output$keggdot_down <- downloadHandler(
    ## specify file name
    filename = function() {
      paste0("keggdotplot.", input$format)
    },
    content = function(file) {
      ggsave(filename = file, device = input$format,
             plot = dotplot(kegg(), x = input$x,
                            color = input$color,
                            showCategory = input$show),
             width = 10, height = 5)
    })

  output$reactomedot_down <- downloadHandler(
    ## specify file name
    filename = function() {
      paste0("reactomedotplot.", input$format)
    },
    content = function(file) {
      ggsave(filename = file, device = input$format,
             plot = dotplot(reactome(), x = input$x,
                            color = input$color,
                            showCategory = input$show),
             width = 10, height = 5)
    })

  output$keggbar_down <- downloadHandler(
    ## specify file name
    filename = function() {
      paste0("keggbarplot.", input$format)
    },
    content = function(file) {
      ggsave(filename = file, device = input$format,
             plot = barplot(kegg(), x = input$x,
                            color = input$color,
                            showCategory = input$show) + ylab(input$x),
             width = 10, height = 5)
    })

  output$reactomebar_down <- downloadHandler(
    ## specify file name
    filename = function() {
      paste0("reactomedbarplot.", input$format)
    },
    content = function(file) {
      ggsave(filename = file, device = input$format,
             plot = barplot(reactome(), x = input$x,
                            color = input$color,
                            showCategory = input$show) + ylab(input$x),
             width = 10, height = 5)
    })

  output$keggcnet_down <- downloadHandler(
    ## specify file name
    filename = function() {
      paste0("keggcnetplot.", input$format)
    },
    content = function(file) {
      ggsave(filename = file, device = input$format,
             plot = cnetplot(kegg(), # x=input$x, color=input$color,
                             foldChange = fc(),
                             showCategory = input$show),
             width = 10, height = 5)
    })

  output$reactomecnet_down <- downloadHandler(
    ## specify file name
    filename = function() {
      paste0("reactomecnetplot.", input$format)
    },
    content = function(file) {
      ggsave(filename = file, device = input$format,
             plot = cnetplot(reactome(), # x=input$x, color=input$color,
                             foldChange = fc(),
                             showCategory = input$show),
             width = 10, height = 5)
    })

  ## download pathway table
  output$keggdf <- downloadHandler(
    ## specify file name
    filename = "keggtable.csv",
    content = function(file) {
      write.csv(kegg() %>% data.frame(), file)
    })

  output$reactomedf <- downloadHandler(
    ## specify file name
    filename = "reactometable.csv",
    content = function(file) {
      write.csv(reactome() %>% data.frame(), file)
    })

  ## maybe create a zip file contain all?
  # output$downloadall <- downloadHandler(
  #   filename = "all.zip",
  #   content = function(file){
  #
  #   }
  # )

  output$sessionInfo <- renderText({
    paste(capture.output(sessionInfo()),
          collapse = "\n")
  })

  # ============================================================================================
  # 生存分析
  rawData <- reactive({
    inFile <- input$file1

    if (!is.null(inFile)) {

      read.csv(inFile$datapath,
               header = input$header,
               sep = input$separator)
    }
    else { read.csv(file.path("data", "ovarian.csv")) }
  }
  )

  # Time and Event Variable Selection
  observe({
    updateSelectInput(session, "var_time", label = "Time Variable",
                      choices = names(rawData()))
  })

  observeEvent(Time(), {
    updateSelectInput(session, "var_event", label = "Event Variable",
                      choices = setdiff(names(rawData()), Time()))

  })

  observeEvent(input$reset_all,
               { shinyjs::reset("separator") })

  Time = reactive({ input$var_time })
  Event = reactive({ input$var_event })

  callModule(nonPara, "nonPara_Panel", data = rawData,
             Time = Time, Event = Event)
  callModule(semiPara, "semiPara_Panel", data = rawData,
             Time = Time, Event = Event)
  callModule(Para, "Para_Panel", data = rawData,
             Time = Time, Event = Event)


  # Output----------
  output$DataTable <- DT::renderDataTable(
    DT::datatable(rawData(), options = list(scrollX = TRUE)))

  # =============================================================================================
  # 上传文件
  data_upload <- reactive({
    # 上传的文件在该目录下
    # C:/Users/Administrator/AppData/Local/Temp/RtmpABfE8H
    # 上传文件的路径input$uploadFile$datapath

    inFile <- input$uploadFile

    if (is.null(inFile))
      return(NULL)

    # print(inFile$name)
    # file.copy(inFile$datapath, paste("app", inFile$name), sep = "/")
    print("upload success")
  })

  output$files <- renderTable(input$uploadFile)

  # ==================================================================================
  # 下载文件  chooseFile
  output$chooseFile <- renderUI({
    # todo  下载文件路径
    selectInput(inputId = "chooseFile", label = "Choose download file: ", choices = list.files("data"), selected = 1)
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      # paste("output", "zip", sep=".")
      input$chooseFile
    },
    content = function(file) {
      # C:/Users/Administrator/AppData/Local/Temp/RtmpABfE8H
      # print(tempdir())   #可以将下载的文件放到临时目录
      # vroom::vroom_write(datasetInput(), file)

      # todo 下载文件路径
      # file.copy("test/out.zip", file)
      file.copy(paste("test", input$chooseFile, sep = "/"), file)
    },
    # contentType = "application/zip"
  )

}

# ==================================================================================

shinyApp(ui = ui, server = server)