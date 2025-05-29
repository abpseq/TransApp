# Install and load required packages with checks
required_packages <- c(
  "shiny", "shinydashboard", "shinyWidgets", "DT", "plotly",
  "GEOquery", "genefilter", "GSEABase", "KEGG.db", "Category",
  "lattice", "SPIA", "annotate", "GOstats", "RColorBrewer",
  "xtable", "Rgraphviz", "hgu133plus2.db", "hgu133plus2cdf",
  "hgu133plus2probe", "illuminaHumanv4.db", "limma", "Biobase",
  "gcrma", "arrayQualityMetrics", "affy", "ReportingTools",
  "hwriter", "KEGGprofile", "DOSE", "shinyjs"
)

# Function to check and install packages
install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Load all required packages
lapply(required_packages, install_and_load)

# Utility function to trim whitespace
trim <- function(x) gsub("^\\s+|\\s+$", "", x)

# Define UI with modern dashboard and custom CSS/JavaScript
ui <- dashboardPage(
  # Dashboard header
  dashboardHeader(title = "TransApp: Transcriptomic Analysis"),
  # Sidebar with menu and inputs
  dashboardSidebar(
    sidebarMenu(
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-bar")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    ),
    selectInput("platform", "Select Array Platform:",
                choices = c("Affymetrix Human Genome U133 Plus 2.0 Array" = 1,
                            "Illumina HumanHT-12 V4.0 expression beadchip" = 2,
                            "Customised" = 3),
                selected = 1),
    textInput("proj_id", "Enter Project ID:", "GSEXXXX"),
    fileInput("targets_file", "Upload Targets File (Optional)", accept = ".txt"),
    uiOutput("group_selection_ui"), # Dynamic UI for group selection
    actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
    downloadButton("download_report", "Download Report"),
    # Custom CSS for styling
    tags$head(
      tags$style(HTML("
        .skin-blue .main-header .logo { 
          background-color: #1a3c6d; 
          font-weight: bold; 
          color: white;
        }
        .skin-blue .main-sidebar { 
          background-color: #2a4b8d; 
        }
        .box { 
          border-radius: 8px; 
          box-shadow: 0 4px 10px rgba(0,0,0,0.2); 
          background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        }
        .btn-primary { 
          background-color: #1a3c6d; 
          border-color: #1a3c6d; 
          transition: all 0.3s ease;
        }
        .btn-primary:hover { 
          background-color: #2a4b8d; 
          transform: scale(1.05); 
        }
        .content-wrapper { 
          background-color: #f4f6f9; 
        }
        .shiny-progress .progress-text { 
          font-size: 16px; 
          color: #1a3c6d; 
        }
      ")),
      # JavaScript for animations and interactivity
      tags$script(HTML("
        $(document).ready(function() {
          // Fade-in effect for boxes
          $('.box').each(function(index) {
            $(this).delay(200 * index).fadeIn(1000);
          });
          // Hover effect for buttons
          $('.btn').hover(
            function() { $(this).css('cursor', 'pointer'); },
            function() { $(this).css('cursor', 'default'); }
          );
        });
      "))
    )
  ),
  # Main dashboard body
  dashboardBody(
    useShinyjs(),
    tabItems(
      tabItem(tabName = "analysis",
              fluidRow(
                box(title = "Analysis Flowchart", status = "primary", solidHeader = TRUE,
                    img(src = "Flowchart.png", height = 300, width = 500)),
                box(title = "Analysis Progress", status = "info", solidHeader = TRUE,
                    verbatimTextOutput("trans_pipe")),
                box(title = "Top Differentially Expressed Genes", status = "success", solidHeader = TRUE,
                    DTOutput("top_genes_table")),
                box(title = "Heatmap", status = "warning", solidHeader = TRUE,
                    plotlyOutput("heatmap_plot", height = "500px"))
              )),
      tabItem(tabName = "help",
              box(title = "Help", status = "primary", solidHeader = TRUE,
                  tags$a(href = "TransApp_help.pdf", "Download Help Document", class = "btn btn-info")))
    )
  )
)

# Define Server
server <- function(input, output, session) {
  options(digits = 3)
  
  # Reactive value to store expression set
  eset_data <- reactiveVal()
  
  # Dynamic group selection UI (GEO2R-like)
  output$group_selection_ui <- renderUI({
    req(input$targets_file)
    targets <- read.table(input$targets_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if ("Group" %in% colnames(targets)) {
      unique_groups <- unique(targets$Group)
      tagList(
        selectInput("group1", "Select Case Group:", choices = unique_groups, selected = unique_groups[1]),
        selectInput("group2", "Select Control Group:", choices = unique_groups, selected = unique_groups[2])
      )
    } else {
      showNotification("Targets file must contain a 'Group' column.", type = "error")
      NULL
    }
  })
  
  # Main analysis logic
  observeEvent(input$run_analysis, {
    req(input$platform, input$proj_id != "GSEXXXX", input$targets_file)
    techid <- as.numeric(input$platform)
    projid <- input$proj_id
    
    # Progress bar for package loading
    withProgress(message = "Loading required packages", value = 0, {
      for (i in 1:10) {
        setProgress(value = i / 10)
        Sys.sleep(0.2)
      }
    })
    
    # Define file paths
    rawdatapath <- paste(projid, "/raw_data/", sep = "")
    targetspath <- input$targets_file$datapath
    processeddatapath <- paste(projid, "/processed_data/", sep = "")
    resultpath <- paste(projid, "/result/", sep = "")
    
    # Create directories if they don't exist
    dir.create(processeddatapath, recursive = TRUE, showWarnings = FALSE)
    dir.create(resultpath, recursive = TRUE, showWarnings = FALSE)
    
    # Read targets file
    targets <- read.table(targetspath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Process data based on platform
    eset.rda_file_name <- paste(processeddatapath, projid, "_eset.rda", sep = "")
    deresultfile <- paste(resultpath, paste(projid, "_All_genes.csv", sep = ""), sep = "")
    
    if (file.exists(eset.rda_file_name)) {
      # Load cached expression set
      load(file = eset.rda_file_name)
      eset_data(eset)
      output$trans_pipe <- renderText({paste("Loaded cached data for", projid)})
    } else {
      withProgress(message = "Processing data", value = 0, {
        if (techid == 1) {
          # Affymetrix processing
          myCovs <- data.frame(targets)
          rownames(myCovs) <- myCovs[, 1]
          metadata <- data.frame(labelDescription = paste(colnames(myCovs), ": ", 
                                                          apply(myCovs, 2, function(x) nlevels(as.factor(x))), 
                                                          " level(s)", sep = ""),
                                 row.names = colnames(myCovs))
          phenoData <- new("AnnotatedDataFrame", data = myCovs, varMetadata = metadata)
          dat <- ReadAffy(sampleNames = myCovs$Name, filenames = myCovs$Celfile, 
                          phenoData = phenoData, celfile.path = rawdatapath)
          eset <- gcrma(dat)
          annotation(eset) <- "hgu133plus2.db"
          eset_filt <- nsFilter(eset, var.filter = FALSE)$eset
          featureNames(eset_filt) <- as.character(unlist(mget(featureNames(eset_filt), hgu133plus2ENTREZID)))
          annotation(eset_filt) <- "org.Hs.eg.db"
          eset_data(eset_filt)
          save(eset_filt, file = eset.rda_file_name)
        } else if (techid == 2) {
          # Illumina processing (placeholder)
          output$trans_pipe <- renderText({"Illumina processing not implemented in this version."})
          return()
        } else if (techid == 3) {
          # Custom data processing
          datafile <- paste(projid, "/rawdata/", projid, "_gexp.txt", sep = "")
          feat.datafile <- paste(projid, "/metadata/", projid, "_feature_annot.csv", sep = "")
          if (!file.exists(datafile) || !file.exists(feat.datafile)) {
            showNotification("Required data or annotation files not found.", type = "error")
            return()
          }
          expr <- as.matrix(read.table(file = datafile, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE))
          feat_annot <- as.matrix(read.table(file = feat.datafile, header = TRUE, sep = "\t", as.is = TRUE))
          rownames(feat_annot) <- feat_annot[, "Probe_ID"]
          
          expr <- expr[, trim(targets$Sample.ID)]
          pData <- targets[, c("Sample.ID", "Group")]
          rownames(pData) <- trim(targets$Sample.ID)
          metadata <- data.frame(labelDescription = c("Patient ID", "Case/control status"), 
                                 row.names = c("Sample.ID", "Group"))
          adf <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)
          eset <- new("ExpressionSet", exprs = expr, phenoData = adf)
          
          nonmissing <- apply(exprs(eset), 1, function(x) sum(!is.na(x)) >= 3)
          eset <- eset[nonmissing, ]
          a.na <- feat_annot[featureNames(eset), "Entrez_ID"]
          a <- trim(a.na[!is.na(a.na) & !grepl("//", a.na)])
          if (length(a) > 0) eset <- eset[names(a), ]
          b <- apply(exprs(eset), 1, IQR, na.rm = TRUE)
          b <- b[!is.na(b)]
          names(b) <- a
          cx <- unique(a)
          maxin <- integer()
          if (length(cx) != length(a)) {
            for (i in seq_along(cx)) {
              index <- which(cx[i] == a)
              maxin[i] <- index[which.max(b[index])]
            }
          }
          if (length(maxin) > 0) {
            eset <- eset[maxin, ]
            featureNames(eset) <- a[maxin]
          }
          null.idx <- which(featureNames(eset) == "")
          if (length(null.idx) > 0) eset <- eset[-null.idx, ]
          annotation(eset) <- "org.Hs.eg.db"
          eset_data(eset)
          save(eset, file = eset.rda_file_name)
        }
        for (i in 1:10) {
          setProgress(value = i / 10)
          Sys.sleep(0.2)
        }
      })
    }
    
    # Log-transformation
    eset <- eset_data()
    if (max(exprs(eset)) > 20) {
      withProgress(message = paste("Log-transforming", projid), value = 0, {
        edat <- apply(exprs(eset), 2, function(x) {
          maxx <- max(x, na.rm = TRUE)
          if (maxx > 20) {
            minx <- min(x, na.rm = TRUE)
            if (minx < 1) x <- x + 1 - minx
            log2(x)
          } else {
            x
          }
        })
        exprs(eset) <- edat
        eset_data(eset)
        for (i in 1:10) {
          setProgress(value = i / 10)
          Sys.sleep(0.2)
        }
      })
    }
    
    # Group selection for differential expression
    if (!is.null(input$group1) && !is.null(input$group2)) {
      eset$Group <- as.character(eset$Group)
      eset$Group[eset$Group == input$group1] <- "Case"
      eset$Group[eset$Group == input$group2] <- "Zcontrol"
    } else {
      eset$Group[eset$Group != "Zcontrol"] <- "Case"
    }
    Grp <- factor(eset$Group)
    
    # Differential expression analysis
    withProgress(message = "Running differential expression analysis", value = 0, {
      des <- model.matrix(~0 + Grp)
      colnames(des) <- levels(Grp)
      f1 <- lmFit(eset, des)
      contrast.matrix <- makeContrasts(CasevsZcontrol = Case - Zcontrol, levels = des)
      ef1 <- contrasts.fit(f1, contrast.matrix)
      ef1 <- eBayes(ef1)
      ef1$Entrez_Id <- featureNames(eset)
      
      top_table1 <- topTable(ef1, number = Inf)
      top_table_all <- data.frame(
        Gene_SYMBOL = as.character(unlist(mget(featureNames(eset), org.Hs.egSYMBOL, ifnotfound = NA))),
        EntrezID = featureNames(eset),
        Log_fold_change = formatC(top_table1[, "logFC"], digits = 3),
        AveExpr = formatC(top_table1[, "AveExpr"], digits = 3),
        p_value = formatC(top_table1[, "P.Value"], digits = 3),
        adj_P_Val = formatC(top_table1[, "adj.P.Val"], digits = 3),
        Gene_name = as.character(unlist(mget(featureNames(eset), org.Hs.egGENENAME, ifnotfound = NA)))
      )
      top_table_all <- top_table_all[!is.na(top_table_all$EntrezID), ]
      write.table(top_table_all, file = deresultfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
      
      # Render top genes table
      output$top_genes_table <- renderDT({
        datatable(top_table_all[1:20, ], options = list(pageLength = 10, autoWidth = TRUE, 
                                                        dom = 'Bfrtip', buttons = c('csv', 'excel')))
      })
      
      for (i in 1:10) {
        setProgress(value = i / 10)
        Sys.sleep(0.2)
      }
    })
    
    # Interactive heatmap
    withProgress(message = "Generating heatmap", value = 0, {
      col.var <- ifelse(eset$Group == "Zcontrol", "green", "red")
      vars <- apply(exprs(eset), 1, var, na.rm = TRUE)
      which.varying <- order(vars, decreasing = TRUE)[1:500]
      heatmap_data <- exprs(eset)[which.varying, ]
      output$heatmap_plot <- renderPlotly({
        heatmaply(heatmap_data, colors = c("green", "red"), showticklabels = FALSE,
                  main = "Heatmap of Top 500 Genes by Variance")
      })
      for (i in 1:10) {
        setProgress(value = i / 10)
        Sys.sleep(0.2)
      }
    })
    
    # Generate HTML report
    withProgress(message = "Generating report", value = 0, {
      html_report <- HTMLReport(
        shortName = paste("Gene_Expression_Report_", projid, sep = ""),
        title = paste("Gene Expression Data Analysis Report for", projid),
        reportDirectory = paste("./", projid, "/report", sep = "")
      )
      publish(hwrite(paste("Project ID: ", projid, sep = ""), style = "font-size:150%", br = TRUE), html_report)
      publish(hwrite(ifelse(techid == 1, "Array: Affymetrix Human Genome U133 Plus 2.0 Array",
                            ifelse(techid == 2, "Array: Illumina HumanHT-12 V4.0 expression beadchip",
                                   "Array: Customised")),
                     style = "font-size:150%", br = TRUE), html_report)
      publish(hwrite("Top 20 Differentially Expressed Genes", style = "font-size:150%", br = TRUE), html_report)
      publish(as.data.frame(top_table_all[1:20, ]), html_report)
      finish(html_report)
      
      output$trans_pipe <- renderText({paste("Analysis completed for", projid)})
      
      # Download report
      output$download_report <- downloadHandler(
        filename = function() {
          paste("Report_", projid, ".html", sep = "")
        },
        content = function(file) {
          file.copy(paste("./", projid, "/report/Gene_Expression_Report_", projid, ".html", sep = ""), file)
        }
      )
      
      for (i in 1:10) {
        setProgress(value = i / 10)
        Sys.sleep(0.2)
      }
    })
    
    # Clean up memory
    gc()
  })
}

# Run the app
shinyApp(ui = ui, server = server)