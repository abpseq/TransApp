library(shiny)
library(shinyIncubator)
#source("trans_pipe.R")
# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  # Application title
  titlePanel("Transcriptomic Data Analysis"),
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    radioButtons("platform", "Select your Array:",
                 c(#"None" = 0,
                   "Affymetrix Human Genome U133 Plus 2.0 Array" = 1,
                   "Illumina HumanHT-12 V4.0 expression beadchip" = 2,
                   "Customised" = 3)),
    
    textInput("proj_id", "Enter Project ID:", "GSEXXXX"),
    submitButton("Run"),
    tags$div(class="header", checked=NA,
             tags$p(""),
             tags$p(""),
             tags$a(href="TransApp_help.pdf", h6("Help",align ="right"))
    )
  ),
  
  mainPanel(  progressInit(),
              img(src = "Flowchart.png", height = 300, width = 500),
              #   verbatimTextOutput("info")
              verbatimTextOutput("trans_pipe")
  )
))
