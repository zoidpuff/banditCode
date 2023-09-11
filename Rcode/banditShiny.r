# Load helper functions
source("banditAnalysis.r")
library(rhdf5)
library(tidyverse)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(ggplot2)

datafile <- '/home/gummi/banditExperiment/data/dataset6012023.data'

# Load the data
datafile <- '/home/gummi/banditExperiment/data/dataset6012023.data'
data <- h5dump(datafile,load=TRUE)

# Compute experiment stats
experimentData <- tabulateExperiments(data,datafile,2)

# Create different datasets
datasets <- list(
  Latter = experimentData[["SplitProbExperiments"]][[2]],
  Former = experimentData[["SplitProbExperiments"]][[1]],
  Entire = experimentData[["ProbExperimentsProbPart"]]
)

safe_names <- make.names(names(datasets[[1]]))
names(safe_names) <- names(datasets[[1]])

# Make sure Mouse column is character
for (i in seq_along(datasets)) {
  datasets[[i]]$Mouse <- as.character(datasets[[i]]$Mouse)
  names(datasets[[i]]) <- make.names(names(datasets[[i]]))  
}

# Define UI for application
ui <- fluidPage(theme = shinytheme("darkly"),
   
   # Application title
   titlePanel("Interactive Barplot of Experiment Data"),
   
   # Sidebar with controls
   sidebarLayout(
      sidebarPanel(
         selectInput("dataset", 
                     "Choose a dataset:", 
                     choices = names(datasets)),
         selectInput("variable", 
                     "Choose a variable to plot:", 
                     choices = names(safe_names)[sapply(datasets[[1]], is.numeric)]), 
         selectInput("fill", 
                     "Choose a fill variable:", 
                     choices = names(safe_names)[sapply(datasets[[1]], is.character)]), 
         selectInput("facet", 
                     "Choose a variable to facet by:", 
                     choices = names(safe_names)[sapply(datasets[[1]], is.character)])
      ),
      
      # Main panel for displaying plots
      mainPanel(
         plotOutput("barplot", height = 900)
      )
   )
)

# Define server logic 
server <- function(input, output) {

   # Reactive expression to get the selected dataset
   selectedData <- reactive({
      datasets[[input$dataset]]
   })

   output$barplot <- renderPlot({
      
      # Generate the plot with the selected dataset
      p <- ggplot(selectedData(), aes_string(x = safe_names[input$fill], y = safe_names[input$variable], fill = safe_names[input$fill])) +
         geom_boxplot() +
         geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
         theme_classic() +
         labs(y = input$variable, title = "Experiment Data") 
      
      # If a facet variable is chosen, add facet_wrap
      if (!is.null(input$facet) && input$facet != "") {
         p <- p + facet_wrap(as.formula(paste("~", safe_names[input$facet])))
      }
      
      print(p)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
