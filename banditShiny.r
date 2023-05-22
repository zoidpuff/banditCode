# Shiny
# Load helper functions
source("banditAnalysis.r")
library(rhdf5)
library(tidyverse)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(ggplot2)

datafile <- '/home/gummi/banditExperiment/dataset.data'

# Load the hdf5 file
data <- h5dump(datafile,load=TRUE)

# Compute experiment stats
experimentData <- tabulateExperiments(data,datafile,2)
experimentData <- experimentData[["ProbExperimentsProbPart"]]
experimentData$Mouse <- as.character(experimentData$Mouse)
names(experimentData) <- make.names(names(experimentData))


# Define UI for application
ui <- fluidPage(theme = shinytheme("darkly"),
   
   # Application title
   titlePanel("Interactive Barplot of Experiment Data"),
   
   # Sidebar with controls
   sidebarLayout(
      sidebarPanel(
         selectInput("variable", 
                     "Choose a variable to plot:", 
                     choices = colnames(experimentData)[sapply(experimentData, is.numeric)]), 
         selectInput("fill", 
                     "Choose a fill variable:", 
                     choices = colnames(experimentData)[sapply(experimentData, is.character)]), 
         selectInput("facet", 
                     "Choose a variable to facet by:", 
                     choices = colnames(experimentData)[sapply(experimentData, is.character)])
      ),
      
      # Main panel for displaying plots
      mainPanel(
         plotOutput("barplot", height = 900)
      )
   )
)

# Define server logic 
server <- function(input, output) {
   
   output$barplot <- renderPlot({
      
      # Generate the plot
      p <- ggplot(experimentData, aes_string(x = input$fill, y = input$variable, fill = input$fill)) +
         geom_boxplot() +
         theme_classic() +
         labs(y = input$variable, title = "Experiment Data") 
      
      # If a facet variable is chosen, add facet_wrap
      if (!is.null(input$facet) && input$facet != "") {
         p <- p + facet_wrap(as.formula(paste("~", input$facet)))
      }
      
      print(p)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
