#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

library(shinyMatrix)
library(rhandsontable)



# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Markov Model Comparison Tool"),
  
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "nstates",
                   label = "Number of states in the model:",
                   value = 3,min=1)
    ,
    
      numericInput(inputId = "oldCycle",
                 label = "Cycle length of original model",
                 value = 12,min=1)
      ,
  
  
      numericInput(inputId = "newCycle",
               label = "Cycle length of converted model (in same units)",
               value = 1,min=1)
      ,
  
    numericInput(inputId = "ncycles",
               label = "Number of cycles to run (using original cycle length)",
               value = 10,min=1),
    
    h5("Original Transition Matrix (A)"),
    rHandsontableOutput("oldM"),
    h5("Converted Transition Matrix (B)"),
    rHandsontableOutput("newM"),
    h5("Initial Distribution of Cohort"),
    rHandsontableOutput("initD")
    #actionButton("calculate","Calculate")
      ),
    
    
  
  mainPanel(
tabsetPanel(
    tabPanel("Intro",
             h5(
               "The purpose of this app is to examine the difference between two Markov matrices over time.",
               textOutput("summary3")
                )
             
             
             
             ),
    tabPanel("Summary",textOutput("summary1"),
             tableOutput("convM"),
             textOutput("summary2"),
             plotOutput("errorPlot"),
             textOutput("summary4")
             ),
    tabPanel("Cohort Plot",plotOutput("cohortPlot")),
    tabPanel("State-by-State plot",plotOutput("plots")),
    tabPanel("Detailed Table",rHandsontableOutput("MarkovDist"))
  
)
  
 )
  
  
  
  )
))
