library(shiny)
library(tidyverse)
library(shinyWidgets)
library(DT)
library(Biostrings)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .container {
        display: flex;
      }
      .item {
        flex-basis: 100px;
        margin: 5px;
      }
      "
    ))
  ),
  titlePanel("AA2Cod"), br(),
  sidebarLayout(
    div(style="width: 500px;", sidebarPanel(width=12,
                                            fluidRow(
                                              div(class="container",
                                                  div(class="item",
                                                      HTML('<strong>Hydrophobic</strong>'), br(), br(),
                                                      materialSwitch('allH', 'all', value = FALSE, status = 'primary', right = TRUE),br(),
                                                      materialSwitch('A', 'A', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('V', 'V', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('I', 'I', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('L', 'L', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('M', 'M', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('F', 'F', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('Y', 'Y', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('W', 'W', value = FALSE, status = 'primary', right = TRUE)
                                                  ),
                                                  div(class="item",
                                                      HTML('<strong>Charged</strong>'), br(), br(),
                                                      materialSwitch('allC', 'all', value = FALSE, status = 'primary', right = TRUE),br(),
                                                      materialSwitch('D', 'D', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('E', 'E', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('K', 'K', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('R', 'R', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('H', 'H', value = FALSE, status = 'primary', right = TRUE)
                                                  ),
                                                  div(class="item",
                                                      HTML('<strong>Polar</strong>'), br(), br(),
                                                      materialSwitch('allP', 'all', value = FALSE, status = 'primary', right = TRUE),br(),
                                                      materialSwitch('S', 'S', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('T', 'T', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('N', 'N', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('Q', 'Q', value = FALSE, status = 'primary', right = TRUE)
                                                  ),
                                                  div(class="item",
                                                      HTML('<strong>Special</strong>'), br(), br(),
                                                      materialSwitch('allS', 'all', value = FALSE, status = 'primary', right = TRUE),br(),
                                                      materialSwitch('G', 'G', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('C', 'C', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('P', 'P', value = FALSE, status = 'primary', right = TRUE),
                                                      materialSwitch('Stop', 'Stop', value = FALSE, status = 'primary', right = TRUE)
                                                  )
                                              ),br(),hr(), htmlOutput('selectaa'), 
                                            ))),
    
    mainPanel(
      fluidRow(column(width =12, br(), DTOutput(outputId = 'outcodontable', width = '100%'))), br(), br(),
      fluidRow(column(width =12, br(), plotOutput(outputId = 'aas', width = '100%'))), br(), br(),
      br(),br()
    )
  )
)