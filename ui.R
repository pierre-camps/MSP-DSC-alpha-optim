###########################
# MSP_DSC
# version 1.0 Fevrier 2018.
# version 2.0 October 2026.
# version 2.1 October 2026. DSC case 5
# version 2.2 October 2026. Tab 'Compare cases' (paired bootstrap)
# version 2.2.1 October 2026. Tab 'Case 5 vs case 3'
###########################

# Packages de base, installés avec R 
library(base)
library(datasets)
library(graphics)
library(grDevices)
library(methods)
library(stats)
library(utils)

# library('boot',lib='~/librairies')
# library('yaml',lib='~/librairies')
# library('withr',lib='~/librairies')
# library('dplyr',lib='~/librairies')
# library('pracma',lib='~/librairies')
# library('ggplot2',lib='~/librairies')
# library('shiny',lib='~/librairies')
# library('shinythemes',lib='~/librairies')
# library('htmlwidgets',lib='~/librairies')
# library('markdown',lib='~/librairies')
# library('DT',lib='~/librairies')

# 
# 
library(shiny)
library(shinythemes)
library(markdown)
library(DT)
library(withr)
library(dplyr)
library(ggplot2)
library(pracma)
library(boot)
library(htmlwidgets)
library(yaml)
library(ggsci)
library(ggrepel)


#fluidPage(theme = shinytheme("darkly"),
fluidPage(theme = shinytheme("slate"),
  # tables scroll inside a frame and keep their column headers visible
  tags$head(tags$style(HTML("
    /* coloured headings in the Readme tab */
    .readme h3 { color: #5bc0de; }
    .readme h4 { color: #f0ad4e; }
    .sticky-table { max-height: 70vh; overflow-y: auto; margin-bottom: 10px; }
    .sticky-table table thead th {
      position: sticky; top: 0; z-index: 2;
      background-color: #272b30;
      box-shadow: inset 0 -1px 0 #7a8288;
    }
  "))),
  navbarPage( "MultiSpecimen Paleointensity",
    id = "tabs",
    tabPanel(
      "Import Dataset",
      HTML('<img src="logoGM.png" width="180" height="69">'),
      headerPanel("MSP-DSC data PLOTTER and ANALYSER"),
      tags$hr(),
      br(),
      p('created by', span("Pierre Camps", style = "color:orangered"), ', CNRS and Univ. Montpellier.'), 
      p(
        span(style = "color:orangered"),
        'Please report bugs to :',
        span("pierre.camps [at] umontpellier.fr", style = "color:orangered")
      ),
      
      h4(
        span("To start, first UPLOAD a datafile", style = "color:orangered")
      ),
      
      fileInput(inputId= 'ifile', 
                          label  = 'Choose file to upload',
                          accept = c( 'text/csv',
                                      'text/comma-separated-values',
                                      '.csv')),
      br(),
      p('Please, see "Readme" Tab for informations on file format')
                
    ), #End tabPanel Open file
    tabPanel("Data explorer", dataTableOutput("contents")
    ),# End TabPanel Data Explorer
    tabPanel("Analysis and Plots",
      sidebarLayout(
        sidebarPanel(
          radioButtons(inputId = 'Protocol', 
                   label = 'Select the protocol type',
                   choices =list("DB","FC","DSC"),
                   selected = "DB",
                   inline = TRUE),
          conditionalPanel(
            condition = "input.Protocol == 'DSC'",
            selectInput(inputId = "caseDSC",
                        label = "MSP-DSC: alpha strategy",
                        choices = c("Case 1: constant alpha" = "1",
                                    "Case 2: single optimized alpha [0, 0.8]" = "2",
                                    "Case 3: per-specimen alpha, bounds [0, 0.8]" = "3",
                                    "Case 4: per-specimen alpha, bounds [a_min, a_max]" = "4",
                                    "Case 5: per-specimen alpha, prior centred on [a_min, a_max]" = "5"),
                        selected = "1"),
            conditionalPanel(
              condition = "input.caseDSC == '3' || input.caseDSC == '4'",
              p(strong("lambda (LOOCV): "), textOutput("lambdaInfo", inline = TRUE))
            ),
            conditionalPanel(
              condition = "input.caseDSC == '4'",
              helpText("a_min and a_max are read from the last two columns of the data file.")
            ),
            conditionalPanel(
              condition = "input.caseDSC == '5'",
              sliderInput(inputId = "kSigma",
                          label = "Width of [a_min, a_max] in standard deviations of the prior",
                          min = 2, max = 8, step = 1, value = 4),
              p(strong("lambda: "), textOutput("lambdaInfo5", inline = TRUE)),
              helpText("alpha is searched in [0, 0.8]. The prior of each specimen is centred on the middle of its interval [a_min, a_max] (last two columns of the data file), with sigma = (a_max - a_min) / k and lambda = 1 / sigma^2. No LOOCV.")
            )
          ),
          conditionalPanel(
            condition = "input.Protocol != 'DB'",
            sliderInput(inputId ="sliderF", 
                    label = "Fraction range", 
                    min = 0, 
                    max = 1, 
                    step = 0.1,
                    value = c(0.2, 0.8)),
            sliderInput(inputId ="Angle", 
                    label = "Threshhold angle between NRM and pTRM", 
                    min = 5, 
                    max = 35, 
                    step = 5,
                    value = 10)
          ),
          conditionalPanel(
            condition = "input.Protocol == 'DSC' && (input.caseDSC == '1' || input.caseDSC == '2')",
            p(strong("alpha: "), textOutput("alphaInfo", inline = TRUE))
          ),
          sliderInput(inputId = "NBoot",
                  label= "Number of bootstrapped values:",
                  min=5000,
                  max=15000,
                  step = 2500,
                  value = 10000)
        ), # End sidebarPAnel
        mainPanel(
           tabsetPanel(id = "resTabs",
           tabPanel("Q values", h4("Table of Q values"), div(class = "sticky-table", tableOutput("text1")), downloadButton('download','Save data as CSV file')),
           tabPanel("Plot MSP", plotOutput("QDB",width = "700px", height = "400px"),
                 downloadButton('downplot','Download the plot'),
                 radioButtons(inputId = 'formatF', 
                              label = '',
                              choices =list("pdf","png","tiff"),
                              selected = "png",
                              inline = TRUE)),
           tabPanel("Plot cook", plotOutput("Pcook",width = "700px", height = "400px"),
                    downloadButton('downplot2','Download the plot'),
                    radioButtons(inputId = 'formatF2', 
                                 label = '',
                                 choices =list("pdf","png","tiff"),
                                 selected = "png",
                                 inline = TRUE)
                    ),
           tabPanel("Plot boot", plotOutput("Pboot",width = "700px", height = "400px"),
                    sliderInput(inputId = "xrangeBoot",
                                label = "x-range of the bootstrap plot (\u00b5T)",
                                min = 0, max = 100, step = 5, value = c(0, 100),
                                width = "700px"),
                    downloadButton('downplot3','Download the plot'),
                    radioButtons(inputId = 'formatF3', 
                                 label = '',
                                 choices =list("pdf","png","tiff"),
                                 selected = "png",
                                 inline = TRUE)
                   ), # End TabPanel Plot Bbot
           tabPanel("Alpha values", h4("Optimized alpha per specimen"),
                    helpText("Bad specimens and Cook outliers are not optimized: their alpha is the ridge centre (0.5 in case 3, (a_min + a_max) / 2 in cases 4 and 5). In case 5, a_min and a_max define the prior; alpha itself is searched in [0, 0.8]."),
                    div(class = "sticky-table", tableOutput("alphaTable")),
                    downloadButton('downloadAlpha','Save alpha values as CSV file')),
           tabPanel("Case 4 vs case 3", plotOutput("Palpha", width = "600px", height = "600px"),
                    downloadButton('downplot4','Download the plot'),
                    radioButtons(inputId = 'formatF4', 
                                 label = '',
                                 choices =list("pdf","png","tiff"),
                                 selected = "png",
                                 inline = TRUE)),
           tabPanel("Case 5 vs case 3", plotOutput("Palpha5", width = "600px", height = "600px"),
                    downloadButton('downplot7','Download the plot'),
                    radioButtons(inputId = 'formatF7', 
                                 label = '',
                                 choices =list("pdf","png","tiff"),
                                 selected = "png",
                                 inline = TRUE)),
           tabPanel("Lambda", plotOutput("Plambda", width = "700px", height = "400px"),
                    downloadButton('downplot5','Download the plot'),
                    radioButtons(inputId = 'formatF5', 
                                 label = '',
                                 choices =list("pdf","png","tiff"),
                                 selected = "png",
                                 inline = TRUE)),
           tabPanel("Compare cases", h4("Comparison of the five MSP-DSC cases"),
                    helpText("The five cases are computed on the same specimens: their paleointensities are correlated and cannot be compared with a test for independent estimates. At each replicate, the same resampled specimens are used for the five cases (paired bootstrap), and the difference of the paleointensities of two cases is recorded. A difference is significant when its bootstrap interval does not include 0 (p-values corrected for the 10 comparisons by the Holm method). The comparison is run only when the five cases have been computed: select each case in 'MSP-DSC: alpha strategy', then click on the button. lambda (cases 3 and 4), the exclusion of the influential specimens and the weights are kept fixed, as in the bootstrap of each case. The number of replicates is the number of bootstrap samples of the left panel; the computation takes a few tens of seconds for 10000 replicates."),
                    verbatimTextOutput("compareStatus"),
                    actionButton("runCompare", "Run the comparison of the five cases"),
                    h5("Paleointensity of each case"),
                    tableOutput("compareEst"),
                    h5("Differences between cases (later case minus earlier case)"),
                    tableOutput("compareDiff"),
                    verbatimTextOutput("compareNote"),
                    plotOutput("Pcompare", width = "1000px", height = "650px"),
                    downloadButton('downplot6','Download the plot'),
                    radioButtons(inputId = 'formatF6',
                                 label = '',
                                 choices =list("pdf","png","tiff"),
                                 selected = "png",
                                 inline = TRUE),
                    br(),
                    downloadButton('downloadCompare','Save the differences as CSV file')),
           tabPanel("Summary", verbatimTextOutput("summary"))
           )# end tabsetPanel
        )# end mainPanel
      )# end sidebarLayout
    ),# end tabPanel ANALYSIS
    tabPanel(
      "Readme",
      HTML('<img src="logoGM.png" width="180" height="69">'),
      br(),
      h4(span("Cautionary notes :", class="text-warning")),
      HTML('<div class="card border-danger mb-3" style="max-width: 120rem;">
              <span class="badge badge-primary">Note 1</span>
              <div class="card-body">
                <p class="card-text">This software is provided on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or implied.</p>
              </div>
              <hr class="my-4">
              <span class="badge badge-info">Note 2</span>
              <div class="card-body">
                <p class="card-text">Papers reporting results obtained with the present source code may cite
                in references or bibliography the paper by Camps et al., SGG, XX, (X), pp-pp, 2026.</p>
              </div>
              <hr class="my-4">
              <span class="badge badge-primary">Note 3</span>
              <div class="card-body">
                <p class="card-text">Please report bugs to pierre.camps [at] umontpellier.fr.</p>
              </div>
              <hr class="my-4">
              <br> 
              <p class="text-success">CODE created by Pierre Camps CNRS and Univ. Montpellier.</p>
              </div>'),
             
             downloadButton("downloadExample", label = "Get example of inputfile"),div(class = "readme", includeMarkdown("ReadMe.md"))
    ) # End TabPanel Readme
  )# END Navbarpage
)# END fluidPage
