library(shiny)
rm(list = ls())
source("ReadDataApp.R")
source("MCMC_HierarchicalThompson.R")

ui <- fluidPage(
#  titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
  verticalLayout(
    includeMarkdown("instructions.md"),
    hr(),
    fluidRow(column(4, fileInput("file1", "Choose CSV file of previous data",
                              multiple = FALSE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv"))),
             column(4,numericInput(inputId="RR", label="Replicate draws", value=5000, min=1)),
             column(4,numericInput(inputId="alpha", label="Share fully randomized", value=.2, min=0, max=1))
             ),
    fluidRow(column(4,actionButton(inputId = "calcbutton", label = "Calculate treatment assignment")),
             column(4,downloadButton("downloadData", "Download assignment probabilities"))),
    hr(),
    fluidRow(column(3,tableOutput("trials")),
             column(3,tableOutput("successrates")),
             column(3,tableOutput("thompsonprobabilities")),
             column(3,tableOutput("actualprobabilities"))
             ) 
  )
)


server <- function(input, output, session) {
  
  v = reactiveValues()
  
  observeEvent(input$calcbutton,{
    req(input$file1)
    
    #loading file
    priordata=ReadDataApp(input$file1$datapath)
    
    #calculating treatment assignment
    if (priordata$nx > 1) {
      v$Pstar=DtchoiceMCMCProbabilities(priordata$Y,priordata$D,priordata$X, #outcomes, treatments, and covariates thus far
                                priordata$k,priordata$nx, #number of treatments and number of strata
                                RR=input$RR)
      
      # summary statistics of data to display
      v$SS=tapply(priordata$Y,list(priordata$D,priordata$X),sum, default=0) #matrix of successes
      v$NN=tapply(priordata$Y,list(priordata$D,priordata$X),length, default=0) #matrix of trials
    } else {
      v$Pstar=DtchoiceThompsonProbabilities(priordata$Y,priordata$D, #outcomes and treatments thus far
                                     priordata$k, #number of treatments
                                     RR=input$RR)
      
      # summary statistics of data to display
      v$SS=tapply(priordata$Y,list(priordata$D),sum, default=0) #matrix of successes
      v$NN=tapply(priordata$Y,list(priordata$D),length, default=0) #matrix of trials
    }
    v$SR=v$SS/v$NN
    v$trials=as_tibble(t(v$NN))
    v$successrates=as_tibble(t(v$SR))
    
    v$Pactual=(1-input$alpha) *v$Pstar + input$alpha*(1/priordata$k)

  })
  
  
  output$trials =  renderTable(
    v$trials,
    align="c",
    digits=0,
    caption="Observations",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
  output$successrates =  renderTable(
    v$successrates,
    align="c",
    caption="Prior success rates",
    caption.placement = getOption("xtable.caption.placement", "top")
  ) 
  output$thompsonprobabilities =  renderTable(
    v$Pstar, 
    align="c",
    caption="Thompson probabilities",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
  output$actualprobabilities =  renderTable(
    v$Pactual, 
    align="c",
    caption="Actual probabilities",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
  
#download optimal design
  output$downloadData <- downloadHandler(
    filename = paste(Sys.Date(), "treatmentprobabilities.csv", sep=""),
    content = function(file) {
      write_csv(v$Pactual, file)
    }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)