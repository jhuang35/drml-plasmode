fluidPage(
  titlePanel("Realistic Evaluations of Finite sample INference using Efficient Estimators (REFINE2)"),
  numericInput("obs", "Number of Simulation Samples:", 100, min = 2, max = 500),
  textInput("path", "Path of the data","./plas_data.csv"),
  sidebarLayout(
    sidebarPanel(
      h4("Simulation Parameters"),
      numericInput("random_seed", "Random Seed:", 1111, min = 1, max = 999999),
      numericInput("num_cf", "Number of Cross-fits:", 5, min = 2, max = 20),
      
      h4("Models for simulation"),
      textInput(inputId = "expForm",
                label = "PS Model", 
                value = "A ~ 1 + VAR_1 + VAR_2"), 
      textInput(inputId = "outForm",
                label = "Outcome Model", 
                value = "Y ~ 1 + A + VAR_1 + VAR_2"),
      
      h4("Models for estimation"),
      textInput(inputId = "expForm.est",
                label = "PS Model", 
                value = "A ~ 1 + VAR_1 + VAR_2"), 
      textInput(inputId = "outForm.est",
                label = "Outcome Model", 
                value = "Y ~ 1 + A + VAR_1 + VAR_2"), 
      
      h4("Estimation Method"),
      selectInput('mtd', 'Estimation Method', c("IPW","GComp","TMLE","AIPW","DCTMLE","DCAIPW"),
                  selected = "IPW"),
      selectInput('is.par', 'Library (not for IPW and GComp)', c("Smooth","Non-Smooth"),
                  selected = "Smooth"),
      actionButton("run", "Run")
    ),
    
    mainPanel(
      textOutput("exp.Form"),textOutput("out.Form"),textOutput("exp.Form.est"),textOutput("out.Form.est"),
      htmlOutput("res.text.1"),htmlOutput("res.text.2"),htmlOutput("res.text.3"),
      plotOutput("plot"), tableOutput('table')
    )
  ) #,
  # fluidRow(
  #   column(12,
  #          tableOutput('table')
  #   )
  # )
)
# shiny::runApp("~/Desktop/HuangGroup/drml-plasmode/REFINE2")
# rsconnect::deployApp()
# rsconnect::configureApp("shinyapp", size="xxlarge")