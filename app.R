library(shiny)
library(readxl)
library(ggplot2)
library(dplyr)

# Source your functions
source("scripts/dsa.R")
source("scripts/tornado.R")

# Load data
m_TP <- read_excel("data/inputs/transition_probabilities.xlsx") |> 
  as.data.frame() |>
  tibble::column_to_rownames(var = "HS") |>
  as.matrix()

# Load parameter descriptions for labeling
param_descriptions <- read_excel("data/inputs/dsa_parameters.xlsx") |>
  select(parameter, description)

# Fixed parameters
v_init_state <- c(H = 1, S1 = 0, S2 = 0, R = 0, D = 0)
n_cycles <- 30

# UI
ui <- fluidPage(
  titlePanel("\"Sick-Sicker\" Markov Model: Deterministic Sensitivity Analysis"),
  
  # Description section
  fluidRow(
    column(12,
           wellPanel(
             h3("About This Model"),
             p("This interactive tool demonstrates a \"Sick-Sicker\" Markov cohort model", 
               " for cost-effectiveness analysis of two competing treatment strategies,",
               " Standard Care and Advanced Therapy, for a generic disease."),
             
             p("The model simulates a cohort of patients transitioning through five health states over 30 annual cycles:"),
             tags$ul(
               tags$li(strong("Healthy (H):"), " Disease-free state"),
               tags$li(strong("Sick Stage 1 (S1):"), " Early disease with moderate symptoms"),
               tags$li(strong("Sick Stage 2 (S2):"), " Advanced disease with severe symptoms"),
               tags$li(strong("Recovery (R):"), " Patients who recover from illness"),
               tags$li(strong("Death (D):"), " Absorbing state")
             ),
             
             h4("Treatment Strategies"),
             tags$ul(
               tags$li(strong("Baseline:"), " No treatment"),
               tags$li(strong("Standard Care:"), " Reduces disease progression to all worse health states by 20%"),
               tags$li(strong("Advanced Therapy:"), " Reduces disease progression to all worse health states by 50%")
             ),
             
             h4("How to Use This Tool"),
             p("Use the tabs below to adjust DSA parameters and view results. Click 'Run Analysis' to calculate the ICER and generate the tornado diagram."),
             p("You can also adjust the Willingness-to-Pay (WTP) threshold using the slider in the Parameters tab to see how cost-effectiveness changes relative to different threshold values.")
           )
    )
  ),
  
  # Tabbed interface
  tabsetPanel(id = "tabs",
              # Parameters tab
              tabPanel("Parameters",
                       value = "params_tab",
                       br(),
                       fluidRow(
                         column(6,
                                wellPanel(
                                  h4("Quality of Life (Utilities)"),
                                  h5("Utility: Healthy"),
                                  numericInput("u_H_base", "Base:", value = 1, min = 0, max = 1, step = 0.01),
                                  sliderInput("u_H_range", "Range:", min = 0.9, max = 1, value = c(0.95, 1), step = 0.01),
                                  
                                  h5("Utility: Sick Stage 1"),
                                  numericInput("u_S1_base", "Base:", value = 0.75, min = 0, max = 1, step = 0.01),
                                  sliderInput("u_S1_range", "Range:", min = 0.5, max = 1, value = c(0.65, 0.85), step = 0.01),
                                  
                                  h5("Utility: Sick Stage 2"),
                                  numericInput("u_S2_base", "Base:", value = 0.5, min = 0, max = 1, step = 0.01),
                                  sliderInput("u_S2_range", "Range:", min = 0.2, max = 0.8, value = c(0.4, 0.6), step = 0.01),
                                  
                                  h5("Utility: Recovery"),
                                  numericInput("u_R_base", "Base:", value = 0.9, min = 0, max = 1, step = 0.01),
                                  sliderInput("u_R_range", "Range:", min = 0.7, max = 1, value = c(0.8, 1), step = 0.01)
                                ),
                                
                                wellPanel(
                                  h4("Annual State Costs (£)"),
                                  h5("Cost: Sick Stage 1"),
                                  numericInput("c_S1_base", "Base:", value = 4000, min = 0, step = 100),
                                  sliderInput("c_S1_range", "Range:", min = 2000, max = 6000, value = c(3000, 5000), step = 100),
                                  
                                  h5("Cost: Sick Stage 2"),
                                  numericInput("c_S2_base", "Base:", value = 15000, min = 0, step = 500),
                                  sliderInput("c_S2_range", "Range:", min = 10000, max = 20000, value = c(12000, 18000), step = 500)
                                ),
                                
                                wellPanel(
                                  h4("Treatment Costs (£)"),
                                  h5("Cost: Standard Care"),
                                  numericInput("c_trt_standard_base", "Base:", value = 3000, min = 0, step = 100),
                                  sliderInput("c_trt_standard_range", "Range:", min = 2000, max = 5000, value = c(2400, 3600), step = 100),
                                  
                                  h5("Cost: Advanced Therapy"),
                                  numericInput("c_trt_advanced_base", "Base:", value = 8000, min = 0, step = 200),
                                  sliderInput("c_trt_advanced_range", "Range:", min = 5000, max = 12000, value = c(6400, 9600), step = 200)
                                )
                         ),
                         
                         column(6,
                                wellPanel(
                                  h4("Risk Ratios: Standard Care"),
                                  h5("RR: S1 → S2"),
                                  numericInput("rr_S1_S2_standard_base", "Base:", value = 0.8, min = 0, max = 1, step = 0.05),
                                  sliderInput("rr_S1_S2_standard_range", "Range:", min = 0.5, max = 1, value = c(0.65, 0.95), step = 0.05),
                                  
                                  h5("RR: S1 → Death"),
                                  numericInput("rr_S1_D_standard_base", "Base:", value = 0.8, min = 0, max = 1, step = 0.05),
                                  sliderInput("rr_S1_D_standard_range", "Range:", min = 0.5, max = 1, value = c(0.65, 0.95), step = 0.05),
                                  
                                  h5("RR: S2 → Death"),
                                  numericInput("rr_S2_D_standard_base", "Base:", value = 0.8, min = 0, max = 1, step = 0.05),
                                  sliderInput("rr_S2_D_standard_range", "Range:", min = 0.5, max = 1, value = c(0.65, 0.95), step = 0.05)
                                ),
                                
                                wellPanel(
                                  h4("Risk Ratios: Advanced Therapy"),
                                  h5("RR: S1 → S2"),
                                  numericInput("rr_S1_S2_advanced_base", "Base:", value = 0.5, min = 0, max = 1, step = 0.05),
                                  sliderInput("rr_S1_S2_advanced_range", "Range:", min = 0.2, max = 0.8, value = c(0.35, 0.7), step = 0.05),
                                  
                                  h5("RR: S1 → Death"),
                                  numericInput("rr_S1_D_advanced_base", "Base:", value = 0.5, min = 0, max = 1, step = 0.05),
                                  sliderInput("rr_S1_D_advanced_range", "Range:", min = 0.2, max = 0.8, value = c(0.35, 0.7), step = 0.05),
                                  
                                  h5("RR: S2 → Death"),
                                  numericInput("rr_S2_D_advanced_base", "Base:", value = 0.5, min = 0, max = 1, step = 0.05),
                                  sliderInput("rr_S2_D_advanced_range", "Range:", min = 0.2, max = 0.8, value = c(0.35, 0.7), step = 0.05)
                                ),
                                
                                wellPanel(
                                  h4("Discount Rate"),
                                  numericInput("discount_rate_base", "Base:", value = 0.03, min = 0, max = 0.1, step = 0.01),
                                  sliderInput("discount_rate_range", "Range:", min = 0, max = 0.1, value = c(0, 0.05), step = 0.01)
                                ),
                                
                                wellPanel(
                                  h4("Willingness-to-Pay Threshold"),
                                  sliderInput("wtp_threshold", "WTP Threshold (£/QALY):", 
                                              min = 0, max = 100000, value = 20000, step = 5000),
                                )
                         )
                       ),
                       
                       fluidRow(
                         column(12, align = "center",
                                br(),
                                actionButton("run", "Run Analysis", class = "btn-primary btn-lg", 
                                             style = "width: 300px; font-size: 18px;")
                         )
                       )
              ),
              
              # Results tab
              tabPanel("Results",
                       value = "results_tab",
                       br(),
                       fluidRow(
                         column(12,
                                wellPanel(
                                  h3("Base Case Cost-Effectiveness"),
                                  verbatimTextOutput("base_results", placeholder = TRUE)
                                )
                         )
                       ),
                       
                       fluidRow(
                         column(12,
                                h3("Tornado Diagram: One-Way Sensitivity Analysis"),
                                p("This diagram shows which parameters have the greatest impact on the ICER when varied individually between their lower and upper bounds."),
                                p(strong("Willingness-to-Pay Threshold: "), 
                                  textOutput("wtp_text", inline = TRUE),
                                  style = "font-size: 16px"),
                                plotOutput("tornado", height = "700px")
                         )
                       )
              )
  )
)

# Server
server <- function(input, output, session) {
  
  # Switch to Results tab when Run Analysis is clicked
  observeEvent(input$run, {
    updateTabsetPanel(session, "tabs", selected = "results_tab")
  })
  
  results <- eventReactive(input$run, {
    
    # Build base_params from base value inputs
    base_params <- list(
      v_utilities = c(
        H = input$u_H_base, 
        S1 = input$u_S1_base, 
        S2 = input$u_S2_base, 
        R = input$u_R_base, 
        D = 0
      ),
      v_costs = c(
        H = 0, 
        S1 = input$c_S1_base, 
        S2 = input$c_S2_base, 
        R = 0, 
        D = 0
      ),
      trt_standard = list(
        name = "Standard Care",
        transitions = list(
          "S1_S2" = input$rr_S1_S2_standard_base,
          "S1_D" = input$rr_S1_D_standard_base,
          "S2_D" = input$rr_S2_D_standard_base
        ),
        cost = input$c_trt_standard_base,
        apply_states = c("S1", "S2")
      ),
      trt_advanced = list(
        name = "Advanced Therapy",
        transitions = list(
          "S1_S2" = input$rr_S1_S2_advanced_base,
          "S1_D" = input$rr_S1_D_advanced_base,
          "S2_D" = input$rr_S2_D_advanced_base
        ),
        cost = input$c_trt_advanced_base,
        apply_states = c("S1", "S2")
      ),
      r = input$discount_rate_base
    )
    
    # Build dsa_params from range inputs
    dsa_params <- data.frame(
      parameter = c(
        "u_H", "u_S1", "u_S2", "u_R",
        "c_S1", "c_S2", "c_trt_standard", "c_trt_advanced",
        "rr_S1_S2_standard", "rr_S1_D_standard", "rr_S2_D_standard",
        "rr_S1_S2_advanced", "rr_S1_D_advanced", "rr_S2_D_advanced",
        "discount_rate"
      ),
      base_value = c(
        input$u_H_base, input$u_S1_base, input$u_S2_base, input$u_R_base,
        input$c_S1_base, input$c_S2_base, input$c_trt_standard_base, input$c_trt_advanced_base,
        input$rr_S1_S2_standard_base, input$rr_S1_D_standard_base, input$rr_S2_D_standard_base,
        input$rr_S1_S2_advanced_base, input$rr_S1_D_advanced_base, input$rr_S2_D_advanced_base,
        input$discount_rate_base
      ),
      lower = c(
        input$u_H_range[1], input$u_S1_range[1], input$u_S2_range[1], input$u_R_range[1],
        input$c_S1_range[1], input$c_S2_range[1], input$c_trt_standard_range[1], input$c_trt_advanced_range[1],
        input$rr_S1_S2_standard_range[1], input$rr_S1_D_standard_range[1], input$rr_S2_D_standard_range[1],
        input$rr_S1_S2_advanced_range[1], input$rr_S1_D_advanced_range[1], input$rr_S2_D_advanced_range[1],
        input$discount_rate_range[1]
      ),
      upper = c(
        input$u_H_range[2], input$u_S1_range[2], input$u_S2_range[2], input$u_R_range[2],
        input$c_S1_range[2], input$c_S2_range[2], input$c_trt_standard_range[2], input$c_trt_advanced_range[2],
        input$rr_S1_S2_standard_range[2], input$rr_S1_D_standard_range[2], input$rr_S2_D_standard_range[2],
        input$rr_S1_S2_advanced_range[2], input$rr_S1_D_advanced_range[2], input$rr_S2_D_advanced_range[2],
        input$discount_rate_range[2]
      )
    )
    
    # Run DSA
    dsa_results <- run_owsa(
      m_TP, v_init_state, n_cycles, base_params, dsa_params,
      strategy1 = "Standard Care",
      strategy2 = "Advanced Therapy"
    )
    
    # Calculate base ICER
    base_icer <- run_comparison(m_TP, v_init_state, n_cycles, base_params,
                                "Standard Care", "Advanced Therapy")$icer
    
    list(dsa_results = dsa_results, base_icer = base_icer)
  })
  
  output$base_results <- renderText({
    res <- results()
    paste0(
      "Base Case ICER: £", format(round(res$base_icer), big.mark = ","), " per QALY gained\n\n",
      "Interpretation: Advanced Therapy costs £", format(round(res$base_icer), big.mark = ","), 
      " for each additional quality-adjusted life year compared to Standard Care."
    )
  })
  
  output$wtp_text <- renderText({
    paste0("£", format(input$wtp_threshold, big.mark = ","), " per QALY (shown as dotted red line)")
  })
  
  output$tornado <- renderPlot({
    res <- results()
    plot_tornado(res$dsa_results, res$base_icer, wtp = input$wtp_threshold, 
                 param_labels = param_descriptions)
  })
}

shinyApp(ui, server)