# https://phaden.shinyapps.io/seght_shiny/

#========================================================================================
# ui
#========================================================================================
ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  shinyFeedback::useShinyFeedback(),
  
  
  titlePanel("Segmented Hypothesis Testing Toolbox"),
  

  
  fluidRow(
    
    column(3, 
          selectInput("txt_max_n_segments",
                       h4("Max N Segments"), 
                       c( "1" =  "1",  "2" =  "2",  "3" =  "3",  "4" =  "4",  "5" =  "5", 
                          "6" =  "6",  "7" =  "7",  "8" =  "8",  "9" =  "9", "10" = "10", 
                         "11" = "11", "12" = "12", "13" = "13", "14" = "14", "15" = "15", 
                         "16" = "16", "17" = "17", "18" = "18", "19" = "19", "20" = "20"),
                       selected = "3")),
    column(2, 
           textInput("txt_alpha_total", h4("Alpha Total"), 
                     value = "0.05")),
    column(2, 
           textInput("txt_alpha_strong", h4("Alpha Strong"), 
                     value = "0.025"))
    
  ),
  
  fluidRow(           
    column(3, actionButton("btn_alpha_weak", "Compute Alpha Weak", width = "200px")),
    column(2, h4("Alpha Weak:")),
    column(2, h4(textOutput("txt_out_alpha_weak"))),
    column(4, span(textOutput("txt_out_error_message"), style="color:red"))
  ),
  
  hr(),
  
  
  #h4("Enter these values to compute power:"),
  
  fluidRow(
    
    column(3, 
           selectInput("lst_stat_procedure_name", h4("Statistical Test"), 
                       c("One Sample t" = "1t", 
                         "Two Sample t" = "2t", 
                         "One Sample z" = "1z", 
                         "Two Sample z" = "2z", 
                         "Pearson r"    = "r"))),
    column(2, 
           textInput("txt_effect_size", h4("Effect Size"), 
                     value = "0.5")),
    column(2, 
           textInput("txt_n_per_segment", h4("N per Segment"), 
                     value = "30"))
  ),
  
  fluidRow(           
    column(3, actionButton("btn_power", "Compute Power", width = "200px")),
    column(2, h4("Power:")),
    column(2, h4(textOutput("txt_out_power")))
  ),
  
  hr(),
  
  
  
  fluidRow(
    column(3, textInput("txt_base_rate", h4("Base Rate"), value = "0.8"))
  ),
  
  fluidRow(
    column(3, actionButton("btn_total_n", "Compute Total N Values", width = "200px")),
    column(2, h4("Expected Total N:")),
    column(2, h4(textOutput("txt_out_total_n"))),
    column(2, h4("Std. Dev. of Total N:")),
    column(2, h4(textOutput("txt_out_std_total_n")))
  ),
  
  fluidRow(column(12, h4("   "))),
  
  fluidRow(
    column(3),
    column(9,
           tableOutput('tbl_outcomes'))
  )
  
) # end ui <- fluidpage(...)


#========================================================================================
# server <- function(...)
#========================================================================================
# Define server logic required to update with computations ----
server <- function(input, output, session) {
  
  #===============================================================
  # Button click event handlers (use observeEvent to allow multiple
  # input dependencies)
  #===============================================================
  
  # Handler for btn_alpha_weak
  
  observeEvent(input$btn_alpha_weak, {
    alpha_total <- as.numeric(input$txt_alpha_total)
    max_n_segments <- as.numeric(input$txt_max_n_segments)
    alpha_strong <- as.numeric(input$txt_alpha_strong)
    alpha_weak <- alpha_weak(max_n_segments, alpha_total, alpha_strong)
    output$txt_out_alpha_weak <- renderText(round(alpha_weak,3))
  })
  
  
  # Handler for btn_power
  observeEvent(input$btn_power, {
    alpha_strong <- as.numeric(input$txt_alpha_strong)
    alpha_total <- as.numeric(input$txt_alpha_total)
    stat_procedure_name <- input$lst_stat_procedure_name
    effect_size <- as.numeric(input$txt_effect_size)
    n_per_segment <- as.numeric(input$txt_n_per_segment)
    max_n_segments <- as.numeric(input$txt_max_n_segments)
    
    outcomes <- segmented_hyp_test_outcomes(max_n_segments,
                                            n_per_segment,
                                            alpha_total,
                                            alpha_strong,
                                            stat_procedure_name, 
                                            effect_size,
                                            base_rates = 1
    )
    
    power <- outcomes$avg_power
    output$txt_out_power <- renderText(round(power,3))
    
    # Recompute alpha weak as well (spec)
    alpha_weak <- alpha_weak(max_n_segments, alpha_total, alpha_strong)
    output$txt_out_alpha_weak <- renderText(round(alpha_weak,3))
  })
  
    
  # Handler for btn_total_n
  observeEvent(input$btn_total_n, {
    alpha_strong <- as.numeric(input$txt_alpha_strong)
    max_n_segments <- as.numeric(input$txt_max_n_segments)
    n_per_segment <- as.numeric(input$txt_n_per_segment)
    alpha_total <- as.numeric(input$txt_alpha_total)
    stat_procedure_name <- input$lst_stat_procedure_name
    effect_size <- as.numeric(input$txt_effect_size)
    base_rate <- as.numeric(input$txt_base_rate)
    
    outcomes <- segmented_hyp_test_outcomes(max_n_segments,
                                            n_per_segment,
                                            alpha_total,
                                            alpha_strong,
                                            stat_procedure_name,
                                            effect_size,
                                            base_rate)
    expected_total_n <- outcomes$exp_n_subj
    sd_total_n <- outcomes$sd_exp_n_subj
    
    # Recompute alpha weak as well (spec)
    alpha_weak <- alpha_weak(max_n_segments, alpha_total, alpha_strong)
    output$txt_out_alpha_weak <- renderText(round(alpha_weak,3))
    
    # Update power (spec)
    power <- outcomes$avg_power
    output$txt_out_power <- renderText(round(power,3))

    output$txt_out_total_n <- renderText(round(expected_total_n,2))
    output$txt_out_std_total_n <- renderText(round(sd_total_n,2))
    
    # Outcome table
    pr_reject <- outcomes$pr_reject_by_segment
    pr_ftr <- outcomes$pr_ftr_by_segment
    segment <- 1:length(pr_reject)
    
    show("tbl_outcomes")
    outcomes_df <- data.frame(segment, pr_reject, pr_ftr)
    colnames(outcomes_df) <- c("Segment", "Pr(reject)", "Pr(fail to reject)")
    output$tbl_outcomes <- renderTable(outcomes_df,
                                       striped = TRUE,
                                       bordered = TRUE,
                                       spacing = "s",
                                       align = "c")
    })
      
    
  
  #===============================================================
  # Code to respond to changes in input boxes, including:
  #  o check validity of values typed into input text boxes
  #  o clear output text boxes whose values are no longer valid
  #===============================================================
  
   handle1Segment <- function(){
    # Special handling when max_n_segments = 1:
    #  set alpha_weak = alpha_strong = alpha_total
    v <- input$txt_max_n_segments
    vnum = as.numeric(v)
    if (vnum==1) {
      output$txt_out_alpha_weak <- renderText(input$txt_alpha_total)
      updateTextInput(session,"txt_alpha_strong", value = input$txt_alpha_total)
      return()
    } else {
      if (input$txt_alpha_total == input$txt_alpha_strong) {
        # reduce alph_strong when you change from 1 to many segments
        newAlphaStrong = as.numeric(input$txt_alpha_total) / 2  # a default
        updateTextInput(session,"txt_alpha_strong", value = paste(newAlphaStrong))
      }
    }
  }

  clearOutcomes <- function(){
    output$txt_out_total_n <- renderText("")
    output$txt_out_std_total_n <- renderText("")
    hide("tbl_outcomes")
  }
  
 #  Change in max_n_segments:
  observeEvent(input$txt_max_n_segments,
               {
                 output$txt_out_alpha_weak <- renderText("")
                 handle1Segment()
                 output$txt_out_power <- renderText("")
                 clearOutcomes()
               })

  #  Change in alpha_total:
  observeEvent(input$txt_alpha_total,
               {
                 output$txt_out_alpha_weak <- renderText("")
                 handle1Segment()
                 output$txt_out_power <- renderText("")
                 clearOutcomes()
               })
  
  #  Change in alpha_strong:
  observeEvent(input$txt_alpha_strong,
               {
                 output$txt_out_alpha_weak <- renderText("")
                 handle1Segment()
                 output$txt_out_power <- renderText("")
                 clearOutcomes()
               })
  
  # Change that affect power: i.e., in statistical procedure, sample size, or effect size.
  inputs_power <- reactive({
    c(input$txt_effect_size, input$txt_n_per_segment, input$lst_stat_procedure_name)
  })

  observeEvent(inputs_power(),
               {
                 output$txt_out_power <- renderText("")
                 clearOutcomes()
               }
              )
  
  observeEvent(input$txt_base_rate,  
               {
                 clearOutcomes()
               }
              )
  
  #===============================================================
  # Check legal input values
  #===============================================================
  
  # txt_alpha_total
  observeEvent(input$txt_alpha_total,
               {
                 req(input$txt_alpha_total)

                 if (input$txt_alpha_total < '0' | input$txt_alpha_total > '1') {
                   showFeedbackWarning(inputId = "txt_alpha_total",
                                       text =  "Alpha Total must be between 0 and 1")
                 } else {
                   hideFeedback(inputId = "txt_alpha_total")
                 }
               })
  
  # txt_alpha_strong
  observeEvent(input$txt_alpha_strong,
               {
                 req(input$txt_alpha_strong)
                 
                 if (input$txt_alpha_strong < '0' | input$txt_alpha_strong > '1') {
                   showFeedbackWarning(inputId = "txt_alpha_strong",
                                       text =  "Alpha Strong must be between 0 and 1")
                 } else {
                   hideFeedback(inputId = "txt_alpha_strong")
                 }
               })
  
  # txt_effect_size
  observeEvent(input$txt_effect_size,
               {
                 req(input$txt_effect_size)
                 
                 if (input$txt_effect_size < '0' | input$txt_effect_size > '1') {
                   showFeedbackWarning(inputId = "txt_effect_size",
                                       text =  "Effect Size must be between 0 and 1")
                 } else {
                   hideFeedback(inputId = "txt_effect_size")
                 }
               })
  
  # txt_n_per_segment
  observeEvent(input$txt_n_per_segment,
               {
                 req(input$txt_n_per_segment)
                 
                 if (input$txt_n_per_segment < '0' | input$txt_n_per_segment > 'a') {
                   showFeedbackWarning(inputId = "txt_n_per_segment",
                                       text =  "N per segment must be greater than 0")
                 } else {
                   hideFeedback(inputId = "txt_n_per_segment")
                 }
               })
  
  # txt_base_rate
  observeEvent(input$txt_base_rate,
               {
                 req(input$txt_base_rate)
                 
                 if (input$txt_base_rate < '0' | input$txt_base_rate > '1') {
                   showFeedbackWarning(inputId = "txt_base_rate",
                                       text =  "Base Rate must be between 0 and 1")
                 } else {
                   hideFeedback(inputId = "txt_base_rate")
                 }
               })
  
  
  # txt_alpha_total > txt_alpha_strong
  observeEvent(c(input$txt_alpha_total, input$txt_alpha_strong),
               {
                 error_message <- ""
                 
                 req(input$txt_alpha_total)
                 req(input$txt_alpha_strong)
                
                 a_t <- as.numeric(input$txt_alpha_total)
                 a_s <- as.numeric(input$txt_alpha_strong)
                   
                 if (!is.na(a_t) & !is.na(a_s)) {
                   if (a_s > a_t) {
                     error_message <- "Alpha strong must be less than alpha total."
                   }
                 }
               
                 
                 output$txt_out_error_message <- renderText(error_message)
               })
                   
  
} # end server


#========================================================================================
# Launch
#========================================================================================


options(warn = -1)

# Install packages if necessary
if (!require(V8)) install.packages('V8')
if (!require(pracma)) install.packages('pracma')
if (!require(stats)) install.packages('stats')
if (!require(SuppDists)) install.packages('SuppDists')
if (!require(shinyjs)) install.packages('shinyjs')
if (!require(shinyFeedback)) install.packages('shinyFeedback')

library(V8)
library(pracma)
library(SuppDists)
library(stats)
library(shinyjs)
library(shinyFeedback)

install.packages("segHT_v1.4.0.tar.gz", repos = NULL, type="source")
library(segHT)

shinyApp(ui = ui, server = server)
