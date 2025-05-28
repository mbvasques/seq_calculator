#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)


# Define UI for application that draws a histogram
ui = page_sidebar(

    # Application title
  title = "Sequencing Calculator",

    # Sidebar 
  sidebar = sidebar(
    width = 400,
    bg = "white",
    textInput("project", "Project:"),
    accordion(
      accordion_panel(
        "Samples",
          checkboxInput("visium", "Visium"),
          numericInput("nsamples",
                       "Number of Samples:",
                       value = 1,
                       min = 1,
                       width = "167px"),
          conditionalPanel(
            condition = ("input.nsamples >= 1"),
            conditionalPanel(
              condition = ("input.visium == false"),
              uiOutput("nsamples1"),
              sliderInput("loss",
                          "Loss (%):",
                          min = 0,
                          max = 100,
                          step = 5,
                          value = 50)),
            conditionalPanel(
              condition = "input.visium == true",
              uiOutput("nsamples2")))),
      accordion_panel(
        "Sequencing",
          selectInput("sequencer",
                      "Sequencer:",
                      choices = c("NovaSeq 6000", "NovaSeq X Plus (10B)",
                                  "NovaSeq X Plus (25B)"),
                      width = "250px"),
          conditionalPanel(
            condition = "input.visium == true",
            numericInput("depth",
                         "Desired Depth (reads/spot):",
                         value = 50000,
                         min = 1,
                         width = "250px")),
          conditionalPanel(
            condition = "input.visium == false",
            numericInput("depth",
                         "Desired Depth (reads/cell):",
                         value = 50000,
                         min = 1,
                         width = "250px")),
          numericInput("rl",
                       "Read Length:",
                       value = 150,
                       min = 1,
                       width = "100px"),
          sliderInput("phix",
                      "PhiX input (%):",
                      min = 0,
                      max = 10,
                      step = 1,
                      value = 1)),
      br(),
      div(style = "margin-left: 200px", 
          actionButton("run", label = "Apply", width = "150px")),
      br()
    )
  ),

    # Main Panel: Show results
    card(card_header("Total Cells"), 
         height = "260px", fill = F,
         nav_panel("Total Cells", tableOutput("outTable")),
         div(style = "position:absolute; left: 20px; top: 180px; width: 200px",
             downloadButton("downloadDataTotal", "Download"))),
    card(card_header("Per Sample"),
         full_screen = T,
         nav_panel("Samples", DT::DTOutput("sampleOutTable")),
         div(style = "left: 20px; width: 200px",
             downloadButton("downloadData", "Download")))
)

# read seq_info.csv
seq_table = read.csv("seq_info.csv", sep = ",", header = T, row.names = 1)

# Define server logic
server = function(input, output) {
  
  
  # add input for n samples 
  n_samples1 = reactive({
    (req(input$nsamples)>0)
    lapply(1:req(input$nsamples), function(i) {
      fluidRow(
        column(width = 1,
               checkboxInput(inputId = paste0("add_",i),
                             label = NULL, 
                             value = ifelse(is.null(input[[paste0("add_",i)]]),
                                            T, input[[paste0("add_",i)]]))),
        column(width = 6,
               textInput(inputId = paste0("sample_name",i),
                         label = paste0("Sample ", i), 
                         value = ifelse(is.null(input[[paste0("sample_name", i)]]),
                                        "", input[[paste0("sample_name", i)]]))),
        column(width = 5,
               numericInput(inputId = paste0("sample_",i),
                            label = paste0("Cells Loaded:"), 
                            value = ifelse(is.null(input[[paste0("sample_", i)]]),
                                           25000, input[[paste0("sample_", i)]]), 
                            min = 1)))
      
    }) 
  })       
  
  output$nsamples1 = renderUI({n_samples1()}) 
  
  # add input for n samples (Visium) 
  n_samples2 = reactive({
    (req(input$nsamples)>0)
    if (input$visium) {
      lapply(1:req(input$nsamples), function(i) {
        fluidRow(
          column(width = 1,
                 checkboxInput(inputId = paste0("vadd_",i),
                               label = NULL, 
                               value = ifelse(is.null(input[[paste0("vadd_",i)]]),
                                              T, input[[paste0("vadd_",i)]]))),
          column(width = 5,
                 textInput(inputId = paste0("vsample_name",i),
                           label = paste0("Sample ", i), 
                           value = ifelse(is.null(input[[paste0("vsample_name", i)]]),
                                          "", input[[paste0("vsample_name", i)]]))),
          column(width = 6,
                 numericInput(inputId = paste0("vsample_",i),
                              label = paste0("Est. Covered Spots:"), 
                              value = ifelse(is.null(input[[paste0("vsample_", i)]]),
                                             3000, input[[paste0("vsample_", i)]]), 
                              min = 1)))
      }) 
    }
  })       
  
  output$nsamples2 = renderUI({n_samples2()})
  
    
    # (assuming)
    # input$cells - number of cells
    # input$loss - % of cells lost
    # input$depth - reads per cell
    # input$rl - read length
    # input$sequencer - sequencer
    # input$phix - phix percentage
  
  # perform calculations for output data table when button is pushed
  observeEvent(input$run, { 
    
    bases_lane = seq_table[input$sequencer, "yield_gb"]*10**9*((100-input$phix)/100)
    
    # per sample calculation
    samples_df = do.call(rbind, lapply(1:input$nsamples, function(i) {
        n_cells = ifelse(input$visium,
                         input[[paste0("vsample_",i)]],
                         input[[paste0("sample_",i)]])
        est_cells = ifelse(input$visium,
                           as.double(n_cells),
                           n_cells*((100-input$loss)/100))
        reads_sample = est_cells*input$depth # this is read PAIRS
        bases_sample = reads_sample*input$rl
        
        n_lanes = bases_sample/bases_lane
        
        
        outtab = data.frame("Added" = ifelse(input$visium,
                                           ifelse(input[[paste0("vadd_",i)]],
                                                  as.character(icon("ok", lib = "glyphicon")),
                                                  " "),
                                           ifelse(input[[paste0("add_",i)]],
                                                  as.character(icon("ok", lib = "glyphicon")),
                                                  " ")),
                            "Sample name" = ifelse(input$visium,
                                 input[[paste0("vsample_name",i)]],
                                 input[[paste0("sample_name",i)]]),
                            "Estimated cells" = est_cells,
                            "Total reads" = reads_sample,
                            "# lanes" = round(n_lanes, digits = 5),
                            "# lanes (up)" = ceiling(n_lanes),
                            "# lanes (down)" = floor(n_lanes),
                            "cost (up)" = ceiling(n_lanes)*seq_table[input$sequencer, "cost_lane"],
                            "cost (down)" = floor(n_lanes)*seq_table[input$sequencer, "cost_lane"],
                            check.names = F)
      }))
    
    
    # total calculations
    total_cells = ifelse(input$visium,
                         sum(sapply(1:req(input$nsamples), function(x) ifelse(input[[paste0("vadd_",x)]],
                                                                              input[[paste0("vsample_",x)]],
                                                                              0))),
                         sum(sapply(1:req(input$nsamples), function(x) ifelse(input[[paste0("add_",x)]],
                                                                              input[[paste0("sample_",x)]],
                                                                              0))))
    est_cells = ifelse(input$visium,
                       as.double(total_cells),
                       total_cells*((100-input$loss)/100))
    reads_sample = est_cells*input$depth # this is read PAIRS
    bases_sample = reads_sample*input$rl
    
    n_lanes = bases_sample/bases_lane
    
    outtab = data.frame("Estimated cells" = est_cells,
                        "Total reads" = reads_sample,
                        "# lanes" = n_lanes,
                        "# lanes (up)" = ceiling(n_lanes),
                        "# lanes (down)" = floor(n_lanes),
                        "cost (up)" = ceiling(n_lanes)*seq_table[input$sequencer, "cost_lane"],
                        "cost (down)" = floor(n_lanes)*seq_table[input$sequencer, "cost_lane"],
                        check.names = F)
    
    downloadtab = data.frame("Sequencer" = input$sequencer,
                             "Estimated cells" = est_cells,
                             "Total reads" = reads_sample,
                             "Depth" = input$depth,
                             "Read Lenght" = input$rl,
                             "Phix (%)" = input$phix,
                             "# lanes" = n_lanes,
                             "# lanes (up)" = ceiling(n_lanes),
                             "# lanes (down)" = floor(n_lanes),
                             "cost (up)" = ceiling(n_lanes)*seq_table[input$sequencer, "cost_lane"],
                             "cost (down)" = floor(n_lanes)*seq_table[input$sequencer, "cost_lane"],
                             check.names = F)
    
    
    # output tables
    output$outTable = renderTable(outtab, width = "100%")
    output$sampleOutTable = DT::renderDT({
      DT::datatable(samples_df,
                    rownames = F,
                    escape = F)
      })
    
    # download tables
    output$downloadDataTotal <- downloadHandler(
      filename = paste0(input$project,"_sequencing_calc.csv"),
      content = function(file) {
        write.csv(downloadtab, file, row.names = F, quote=F)
      }
    )
    
    output$downloadData <- downloadHandler(
      filename = paste0(input$project,"_samples_calc.csv"),
      content = function(file) {
        write.csv(samples_df, file, row.names = FALSE)
      }
    )
    
    })
}


# Run the application 
shinyApp(ui = ui, server = server)
