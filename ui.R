
# Define UI for application that draws a histogram
library(shinythemes)
# install.packages("shinythemes")

shinyUI(fluidPage(

  # Application title
  # tags$head(
  #   tags$style(
  #     ".title {margin: auto; width: 1000px}"
  #   )
  # ),
  # tags$div(class="title", titlePanel("Indentify Cell locations from gene markers - Single Cell Transcriptomic Challenge")),
  # column(3, offset = 4, titlePanel("Indentify Cell locations from gene markers - Single Cell Transcriptomic Challenge")),
  # titlePanel(
  #   ("Indentify cell locations by scRNA-Seq data with gene markers")
  #   # h2("Please find the manual here: https://github.com/thanhbuu04/SCTCwhatateam")
  #   ),

  # tags$hr(),

  fluidRow(
    mainPanel(width = 10,
              h1(strong("Indentify cell locations by scRNA-Seq data with gene markers"), align = "center"),
              h2("Buu Truong, Xiaomei Li, Hoang VV. Pham, Thin Nguyen, Thuc D. Le", align = "center"),
              # h3("Manual and Tutorial"),
              h3(tags$a(href="https://www.dropbox.com/s/abmpfmsn0y1eaqi/Shiny%20app%20tutorial1.docx?dl=0", "Vignette"), align = "center")
    )),

  tags$hr(),

  fluidRow(
    column(5,
           wellPanel(
             titlePanel("Input data (.csv)"),
             fileInput("bdtnp", "Raw reference gene expression for locations (Locations x Genes)",accept = c(".csv")),
             hr(),
             fileInput("binarized_bdtnp", "Binarized reference gene expression for locations (Locations x Genes)", multiple = F, accept = c(".csv")),
             hr(),
             fileInput("dge_raw", "scRNA-Seq raw file (Cells x Genes)", multiple = F, accept = c(".csv")),
             hr(),
             fileInput("dge_normalized", "scRNA-Seq normalized file (Cells x Genes)", multiple = F, accept = c(".csv")),
             hr(),
             fileInput("geometry", "3D coordinates for each locations (Locations x (x,y,z) coordinates)", multiple = F, accept = c(".csv")),
             hr(),
             fileInput("seed_list", "Seed list (Locations x Coordinates)", multiple = F, accept = c(".csv")))),

    column(6,
           wellPanel(
             titlePanel("Preprocessing"),
             selectInput("preprocessing",
                         "",
                         choices = list("Use MAGIC" = "P1_1",
                                        "Use Raw file" = "P1_2",
                                        "Use Normalized file" = "P1_3"
                         ),
                         selected = "P1_1"))),

    column(6,
           wellPanel(
             titlePanel("Gene Selection Method"),
             selectInput("gselection", "",
                         choices =
                           list("Select Genes with highest variance" = "P2_1", "Select Genes with highest variance with Seeds" = "P2_2",
                                "Select Genes by high ranking from Google PageRank" = "P2_3", "Select Genes by MAD Expression" = "P2_4",
                                "Select Genes by MAD Distance" = "P2_5","Select Genes by Influence" = "P2_61",
                                "Select Genes by Influence With Seeds" = "P2_62","Select Genes by high ranking by Influence" = "P2_63",
                                "Select Genes by Efficient Reverse Stepwise Linear Regression" = "P2_7","Select Genes by Efficient Forward Stepwise Linear Regression" = "P2_8", "Select Genes by Efficient General Stepwise Linear Regression" = "P2_9"),
                         selected = "P2_1"))),
    column(6,
           wellPanel(
             titlePanel("Prediction"),
             selectInput("prediction", "",
                         choices =
                           list("By MCC" = "P3_1", "By MCC-LOF" = "P3_2", "By Correlation" = "P3_3"),
                         selected = "P3_1"))),

    column(6,
           wellPanel(
             titlePanel("Get set genes"),
             selectInput("gene_list_res", "",
                         choices =
                           list("20 genes" = "P5_1", "40 genes" = "P5_2",
                                "60 genes" = "P5_3"),
                         selected = "P5_1")))

    # column(6,
    #        wellPanel(
    #          titlePanel("Visualization"),
    #          checkboxGroupInput("visualization", "",
    #                             choices =
    #                               list("2D Cell Positions" = "P4_1"
    #                                    # "2D Gene Pattern" = "P4_2",
    #                                    # "3D Cell Positions" = "P4_3"
    #                                    # "2D Gene Pattern" = "P4_4"
    #                               ))))

  ),

  tags$hr(),


  # fluidRow(
  #   column(5,
  #          wellPanel(
  #            textInput("email", "Email")
  #          ))
  # ),
  #
  # tags$hr(),
  #
  fluidRow(
    titlePanel("Summary"),
    mainPanel(width = 10,
              textOutput("cont1"),
              textOutput("cont2"),
              textOutput("cont3"),
              textOutput("cont4")
              # textOutput("cont5")
              # textOutput("cont_email")
    )),

  tags$hr(),
  actionButton(
    inputId = "submit_bt",
    label = "Submit"
  ),

  tags$hr(),
  verbatimTextOutput("console"),

  downloadButton("downloadData", "  Download", icon("download"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  # uiOutput("downloadData")

  # div(style="display:inline-block",, style="float:right")
  # actionButton("submit_button", "Submit")
))

# Run the application
# shinyApp(ui = ui, server = server)
# runApp(list(ui=ui,server=server), host = "130.220.238.199", port=3838,launch.browser=T)
# runApp(list(ui=ui,server=server),launch.browser=T)
# runApp(launch.browser = TRUE)
# runApp(list(ui=ui,server=server),launch.browser=T)

