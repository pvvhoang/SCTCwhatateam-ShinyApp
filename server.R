#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# install.packages("shinyjs")

library(reticulate)

#============================
# Hoang 29 Aug 2019
# Add constants

# c_conda <- "/home/buu/anaconda3/bin/conda"
# c_python <- "/home/buu/anaconda3/bin/python"
# c_linmethods <- "/home/buu/WhatATeam/Shiny_SCTC/codes/linmethods.py"
# c_shiny_run <- "/home/buu/WhatATeam/Shiny_SCTC/Shiny_run"
# c_run_experiment <- "/srv/shiny-server/SCTCWhatATeam/runExperiments.R"
# c_codes <- "/home/buu/WhatATeam/Shiny_SCTC/codes/"
# c_shiny_run_dash <- "/home/buu/WhatATeam/Shiny_SCTC/Shiny_run/shiny_run_"

c_conda <- "C:/Users/phavy022/AppData/Local/Continuum/anaconda3/Scripts/conda.exe"
c_python <- "C:/Users/phavy022/AppData/Local/Continuum/anaconda3/python.exe"
c_linmethods <- "C:/Users/phavy022/MyDoc/13Dream/WhatATeam/Shiny_SCTC/codes/linmethods.py"
c_shiny_run <- "C:/Users/phavy022/MyDoc/13Dream/WhatATeam/Shiny_SCTC/Shiny_run"
c_run_experiment <- "C:/Users/phavy022/MyDoc/13Dream/Script/runExperiments.R"
c_codes <- "C:/Users/phavy022/MyDoc/13Dream/WhatATeam/Shiny_SCTC/codes/"
c_shiny_run_dash <- "C:/Users/phavy022/MyDoc/13Dream/WhatATeam/Shiny_SCTC/Shiny_run/shiny_run_"
#============================

#============================
# Hoang 29 Aug 2019
# Update, using constants
# use_condaenv(condaenv = "conda_env", conda = "/home/buu/anaconda3/bin/conda")
# use_condaenv(condaenv = "rstudio", conda = c_conda)
use_condaenv(condaenv = "conda_env", conda = c_conda)
# use_python("/home/buu/anaconda3/bin/python")
use_python(c_python)
#============================

py_discover_config(required_module = "magic")

#============================
# Hoang 29 Aug 2019
# Update, using constants
# source_python("/home/buu/WhatATeam/Shiny_SCTC/codes/linmethods.py")
source_python(c_linmethods)
#============================

library(shiny)
library(ggplot2)
library(SCTCwhatateam)

data_choices = c("Use MAGIC" = "P1_1", 
                 "Use Raw file" = "P1_2",
                 "Use Normalized file" = "P1_3",
                 
                 "Select Genes with highest variations" = "P2_1", "Select Genes with highest variations with Seeds" = "P2_2", 
                 "Select Genes by high ranking from Google PageRank" = "P2_3", "Select Genes by MAD Expression" = "P2_4", 
                 "Select Genes by MAD Distance" = "P2_5","Select Genes by Influence" = "P2_61",
                 "Select Genes by Influence With Seeds" = "P2_62","Select Genes by high ranking by Influence" = "P2_63",
                 "Select Genes by LinRev" = "P2_7","Select Genes by LinFwd" = "P2_8", "Select Genes by LinGen" = "P2_9",
                 
                 "By MCC" = "P3_1", "By MCC-LOF" = "P3_2", "By Correlation" = "P3_3",
                 
                 "2D Cell Positions" = "P4_1", "plot2DgenePattern" = "P4_2", 
                 "plot3DCellPositions" = "P4_3", "plot3DgenePattern" = "P4_4",
                 
                 "20 genes" = "P5_1", "40 genes" = "P5_2", "60 genes" = "P5_3"
)


run_count_f = function() {
  # writeLines("yes")
  
  check_file = file.exists("count_run_shiny.txt")
  if (check_file == F) {
    sink("count_run_shiny.txt")
    writeLines("1")
    sink()
  }
  
  run_count = read.table("count_run_shiny.txt")[,1]
  
  
  return(run_count)
}
# options(shiny.sanitize.errors = FALSE)
# rsconnect::showLogs()
# option = c()

shinyServer(function(input, output) {
  options(shiny.maxRequestSize=100*1024^2) 
  output$cont1 = renderText({  paste("Preprocessing Method:", names(data_choices)[which(data_choices == input$preprocessing)])})
  output$cont2 = renderText({  paste("Gene Selection Method:", names(data_choices)[which(data_choices == input$gselection)])})
  output$cont3 = renderText({  paste("Prediction Method:", names(data_choices)[which(data_choices == input$prediction)])})
  output$cont4 = renderText({  
    print(input$gene_list_res) 
    tmp =  paste(names(data_choices)[match(input$gene_list_res, data_choices)], collapse  =", ")
    paste("Gene list:", tmp)})
  
  # output$cont5 = renderText({  
  #   print(input$visualization) 
  #   tmp =  paste(names(data_choices)[match(input$visualization, data_choices)], collapse  =", ")
  #   paste("Visualization Method (NOT required):", tmp)})
  

  
  # output$cont_email = renderText({  paste("Results will be sent to:", input$email)})
  
  observeEvent(
    eventExpr = input[["submit_bt"]],
    handlerExpr = {
      # setwd("/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/SCTC/SCTC_shiny/codes")
      #=======================
      # Hoang 3 Sep 2019
      # Update with using constant
      # setwd("/home/buu/WhatATeam/Shiny_SCTC/Shiny_run")
      setwd(c_shiny_run)
      #=======================
      
      # run_count = 1
      check_file = file.exists("count_run_shiny.txt")
      if (check_file == F) {
        sink("count_run_shiny.txt")
        writeLines("1")
        sink()      
      }
      
      run_count = read.table("count_run_shiny.txt")[,1]
  


      print(run_count)
      dir.create(paste("shiny_run_", run_count, sep =""))
      print(getwd())
      
      sink("./count_run_shiny.txt")
      writeLines(as.character(run_count+1))
      print(run_count+1)
      sink()
      setwd(paste("./shiny_run_", run_count, sep =""))
      
      
      id = showNotification(
        paste("Thanks for your request! Please wait until the program finish. You will see 'Finish' Notification. Please take note your process ID is ", run_count, sep=""),
        duration = 2000,
        type = "error"
      )
      # removeUI(selector = "#submit_bt", immediate = T)
      # stopApp()
      
      
      ########################### APP PROCESSING #################################
      # setwd("/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/SCTC/SCTC_shiny/SCTCWhatATeam")
      #=======================
      # Hoang 3 Sep 2019
      # Update with using constant
      # source("/srv/shiny-server/SCTCWhatATeam/runExperiments.R") 
      source(c_run_experiment) 
      #=======================
      


      bdtnp = binarized_bdtnp = dge_normalized = geometry = binarized_bdtnp = seed_list = NULL
      # bdtnp = ifelse(!is.null(input$bdtnp$datapath), read.csv(input$bdtnp$datapath, row.names = 1), NULL)
      # binarized_bdtnp = ifelse(!is.null(input$binarized_bdtnp$datapath), read.csv(input$binarized_bdtnp$datapath, row.names = 1), NULL)
      # dge_normalized = ifelse(!is.null(input$dge_normalized$datapath), read.csv(input$dge_normalized$datapath, row.names = 1), NULL)
      # geometry = ifelse(!is.null(input$geometry$datapath), read.csv(input$geometry$datapath, row.names = 1), NULL)
      # binarized_bdtnp = ifelse(!is.null(input$binarized_bdtnp$datapath), read.csv(input$binarized_bdtnp$datapath, row.names = 1), NULL)
      # seed_list = ifelse(!is.null(input$seed_list$datapath), read.csv(input$seed_list$datapath, stringsAsFactors = F), NULL)
      
      # seed_list = read.csv(input$seed_list$datapath, stringsAsFactors = F)
      
      
      
      # setwd("/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/SCTC/SCTC_shiny/codes")
      #=======================
      # Hoang 3 Sep 2019
      # Update with using constant
      # setwd("/home/buu/WhatATeam/Shiny_SCTC/codes/")
      setwd(c_codes)
      #=======================

      bdtnp = read.csv(input$bdtnp$datapath, row.names = 1)
      binarized_bdtnp = read.csv(input$binarized_bdtnp$datapath, row.names = 1)
      dge_normalized = read.csv(input$dge_normalized$datapath, row.names = 1)
      dge_raw = read.csv(input$dge_raw$datapath, row.names = 1)
      geometry = read.csv(input$geometry$datapath, row.names = 1)
      seed_list = read.csv(input$seed_list$datapath, stringsAsFactors = F)[,1]

      


      # bdtnp = read.csv("../Data/bdtnp.csv", row.names = 1)
      # binarized_bdtnp = read.csv("../Data/binarized_bdtnp.csv", row.names = 1)
      # dge_normalized = read.csv("../Data/dge_normalized1.csv", row.names = 1)
      # 
      # # dge_normalized1 = read.csv("../Data/dge_normalized1.csv", row.names = 1)
      # # write.csv(t(dge_normalized), file="../Data/dge_normalized1.csv")
      # geometry = read.csv("../Data/geometry.csv", row.names = 1)
      # seed_list = read.csv("../Data/seed_list.csv", stringsAsFactors = F)[,1]
      # preprocessing = "P1_3"
      # gselection = "P2_7"
      # prediction = "P3_1"
      # gene_list_res = "P5_1"
      # visualization = "P4_1"
      # run_count = 1

      dge_normalized = t(dge_normalized)
      preprocessing = input$preprocessing
      gselection = input$gselection
      prediction = input$prediction
      gene_list_res = input$gene_list_res
      # visualization = input$visualization

      
      
      #=======================
      # Hoang 3 Sep 2019
      # Update with using constant
      # setwd(paste("/home/buu/WhatATeam/Shiny_SCTC/Shiny_run/shiny_run_", run_count, sep =""))
      setwd(paste(c_shiny_run_dash, run_count, sep =""))
      #=======================

      # check = rep("F", 10)
      # names(check) = c('bdtnp', "binarized_bdtnp", "dge_normalized", 'geometry', 'seed_list', 'preprocessing', 'gselection', 'prediction', 'gene_list_res', 'visualization')
      # if (!is.null(bdtnp)) check[1] = "T"
      # if (!is.null(binarized_bdtnp)) check[2] = "T"
      # if (!is.null(dge_normalized)) check[3] = "T"
      # if (!is.null(geometry)) check[4] = "T"
      # if (!is.null(seed_list)) check[5] = "T"
      # check[6] = preprocessing
      # check[7] = gselection
      # check[8] = prediction
      # check[9] = gene_list_res
      # check[10] = visualization

      # write.table(data.frame(check), file="check_input.txt", sep ="\t", row.names=F, col.names=F, quote=F)

      # expri = list(1,2,3)
      expri = runExperiments(bdtnp, binarized_bdtnp, dge_raw, dge_normalized, geometry, seed_list, preprocessing, gselection, prediction, gene_list_res, visualization)
      
      # id = showNotification(
      #   paste(expri, collapse=" "),
      #   duration = 120,
      #   type = "error"
      #   )

      # print("Finished all experiments")
      # 
      # load("/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/SCTC/SCTC_shiny/codes/shiny_run_4/save5.RData")
      # load("/home/buu/WhatATeam/Shiny_SCTC/save5.RData")
      # expri = list(output_genes, list2, visual_output)
      # write.table(expri[[1]], file = "gene_set.txt", row.names = F, col.names = F, sep ="\t", quote=F)
      
      
      # write.table(expri[[1]], file = "gene_set.txt", row.names = F, col.names = F, sep ="\t", quote=F)
      # write.table(expri[[2]], file = "cell_locations.txt", row.names = F, col.names = F, sep = "\t", quote=F)
      
      filename=paste("Result_Process_ID_", run_count, ".zip", sep ="")
      # zip(zipfile = filename, files=c("gene_set.txt", 'cell_locations.txt'))
      
      # writeLines("test")
      # output$downloadData = renderUI({
      #   downloadButton('downloadData01', label = "Download results")
      # })
      
      
      id = showNotification(
        paste("FINISH!", sep=""),
        duration = 2000,
        type = "error"
      )
      
      output$downloadData <- downloadHandler(
        filename = function() {
          #=======================
          # Hoang 5 Feb 2020
          # Modify
          
          # print(paste("Result_Process_ID_", run_count, ".zip", sep =""))
          
          paste("Result_Process_ID_", run_count, ".tar", sep ="")
          #=======================
        },
        content = function(filename)  {
          # write.csv("hello", file=filename)
          # writeLines("hello", filename)
          # writeLines(expri[[1]], file = filename, row.names = F, col.names = F, sep ="\t", quote=F)
          write.table(data.frame(expri[[1]]), file = "gene_set.txt", row.names = F, col.names = F, sep ="\t", quote=F)
          write.table(expri[[2]], file = "cell_locations.txt", row.names = F, col.names = F, sep = "\t", quote=F)
          # tmp = expri[[3]]
          # tmp
          # save(tmp, file="visualization.RData")
          # print(head(expri[[1]]))
          # print(getwd())
          # ggsave("visualization.pdf")
          # filename = paste("Result_Process_ID_", run_count, ".zip", sep ="")
          # zip(zipfile = filename, files=c("gene_set.txt", 'cell_locations.txt', "visualization.pdf"))
          # zip(zipfile = filename, files=c("gene_set.txt", 'cell_locations.txt', 'visualization.pdf'))
          #=======================
          # Hoang 5 Feb 2020
          # Modify
          
          # zip(zipfile = filename, files=c("gene_set.txt", 'cell_locations.txt'))
          
          tar(paste("Result_Process_ID_", run_count, ".tar", sep =""), ".")
          file.copy(paste("Result_Process_ID_", run_count, ".tar", sep =""), filename)
          #=======================
          
          #=======================
          # Hoang 5 Feb 2020
          # Modify
          
        # }, contentType = "application/zip"
          
        }
        #=======================
        # }, contentType = "text/plain"
        # }
      )
    }, autoDestroy = T
  )
  
})

# runApp(list(ui=ui,server=server),launch.browser=T)