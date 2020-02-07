
library(SCTCwhatateam)
library(reticulate)
library(Rmagic)
# load("~/TMTB/Computational_Research/UniSA/Dream_Challenges/SCTC/DC_package/package/SCTCwhatateam/data/SCTCpreprocessedData.RData")

#============================
# Hoang 3 Sep 2019
# Add constants
# c_rlist <- "/home/buu/WhatATeam/Shiny_SCTC/Data/rlist.rda"
# c_linmethods_py <- "/home/buu/WhatATeam/Shiny_SCTC/codes/linmethods.py"
c_rlist <- "C:/Users/phavy022/MyDoc/13Dream/WhatATeam/Shiny_SCTC/Data/rlist.rda"
c_linmethods_py <- "C:/Users/phavy022/MyDoc/13Dream/WhatATeam/Shiny_SCTC/codes/linmethods.py"
#============================

make_dist_map = function(dge_normalized, genes84=NULL, thres=0.23) {
  thres=0.23
  tmp = match(genes84, rownames(dge_normalized))
  nml_extract = dge_normalized[tmp, ]
  
  threshold = apply(nml_extract, 1, function(x) {
    x = x[-which(x==0)]
    return(quantile(x, thres))
  })
  
  distmap_mat = matrix(0L, nrow=nrow(nml_extract), ncol=ncol(nml_extract))
  for (i in 1:nrow(nml_extract))
    for (j in 1:ncol(nml_extract))
      distmap_mat[i,j] = ifelse(nml_extract[i,j] >= threshold[i], 1, 0)  
  rownames(distmap_mat) = rownames(nml_extract)
  colnames(distmap_mat) = colnames(nml_extract)
  return(data.frame(distmap_mat))
}

prediction_method = function(genes, prediction, list_option, bicells84=NULL, biloc84=NULL, dge_normalized=NULL, bdtnp=NULL) {
  print(prediction)
  print(list_option)
  # list_option = "20"
  # list_option = gene_list_res
  gene_option = switch (EXPR = list_option,
                        "P5_1" = genes$ls20,
                        "P5_2" = genes$ls40,
                        "P5_3" = genes$ls60
  )
  
  writeLines("Working on prediction")
  list2 = switch(EXPR=prediction,
                 "P3_1" = ListbyMCC(bicells84 = bicells84, biloc84 = biloc84, geneList = gene_option),
                 "P3_2" = ListbyMCCLOF(bicells84 = bicells84, biloc84 = biloc84, geneList = gene_option),
                 "P3_3" = Listbycor(normcells = t(as.matrix(dge_normalized[gene_option,])), normloc = as.matrix(bdtnp[, gene_option]), method = "pearson") ################## Xiaomei
  )
  return(list2)
}

fix_names <- function(data) {
  data = as.matrix(data)
  genes = rownames(data)
  genes = gsub("-",".",genes,fixed = T)
  genes = gsub("(spl)",".spl.",genes,fixed = T)
  # genes = gsub("","'",genes,fixed = T)
  rownames(data) = genes
  return(as.data.frame(data))
}

dist_preproc = function(mat, gene_select = NULL) {
  id = match(gene_select, rownames(mat))
  ans = as.matrix(dist(mat[id, ]))
  rownames(ans) = colnames(ans) = rownames(mat)
  
  if (!is.null(gene_select))
    ans = ans[gene_select, gene_select]
  return(ans)
}
# dge_normalized = fix_names(dge_normalized)

preproc_MAGIC1 = function(data, gene_markers) {
  #============================
  # Hoang 04 Feb 20
  # Add for debug
  #saveRDS(data, file = "C:/Users/phavy022/MyDoc/Temp/data.rds")
  #saveRDS(gene_markers, file = "C:/Users/phavy022/MyDoc/Temp/gene_markers.rds")
  #============================
  after_magic = Rmagic::magic(data, gene_markers, t = 'auto', seed = 1234)
  # mat_magic = as.matrix(after_magic)
  # data_new = mat_magic[gene_markers, ]
  Eucl = as.matrix(dist(t(after_magic$result)))
  return(Eucl)
}

# genelist = genes84
runExperiments = function(bdtnp, binarized_bdtnp, dge_raw, dge_normalized, geometry, seed_list=NULL, preprocessing, gselection, prediction, gene_list_res, visualization) {
  

  genes84 = colnames(bdtnp)
  dge_binarized_distMap = make_dist_map(dge_normalized, genes84)
  bicells84 = t(dge_binarized_distMap[genes84,])
  biloc84 = binarized_bdtnp[,genes84]
  

  numofbins(biloc84,genes84)
  numofbins(bicells84,genes84)

  Eucl = switch(EXPR = preprocessing,
                "P1_1" = preproc_MAGIC1(dge_raw, genes84),
                "P1_2" = dist_preproc(t(dge_raw), genes84),
                "P1_3" = dist_preproc(dge_normalized[genes84, ])
                # "P1_2" = preproc_MAGIC_dremi(bdtnp, genes84)
  )


  if (("P2_61" %in% gselection) | ("P2_62" %in% gselection) | ("P2_63" %in% gselection)) {
    writeLines("Working on Dremi")
    dremi = preproc_MAGIC_dremi(dge_raw, genes84)
  }
  
  Eucl1 = Eucl

  # save.image("save1.RData")
  # return()
  writeLines("Working on rlist=ListbyMCC")
  start_time <- Sys.time()
  # rlist = ListbyMCC(bicells84, biloc84, genes84)
  
  #=======================
  # Hoang 3 Sep 2019
  # Update with using constant
  # load("/home/buu/WhatATeam/Shiny_SCTC/Data/rlist.rda")
  load(c_rlist)
  #=======================
  
  end_time <- Sys.time()
  print(paste("Running time:", end_time - start_time))
  # save.image("save1_1.RData")
  
  writeLines("Working on gene selection")
  start_time <- Sys.time()
  
  # if (('P2_7' %in% gselection)  | ("P2_8" %in% gselection) | ("P2_9" %in% gselection)) {
  #   write.csv(bdtnp, file= "bdtnp.csv")
  #   write.csv(dge_normalized, file = "dge_normalized.csv")
  # }
  # return(colnames(data.frame(t(dge_normalized))))
  
  #=======================
  # Hoang 3 Sep 2019
  # Update with using constant
  # source_python("/home/buu/WhatATeam/Shiny_SCTC/codes/linmethods.py")
  source_python(c_linmethods_py)
  #=======================
  # tmp = linGen(data.frame(bdtnp), data.frame(t(dge_normalized)))
  
  genes = switch(EXPR = gselection,
                 'P2_1' = selectFarGenes(Eucl),
                 'P2_2' = selectFarGenesWithSeeds(seed_list, Eucl = Eucl),
                 'P2_3' = selectHighRakingGenes(Eucl),
                 'P2_4' = selectGenesByMADExpression(dge_normalized,Eucl),
                 'P2_5' = selectGenesByMADDistance(Eucl),
                 'P2_61' = selectGenesByInfluence(dremi, Eucl1),
                 'P2_62' = selectGenesByInfluenceWithSeeds(seed_list, dremi),
                 'P2_63' = selectHighRakingGenesByInfluence(dremi, Eucl1),
                 # 'P2_7' = linRev("./bdtnp.csv", "./dge_normalized.csv"),
                 # 'P2_8' = linFwd("./bdtnp.csv", "./dge_normalized.csv"),
                 # 'P2_9' = linGen("./bdtnp.csv", "./dge_normalized.csv")
                 'P2_7' = linRev(data.frame(bdtnp), data.frame(t(dge_normalized))),
                 'P2_8' = linFwd(data.frame(bdtnp), data.frame(t(dge_normalized))),
                 'P2_9' = linGen(data.frame(bdtnp), data.frame(t(dge_normalized)))
  )

  
  # return(paste(genes$ls40, collapse=", "))

  
  end_time <- Sys.time()
  print(paste("Total running time:", end_time - start_time))
  
  
  numofbins(biloc84, genes$ls40)
  numofbins(bicells84,genes$ls40)
  
  # return("hello")
  writeLines("working on prediction")
  # load("/home/buu/WhatATeam/Shiny_SCTC/Data/list2.rda")
  list2 = prediction_method(genes, prediction, gene_list_res, bicells84, biloc84, dge_normalized, bdtnp)
  #=======================
  # Hoang 5 Feb 2020
  # Modify
  
  # save.image("save2.RData")
  #=======================
  
  # writeLines("working on QEvaluation")
  # QEvaluation(realLoc = rlist, predloc = list2, geometry = geometry, plotDist = TRUE)
  
  # visual_output = list()
  # if ("P4_1" %in% visualization)
  #   visual_output = append(visual_output, list(plot2DCellPositions(geometry = geometry, locidx = list2[,1])))
  # if ("P4_2" %in% visualization)
  #   visual_output = append(visual_output, list(plot2DgenePattern(geometry = geometry, locidx = list2[,1]))) ################## Xiaomei
  # if ("P4_3" %in% visualization)
  #   visual_output = append(visual_output, list(plot3DCellPositions(geometry = geometry, locidx = list2[,1])))
  # if ("P4_4" %in% visualization)
  #   visual_output = append(visual_output, list(plot3DgenePattern(geometry = geometry, locidx = list2[,1]))) ################## Xiaomei
  
  # dev.off()
  
  output_genes = switch(EXPR = gene_list_res,
                        "P5_1" = genes$ls20,
                        "P5_2" = genes$ls40,
                        "P5_3" = genes$ls60
  )
  # save.image("save5.RData")
  # 
  # write.table(output_genes, file = "gene_set.csv", row.names = F, col.names = F, sep ="\t")
  # write.table(list2, file = "cell_locations.txt", row.names = F, col.names = F, sep = "\t")
  # save(visual_output, file="visualization.RData")
  # zip(zipfile = paste("Result_", run_count, ".zip", sep =""), files=c("gene_set.csv", 'cell_locations.txt', 'visualization.RData'))
  # 
  # 

  writeLines("finished")
  return(list(output_genes, list2))
  # return(list(output_genes, list2, visual_output))
}
