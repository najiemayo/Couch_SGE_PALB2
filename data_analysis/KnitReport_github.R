library(rmarkdown)
### note: some directory names were removed for security reason.
setwd("###/BRCA2/rpgm/") ## where the rmd file reside
render_report = function(Gene, Exon, Output_dir, loessSubset, EventCount_Filter, vset, C_fc) {
  rmarkdown::render(
    "Ind_Exon.Rmd", params = list(
      Exon = Exon, 
      Output_dir = Output_dir, 
      loessSubset = loessSubset, 
      EventCount_Filter = EventCount_Filter,
      vset = vset,
      Gene = Gene,
      C_fc = C_fc
    ),
    output_file = paste0(Output_dir, Gene, "_E", Exon, "_",  vset, "_", loessSubset, "_", EventCount_Filter,  ".html")
  )
}

render_report2 = function(Gene, Exons, Input_dir, Output_dir, loessSubset, EventCount_Filter, vset, Score_method, Lab_curation, Ratio, MM, Outlierremove, clinvarfile, extraFilter = "No") {
  #define main directory
  main_dir <- paste0("###results/11_2023/", Output_dir)
  
  #create directory if it doesn't exist
  if (!dir.exists(main_dir)) {
    dir.create(main_dir)
  }
  
  #define sub directory
  sub_dir <- "reports"
  
  #define full directory
  full_dir <- file.path(main_dir, sub_dir)
  
  #create directory if it doesn't exist
  if (!dir.exists(full_dir)) {
    dir.create(full_dir)
  }
  
  rmarkdown::render(
    "norm_classification_parm.Rmd", params = list(
      Input_dir = Input_dir,
      Output_dir = main_dir,
      loessSubset = loessSubset, 
      EventCount_Filter = EventCount_Filter,
      vset = vset,
      Score_method = Score_method,
      Lab_curation = Lab_curation,
      Exons = Exons,
      Ratio = Ratio,
      Gene = Gene,
      MM = MM,
      Outlierremove = Outlierremove,
      clinvarfile = clinvarfile,
      extraFilter = extraFilter
    ),
    output_file = paste0( paste0(main_dir,"/reports/" ), 
                          Gene, "_", vset, "_", loessSubset, "_", EventCount_Filter, "_FS", Score_method, 
                          "Lab_curation_", Lab_curation, "_E", gsub(", ", "", Exons),"_", Ratio,
                          "_MM_", MM,
                          "_Outlierremove_", Outlierremove,"_extraFilter_", extraFilter, ".html")
  )
}

Gene <- "BRCA2"
Es <-
  c(
    "15",
  "16",
  "17", "18C",
  "18N",  "19",
  "20_sg1",
  "21_sg4", "22_sg5",
  "23_sg4_3",
  "24_sg3",
  "25C_sg5",
  "25N_sg3", 
  "26_sg4_3")
loessSubsetFiters <-  c("ModFindlay2")
EventCount_Filters <- c("libD510")
Score_methods <- c("C")
vsets <- c("SNV")
Lab_curations <- c("No")
for(Exon in Es){
  for(loessSubset in loessSubsetFiters){
    for(EventCount_Filter in EventCount_Filters){
      for(vset in vsets){
        render_report(Exon= Exon, 
                      Output_dir = "####/Ind_full/", 
                      loessSubset = loessSubset, 
                      EventCount_Filter = EventCount_Filter,
                      vset = vset,
                      Gene = Gene,
                      C_fc = 1) 
      }
    }
  }
}


## after get all the reports, run this to get the one metric file
##
setwd("/research/bsi/projects/breast/s108235.tripneg_Couch/projects/s108235.tripneg/breast_requests/Predisposition_genes/Chunling/BRCA2/results/11_2023/Ind/")
library(data.table)
mylist <- list()
for (fn in list.files(pattern = glob2rx("BRCA2_E*_metrics.csv"), full.names = TRUE)){
  fi <- read.csv(fn)
  mylist[[length(mylist)+1]] <- fi 
}

df <- rbindlist(mylist, use.names = TRUE)
write.csv(df, "Exon_metrics.csv", row.names = FALSE)

## cross exons
Gene <- "BRCA2"
Ratios <- c("D14vLib", "D14vD5")
loessSubsetFiters <-  c( "ModFindlay2")
EventCount_Filters <- c("libD510")
Score_methods <-"C" 
vsets <- c("SNV")
Lab_curations <- c("No")
MMs <- c("Normal") 
Outlierremoves <- c("No") 
exfs <- c("Yes") 
fcs <- c(1) 

setwd("####/BRCA2/rpgm/")
clinvardir <- "###/data/10_2023_curation/"
dtdir <- "###/results/11_2023/"

exonsset <- c("15, 16, 17, 18C, 18N, 19, 20, 21, 22, 24, 25C, 25N",
              "15, 16, 17, 18C, 18N, 19, 20, 21, 22, 23, 24, 25C, 25N, 26")

for(exons in exonsset){
  fn <- ifelse(exons == "15, 16, 17, 18C, 18N, 19, 20, 21, 22, 23, 24, 25C, 25N, 26", "1526",
               ifelse(exons == "15, 16, 17, 18C, 18N, 19, 20, 21, 22, 24, 25C, 25N", "15222425",
                      ifelse(exons == "15, 16, 17, 18C, 18N, 19", "1519",
                             ifelse(exons == "23", "23",
                                    ifelse(exons == "26", "26", 
                                           ifelse(exons == "23, 26", "2326", "15202425"))))))
  for(loessSubset in loessSubsetFiters){
    for(EventCount_Filter in EventCount_Filters){
      for(Score_method in Score_methods){
        for(vset in vsets){
          for(Lab_curation in Lab_curations){
            for(Ratio in Ratios){
              for(MM in MMs){
                for(Outlierremove in Outlierremoves){
                  for(fc in fcs){
                    for(exf in exfs){
                      
                      render_report2( 
                        Input_dir = paste0(dtdir, "Ind_full/"), 
                        Output_dir = paste0("xexon_full", fn, "_fc", fc, "_exf_", exf, "_nonsense_syn_missense_exclude_single/"), 
                        loessSubset = loessSubset, 
                        EventCount_Filter = EventCount_Filter,
                        vset = vset,
                        Score_method = Score_method,
                        Lab_curation = Lab_curation,
                        Exons = exons,
                        Ratio  = Ratio,
                        Gene = Gene,
                        MM = MM,
                        Outlierremove = Outlierremove,
                        clinvarfile = paste0(clinvardir, "/ClinVar_BRCA2_DBD_control_111523_nonsense_syn_missense_exclude_single.txt"),
                        extraFilter = exf)
                      
                    }
                  }
                }
              }
              
            }
          }
        }
      }
    }
  }
}

