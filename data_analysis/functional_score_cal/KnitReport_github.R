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


# PALB2
Gene <- "PALB2"
Es <- c("10", "10_2")
for(Exon in Es){
  for(loessSubset in loessSubsetFiters){
    for(EventCount_Filter in EventCount_Filters){
      for(vset in vsets){
        render_report(Exon= Exon, 
                      Output_dir = "/research/bsi/projects/breast/s108235.tripneg_Couch/projects/s108235.tripneg/breast_requests/Predisposition_genes/Chunling/PALB2/results/11_2023/Ind/", 
                      loessSubset = loessSubset, 
                      EventCount_Filter = EventCount_Filter,
                      vset = vset,
                      Gene = Gene) 
      }
    }
  }
}
