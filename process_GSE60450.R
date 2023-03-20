library(GEOquery)
geo <- getGEO("GSE60450")
pData(geo[[1]])

getGEOSuppFiles("GSE60450")

counts <- readr::read_tsv("GSE60450/GSE60450_Lactation-GenewiseCounts.txt.gz") %>% 
  mutate(EntrezGeneID = as.character(EntrezGeneID))
View(counts)

library(org.Mm.eg.db)
library(dplyr)
anno <- AnnotationDbi::select(org.Mm.eg.db, keys = as.character(counts$EntrezGeneID), 
               columns = c("SYMBOL","GENENAME"))

left_join(counts, anno, by = c("EntrezGeneID" = "ENTREZID")) %>% 
  dplyr::select(EntrezGeneID, Length, SYMBOL,GENENAME,contains("MCL1")) %>% 
  filter(!duplicated(EntrezGeneID)) %>% 
  rename_at(vars(contains("MCL1")), function(x) strtrim(x,7)) %>% 
  readr::write_csv("GSE60450_Lactation_forAnalysis.csv")
  
