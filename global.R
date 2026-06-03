library(shiny)
library(htmlwidgets)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
library(tidyverse)
library(wordcloud2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggtangle)
library(ggplot2)

##---- Load RPIs Data ----
load(file = "data/data.RData")

##---- Data Preprocessing
prepare_data_with_links <- function(df) {
  df$lncRNA_Name_link <- with(df, sprintf(
    '<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank">%s</a>',
    lncRNA_Name, lncRNA_Name
  ))
  
  df$lncRNA_RNALocate_link <- with(df, ifelse(Entry != "", 
    sprintf(
      '<a href="http://rnalocate.org/show_search?searchType=exact&dataset=Symbol&Keyword=%s&category=lncRNA&species=Homo+sapiens&sources=experiment&score1=0.0&score2=1.0" target="_blank">View</a>',
      lncRNA_Name),
    ""
  ))
  
  df$Protein_name_link <- with(df, ifelse(Entry != "", 
    sprintf('<a href="https://www.uniprot.org/uniprotkb/%s" target="_blank">%s</a>', Entry, Protein_name),
    Protein_name
  ))
  
  df$Protein_Domains_link <- with(df, ifelse(Entry != "", 
    sprintf('<a href="https://www.ebi.ac.uk/interpro/protein/UniProt/%s" target="_blank">View</a>', Entry),
    ""
  ))

  df$AlphaFoldDB_link <- with(df, ifelse(Entry != "", 
    sprintf('<a href="https://alphafold.ebi.ac.uk/entry/%s" target="_blank">View</a>', Entry),
    ""
  ))
  
  df$KEGG_link <- with(df, ifelse(KEGG != "", 
    sprintf('<a href="https://www.kegg.jp/entry/%s" target="_blank">%s</a>', KEGG, KEGG),
    ""
  ))
  
  df$Data_link <- with(df, ifelse(
    Data != "" & grepl("^GSE", Data),
    sprintf('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" target="_blank">%s</a>', Data, Data),
    ifelse(
      Data != "" & grepl("^ENCSR", Data),
      sprintf('<a href="https://www.encodeproject.org/%s/" target="_blank">%s</a>', Data, Data),
      ifelse(
        Data != "" & grepl("^[0-9]+$", Data),
        sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank">%s</a>', Data, Data),
        ""
      )
    )
  ))
  
  df
}
