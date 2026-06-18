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
library(igraph)
library(waiter)
library(shadowtext)

theme_set(theme_bw(base_family = "DejaVu Sans"))
options(ggplot2.device = "png")

##---- Load RPIs Data ----
load(file = "data/data.RData")
load(file = "data/lncRNA_bed_data.RData")
load(file = "data/protein_binding_data.RData")

RNA_list <- unique(data$lncRNA_Name)
RNA_list <- sort(RNA_list)
RNA_list <- c("", RNA_list)
protein_list <- unique(data$Protein_name)
protein_list <- sort(protein_list)
protein_list <- c("", protein_list)
cell_list <- unique(data$Cell_Line)
cell_list <- sort(cell_list)
cell_list <- c("", cell_list)

##---- Data Preprocessing ----
prepare_data_with_links <- function(df) {
  if (!"Entry" %in% colnames(df)) df$Entry <- ""
  if (!"KEGG" %in% colnames(df)) df$KEGG <- ""
  
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
  
  generate_single_link <- function(data_id) {
    if (data_id == "" || is.na(data_id)) return("")
    
    if (grepl("^GSE", data_id)) {
      return(sprintf('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" target="_blank">%s</a>', data_id, data_id))
    } else if (grepl("^ENCSR", data_id)) {
      return(sprintf('<a href="https://www.encodeproject.org/%s/" target="_blank">%s</a>', data_id, data_id))
    } else if (grepl("^[0-9]+$", data_id)) {
      return(sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s/" target="_blank">%s</a>', data_id, data_id))
    } else {
      return(data_id) 
    }
  }
  
  df$Data_link <- sapply(df$Data, function(x) {
    if (x == "" || is.na(x)) return("")
    
    if (grepl(",", x)) {
      data_items <- trimws(unlist(strsplit(x, ",")))
      links <- sapply(data_items, generate_single_link)
      return(paste(links, collapse = ", "))
    } else {
      return(generate_single_link(x))
    }
  })
  
  df
}


generate_wordcloud <- function(data) {
  protein_freq <- as.data.frame(table(data$Protein_name))
  colnames(protein_freq) <- c("word", "freq")
  
  if (nrow(protein_freq) > 100) {
    protein_freq <- protein_freq[order(-protein_freq$freq), ][1:100, ]
  }
  
  wordcloud2(protein_freq, size = 0.8,
             color = "random-light", backgroundColor = "white")
}


########################### UI ###########################
ui <- shinyUI(
  fluidPage(
    ##-- Favicon ----
    tags$head(
      tags$link(rel = "shortcut icon", href = "img/logo.ico"),
      tags$link(rel="stylesheet", type = "text/css",
                href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"),
      tags$link(rel="stylesheet", type = "text/css",
                href = "https://fonts.googleapis.com/css?family=Open+Sans|Source+Sans+Pro"),
      tags$script(HTML("
        $(document).on('click', '.protein-network-link', function(e) {
          e.preventDefault();
          var protein = $(this).data('protein');
          Shiny.setInputValue('selected_protein_network', protein, {priority: 'event'});
        });
      "))
    ),
    ##-- Logo ----
    list(tags$head(HTML('<link rel="icon", href="img/logo.png",type="image/png" />'))),
    div(style="padding: 1px 0px; width: '50%'",
        titlePanel(
          title ="", windowTitle = "HuRInterDB"
        )
    ),

    use_waiter(),
  
    waiter_show_on_load(
      html = tagList(
        div(style = "text-align: center;",
            div(class = "loader", style = "border: 5px solid #f3f3f3; border-top: 5px solid #808080; border-radius: 50%; width: 50px; height: 50px; animation: spin 1s linear infinite; margin: 0 auto;"),
            br(),
            h3("Loading......", style = "color: #808080;")
        )
      ),
      color = "rgba(255,255,255,0.9)"
    ),
    
    ##-- Header ----
    navbarPage(title = div(img(src = "img/logo.png", height = "50px"), style = "padding-left:40px;"),
              id = "navbar", 
              selected = "home",  
              theme = "styles.css", 
              fluid = TRUE,
              
###------ Home  ------###
              tabPanel(title = "Home", value = "home", icon = icon("house"),
                     br(), hr(), 
                     tags$style("
                            .search-box {background-color: #34b3db; padding: 10px;border-radius: 5px; color: white; margin: 8px 0;}
                            .stat-card {background-color: #b0cd97; padding: 20px;border-radius: 8px; color: white; text-align: center;
                                          height: 100px; display: flex; flex-direction: column; justify-content: center; margin: 10px 0;}
                            .resource-container {background-color: #f0f0f0; padding: 20px; border-radius: 8px; margin: 10px 0;}
                            .db-link {font-size: 16px; color: #0066cc; text-decoration: underline;}
                            .db-link:hover {color: #003399; font-weight: bold;}
                     "),
                     ##----- Head -------
                     HTML("<h1><center>Welcome to <b>HuRInterDB</b>!</center></h1>"),
                     ##----- Search RNA-------
                     fluidRow(
                       column(6,
                              div(class = "search-box", h3("Search lncRNA in Whole Database:"),
                                   searchInput(inputId = "RNA_search", label = NULL, placeholder = "HOTAIR", 
                                                 btnSearch = icon("search"), btnReset = icon("remove"), width = "100%")
                              )
                       ),
                       ##----- Search Protein-------
                       column(6,
                              div(class = "search-box", h3("Search Protein in Whole Database:"),
                              searchInput(inputId = "protein_search", label = NULL, placeholder = "CTCF", 
                                          btnSearch = icon("search"), btnReset = icon("remove"), width = "100%")
                              )
                       )
                     ),
                     ##-- Statistics ----
                     fluidRow(
                       column(4, div(class = "stat-card",
                                     tags$i(class = "fas fa-dna fa-2x", style = "margin-bottom: 1px;"),
                                     h2("55,925 lncRNA", style = "margin: 0; font-size: 24px;"))
                       ),
                       column(4, div(class = "stat-card",
                                     tags$i(class = "fas fa-project-diagram fa-2x", style = "margin-bottom: 1px;"),
                                     h2("7,864 Protein", style = "margin: 0; font-size: 24px;"))
                       ),
                       column(4, div(class = "stat-card",
                                     tags$i(class = "fas fa-atom fa-2x", style = "margin-bottom: 1px;"),
                                     h2("2,159,916 Interaction", style = "margin: 0; font-size: 24px;"))
                       )
                     ),
                     ## ----  Resources  ----   
                     h3("Resources"),
                     fluidRow(
                       column(12, div(class = "resource-container",
                            HTML("<h4><b>HuRInterDB</b> is a specialized database dedicated to human lncRNA-protein interactions, providing a comprehensive and curated resource for experimentally validated interactions between lncRNAs and proteins. 
                            It serves as a valuable tool for researchers investigating the molecular mechanisms underlying lncRNA functions in cellular processes and diseases. 
                            In addition to <b>HuRInterDB</b>, numerous other databases have been developed to support lncRNA-related research, offering diverse data types such as lncRNA expression profiles, functional annotations, disease associations, and interaction networks. 
                            These resources collectively enhance our understanding of the complex regulatory roles of lncRNAs in biological systems.
                            </h4>"),
                            lapply(1:3, function(row) {
                            fluidRow(
                              lapply(1:4, function(col) {
                                   idx <- (row - 1) * 4 + col
                                   db_list <- list(
                                     list(name = "FANTOM", url = "https://fantom.gsc.riken.jp"),
                                     list(name = "RNAcentral", url = "https://rnacentral.org/"),
                                     list(name = "LncRNADisease", url = "http://www.rnanut.net/lncrnadisease/"),
                                     list(name = "Lnc2Cancer", url = "http://bio-bigdata.hrbmu.edu.cn/lnc2cancer/"),
                                     list(name = "Mfold", url = "http://www.unafold.org/mfold/applications/rna-folding-form.php"),
                                     list(name = "RNAfold", url = "http://rna.tbi.univie.ac.at/"),
                                     list(name = "Rfam", url = "https://rfam.org/"),
                                     list(name = "Pfam", url = "http://pfam.xfam.org/"),
                                     list(name = "KEGG", url = "https://www.kegg.jp"),
                                     list(name = "Gene Ontology", url = "http://geneontology.org"),
                                     list(name = "OMIM", url = "https://omim.org"),
                                     list(name = "BioGRID", url = "https://thebiogrid.org")
                                   )
                                   
                                   if (idx <= length(db_list)) {
                                     item <- db_list[[idx]]
                                     column(3,
                                            a(item$name, href = item$url, 
                                              target = "_blank", class = "db-link",
                                              style = "display: block; margin: 8px 0;")
                                     )
                                   } else {
                                     column(3)
                                   }
                              })
                            )
                            })
                       )
                     )
              )
),

###------ Search  ------###
              tabPanel(title = "Search", value = "searchs", icon = icon("magnifying-glass"),
                     br(), hr(),
                     column(width = 10,style = "padding-top: 0px;",
                            column(2, selectInput(inputId = "filter_gene", label = "lncRNA", choices = RNA_list)),
                            column(2, selectInput(inputId = "filter_protein", label = "Protein", choices = protein_list)),
                            column(2, selectInput(inputId = "filter_cellline", label = "Cell Line", choices = cell_list)),
                            column(2, selectInput(inputId = "filter_method", label = "Method", 
                                                        choices = c("","RNA Pulldown","ChIRP-MS","RAP-MS","HyPR-MS","HPLC-MS","CARPID","SILAC-MS","TREX",
                                                               "RIP-seq","eCLIP-seq","CLIP-seq","PAR-CLIP","HITS-CLIP","LACE-seq","ARTR-seq",
                                                               "PRIM-seq")))
                     ),
                     column(width = 2, style = "padding-top: 55px;",
                            actionBttn(inputId = "apply_filter", label = "Select", style = "fill", 
                                          color = "success", icon = icon("check"), size = "sm") 
                     ),
                     column(width = 12, 
                            titlePanel(h3("Search results:", style = "color: #0277bd; text-align: left;")) ,
                            withSpinner(dataTableOutput("result_table"), type = 6, color = "#0277bd"),
                            downloadButton("download_data", "Download Results"),
                            br(), br(),
                            conditionalPanel(
                                condition = "output.protein_network_visible == true",
                                column(width = 12,
                                    h3("🔗 LncRNA Interaction Network for Selected Protein:", style = "color: #0277bd;"),
                                    withSpinner(plotOutput("protein_lncRNA_network", height = "400px", width = "100%"), type = 6, color = "#0277bd"),
                                    downloadButton("download_protein_network", "Download Network (PNG)", class = "btn-sm btn-primary")
                                )
                            )
                     )
              ),
              
###------ Analysis  ------###
              tabPanel(title = "RPIanalysis", value = "RPIanalysis", icon = icon("atom"),
                     br(), hr(), 
                     tags$head(
                       tags$style(HTML("
                            .inputBox {
                            background-color: #E1F5FE; 
                            padding: 10px; 
                            margin-bottom: 10px;
                            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                            }
                            .resultBox {
                            background-color: #E1F5FE; 
                            padding: 10px; 
                            border-radius: 10px; 
                            margin-top: 10px; 
                            margin-bottom: 10px;
                            }
                       "))
                     ),
                     ##---- Head ----
                     titlePanel(h2("🔍 On-line lncRNA-protein interaction analysis", style = "color: #0277bd; text-align: center;")),
                     ##---- Input ----                                   
                     column(width = 10,style = "padding-top: 0px;",
                            column(2, selectInput(inputId = "lncrna_input", label = "lncRNA", choices = RNA_list)),
                            column(2, selectInput(inputId = "cellline_input", label = "Cell Line", choices = cell_list)),
                            column(2, selectInput(inputId = "method_input", label = "Method", 
                                                 choices = c("","RNA Pulldown","ChIRP-MS","RAP-MS","HyPR-MS","HPLC-MS","CARPID","SILAC-MS","TREX",
                                                               "RIP-seq","eCLIP-seq","CLIP-seq","PAR-CLIP","HITS-CLIP","LACE-seq","ARTR-seq",
                                                               "PRIM-seq")))
                     ),
                     column(width = 2, style = "padding-top: 55px;",
                            actionBttn(inputId = "analyze_btn", label = "Continue", style = "fill", 
                                   color = "success", icon = icon("arrow-right"), size = "sm") 
                     ),
                     
                     ##-- Outputs ----
                     # Fig result
                     column(width = 12,
                            div(class = "resultBox", 
                                h3("📊 Interaction Analysis Results"),
                                tabsetPanel(
                                  tabPanel("Protein Wordcloud", 
                                           withSpinner(wordcloud2Output("wordcloud", width = "100%", height = "400px"), type = 6, color = "#0277bd"),
                                           downloadButton("download_wordcloud", "Download Wordcloud", class = "btn-sm btn-primary")),
                                  tabPanel("RBP binding", 
                                           withSpinner(plotOutput("lolliplot", width = "100%", height = "400px"), type = 6, color = "#0277bd"),
                                           downloadButton("download_lolliplot", "Download RBP binding(PNG)", class = "btn-sm btn-primary")),
                                  tabPanel("PPI", 
                                           withSpinner(plotOutput('network', width = "100%", height = "400px"), type = 6, color = "#0277bd"),
                                           downloadButton("download_network", "Download PPI(PNG)", class = "btn-sm btn-primary")),
                                  tabPanel("GO Enrichment", 
                                           withSpinner(plotOutput("godotplot", width = "80%", height = "400px"), type = 6, color = "#0277bd"),
                                           downloadButton("download_godotplot", "Download GO(PNG)", class = "btn-sm btn-primary"))
                                )
                            ),
                            # Table result
                            h3("📋 List of Interacting Proteins"),
                            br(),
                            column(7,
                                   withSpinner(DT::dataTableOutput("analysis_table"), type = 6, color = "#0277bd"),
                                   downloadLink(outputId = "download_table_csv", icon = icon("download"), label = "Download (CSV)")),
                            column(5,
                                   withSpinner(DT::dataTableOutput("analysis_table2"), type = 6, color = "#0277bd"),
                                   downloadLink(outputId = "download_table_csv2", icon = icon("download"), label = "Download (CSV)"))
                     )
              ),
###------ Download  ------###
             tabPanel(title = "Download", value = "download", icon = icon("download"),
                     br(), hr(), 
                     titlePanel(h2("Download interaction by Protein", style = "color: #0277bd; text-align: left;")),
                     ##-- Protein download ----
                     wellPanel(style = "background: #E1F5FE",
                               fluidRow(
                                 column(8, align = "left", p("•  Download lncRNA-protein interaction information for ",strong("all")," proteins.", style = "font-size: 20px;")),
                                 column(4, align = "right", downloadButton("downloas_all_data","Download"))
                               ),
                               hr(),
                               fluidRow(
                                 column(8, align = "left", p("•  Download lncRNA-protein interaction information for ",strong("TFs")," only.", style = "font-size: 20px;")),
                                 column(4, align = "right", downloadButton("downloas_TF_data","Download"))
                               ),
                               hr(),
                               fluidRow(
                                 column(8, align = "left", p("•  Download lncRNA-protein interaction information for typical ",strong("RBPs")," only.", style = "font-size: 20px;")),
                                 column(4, align = "right", downloadButton("downloas_RBP_data","Download"))
                               )
                     )
              ),
###------ About  ------###
              tabPanel(title = "About", value = "about", icon = icon("ghost"),
                     br(), hr(),              
                     # Tutorial
                     titlePanel(h2("Tutorial", style = "color: #0277bd; text-align: left;")),
                     wellPanel(style = "background: #E1F5FE",
                               p("Welcome to this step-by-step tutorial on using HuRInterDB for retrieving and analyzing data.", style = "margin-top: 12px; font-size: 18px;"),
                               p(strong("How to search your intersted gene?"), style = "margin-top: 12px; font-size: 18px;"),
                               p("Users can perform data searches from either the Home Page or the Search Page. On the Home Page, 
                                  a prominent search bar allows quick queries by entering keywords.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("Alternatively, the Search Page provides advanced options for a more refined search. 
                                  Once a search is submitted, the system processes the request and returns matching results in a clear, organized table format.
                                  The search box at the top of the table supports fuzzy text matching across all columns.
                                  The table also supports interactive features like sorting by column and pagination. 
                                  The nine dropdown menus below can be used to filter the content of each field precisely.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("The Search page also supports exporting the current query results. All files are exported in CSV format and can be opened with Excel.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               br(),
                               p(strong("How to use RPI on-line analysis module?"), style = "margin-top: 12px; font-size: 18px;"),
                               p("The RPIs Analysis section provides an online tool for analyzing RNA-protein interactions.
                                 Enter the lncRNA name and click the 'Continue' button,it will start the analysis process.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("The analysis results consisting of four visualizations and one detailed results table to support in-depth interpretation.The four generated images provide intuitive insights into the data:", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("•    Word Cloud: Highlights the most frequent terms or keywords from the dataset, with font size indicating term prominence.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("•    Protein Binding Diagram: Illustrates the predicted or known binding interactions between proteins and ligands or other molecules.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("•    PPI Network displays the interaction network among proteins, showing functional relationships and key hub proteins.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("•    GO Enrichment Dot Plot: Visualizes the Gene Ontology (GO) enrichment results, with dots representing biological processes, molecular functions, or cellular components—positioned and colored by significance and enrichment score.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;"),
                               p("In addition to the visual outputs, a detailed results table is provided, listing all enriched terms, associated genes, p-values, enrichment scores, and other relevant metrics. 
                                  This table can be easily downloaded in CSV or Excel format for further analysis or reporting.
                                  Together, these outputs enable users to quickly interpret complex biological data and export detailed findings for downstream applications.", style = "margin-top: 12px; margin-left: 15px; font-size: 18px;")
                     ),
                     # Help
                     titlePanel(h2("Help", style = "color: #0277bd; text-align: left;")),
                     wellPanel(style = "background: #E1F5FE",
                               p("•  If you find that the search box on the home page isn't working, try clicking the delete icon next to the search icon, then re-enter the name of the lncRNA or protein you wish to search for.", style = "margin-top: 12px; font-size: 18px;")
                     ),
                     # Contact Us
                     titlePanel(h2("Contact Us", style = "color: #0277bd; text-align: left;")),
                     wellPanel(style = "background: #E1F5FE",
                               p("Welcome researchers from all over the world to provide valuable advice.Please feel free to contact us if you have any questions, comments or suggestions.", style = "margin-top: 12px; font-size: 18px;"),
                               p(strong("Jian Yan"),": jian.yan@nwu.edu.cn", style = "margin-top: 12px; font-size: 18px;"),
                               p(strong("Wenju Sun"),": wenju.sun@nwu.edu.cn", style = "margin-top: 12px; font-size: 18px;"),
                               p(strong("Qianwen Xie"),": xieqianwen@stumail.nwu.edu.cn", style = "margin-top: 12px; font-size: 18px;")
                     )
              )
    ),
    ##-- Footer ----
    div(class = "footer",
        includeHTML("www/footer.html")
    )
  )
)

########################### Server ###########################
server <- shinyServer(function(input, output, session){
          waiter_hide()
    
##---- HOME ----
    current_results <- reactiveVal(NULL)
    ##-- RNA search  ----
    observeEvent(input$RNA_search, {
        if(input$RNA_search != "") {
            query <- input$RNA_search
            filtered <- filter(data, lncRNA_Name == query)
            filtered <- prepare_data_with_links(filtered)
            current_results(filtered)
            updateTabsetPanel(session = session, inputId = "navbar", selected = "searchs")
        }
    })
    ##-- Protein search  ----
    observeEvent(input$protein_search, {
        if(input$protein_search != "") {
            query <- input$protein_search
            filtered <- filter(data, Protein_name == query)
            filtered <- prepare_data_with_links(filtered)
            current_results(filtered)
            updateTabsetPanel(session = session, inputId = "navbar", selected = "searchs")
        }
    })

##---- SEARCH ----
    # === Search Page: Apply filters ===
    observeEvent(input$filter_gene, {
        if(input$filter_gene != "") {
            related_proteins <- data %>%
                filter(lncRNA_Name == input$filter_gene) %>%
                pull(Protein_name) %>%
                unique() %>%
                sort()
            
            related_proteins <- c("", related_proteins)
            
            updateSelectInput(session, "filter_protein", 
                            choices = related_proteins,
                            selected = "")
        } else {
            updateSelectInput(session, "filter_protein", 
                            choices = protein_list,
                            selected = "")
        }
    })

    observeEvent(input$filter_gene, {
        if(input$filter_gene != "") {
            related_cells <- data %>%
                filter(lncRNA_Name == input$filter_gene) %>%
                pull(Cell_Line) %>%
                unique() %>%
                sort()
            
            related_cells <- c("", related_cells)
            
            updateSelectInput(session, "filter_cellline", 
                            choices = related_cells,
                            selected = "")
        } else {
            updateSelectInput(session, "filter_cellline", 
                            choices = cell_list,
                            selected = "")
        }
    })
    
    search_data <- eventReactive(input$apply_filter, {
        req(data)
        
        empty_result <- function() {
            data.frame(
                lncRNA_Name = character(),
                Protein_name = character(),
                Cell_Line = character(),
                Method = character(),
                Data = character(),
                Score = numeric(),
                Entry = character(),
                KEGG = character()
            )
        }

        df <- data
        
        gene_f <- trimws(input$filter_gene)
        prot_f <- trimws(input$filter_protein)
        cell_f <- trimws(input$filter_cellline)
        meth_f <- input$filter_method

        if (gene_f != "") {
            df <- filter(df, lncRNA_Name == gene_f)
            if (nrow(df) == 0) {
                showNotification(paste("Not found lncRNA:", gene_f), type = "warning", duration = 3)
                return(empty_result())
            }
        }
        
        if (prot_f != "") {
            df <- filter(df, Protein_name == prot_f)
            if (nrow(df) == 0) {
                showNotification(paste("Not found Protein:", prot_f), type = "warning", duration = 3)
                return(empty_result())
            }
        }

        if (cell_f != "") {
            df <- filter(df, Cell_Line == cell_f)
            if (nrow(df) == 0) {
                showNotification(paste("Not found Cell Line:", cell_f), type = "warning", duration = 3)
                return(empty_result())
            }
        }
        
        if (meth_f != "" && meth_f != "All") {
            df <- df[df$Method == meth_f, ]
            if (nrow(df) == 0) {
                showNotification(paste("Not found Method:", meth_f), type = "warning", duration = 3)
                return(empty_result())
            }
        }
        
        if (nrow(df) == 0) {
            showNotification("No matching data found.", type = "warning", duration = 3)
            return(empty_result())
        }
                
        df <- prepare_data_with_links(df)
        df
    }, ignoreNULL = FALSE)

    observeEvent(input$apply_filter, {
        current_results(search_data())
    })

    # Render table in Search page
    output$result_table <- DT::renderDataTable({
        res <- current_results()
        
        if (is.null(res) || nrow(res) == 0) {
            return(
                datatable(
                    data.frame(NOTE = "No matching results were found. Please adjust your search criteria."),
                    options = list(dom = 't', pageLength = 1),
                    rownames = FALSE
                )
            )
        }
        
        required_cols <- c("lncRNA_Name_link", "lncRNA_RNALocate_link", "Protein_name", 
                           "Protein_Domains_link", "AlphaFoldDB_link", 
                           "KEGG_link", "Cell_Line", "Method", "Score", "Data_link")
        
        existing_cols <- required_cols[required_cols %in% colnames(res)]
        
        if(length(existing_cols) == 0) {
            return(datatable(data.frame(Message = "No data available"), rownames = FALSE))
        }
        
        display_df <- unique(res[, existing_cols, drop = FALSE])
        display_df <- subset(display_df, !is.na(Method))
        
        if(nrow(display_df) == 0) {
            return(datatable(data.frame(NOTE = "No matching data found."), rownames = FALSE))
        }

        if ("Score" %in% colnames(display_df)) {
            display_df$Score <- round(display_df$Score, 2)
        }
        
        colnames(display_df) <- c("lncRNA", "lncRNA Localization", "Protein", "Protein Domains", 
                                  "Protein Structure", "Protein KEGG", "Cell Line", "Method", "Score", "Data")

        datatable(display_df, escape = FALSE,
                options = list(pageLength = 10, lengthMenu = c(10, 25, 50)),
                rownames = FALSE, filter = "top"
        )
    })
    
    # Download filtered data
    output$download_data <- downloadHandler(
        filename = function() "search_results.csv",
        content = function(file) {
            res <- current_results()
            if (!is.null(res) && nrow(res) > 0) {
                cols_to_export <- c("lncRNA_Name", "Protein_name", "Entry", "KEGG", 
                                    "Cell_Line", "Method", "Data")
                existing_cols <- cols_to_export[cols_to_export %in% colnames(res)]
                write.csv(res[existing_cols], file, row.names = FALSE)
            } else {
                write.csv(data.frame(), file, row.names = FALSE)
            }
        }
    )

    # ==================== Listen for protein click events ====================
    selected_protein_network <- reactiveVal(NULL)
    observeEvent(input$selected_protein_network, {
        selected_protein_network(input$selected_protein_network)
        showNotification(paste("Selected protein:", input$selected_protein_network), 
                        type = "message", duration = 3)
    })
    
    output$protein_network_visible <- reactive({
        !is.null(selected_protein_network()) && selected_protein_network() != ""
    })
    outputOptions(output, "protein_network_visible", suspendWhenHidden = FALSE)
    
    # 网络统计信息
    network_stats <- reactiveValues(
        original_count = 0,
        displayed_count = 0,
        sampled = FALSE,
        protein_name = ""
    )

    # 生成蛋白质关联的 lncRNA 网络图
    protein_lncRNA_network <- reactive({
        req(selected_protein_network())
        
        protein_name <- selected_protein_network()
        protein_data <- data %>% filter(Protein_name == protein_name)
        
        if(nrow(protein_data) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = paste("No lncRNA interactions found for protein:", protein_name),
                           size = 6, color = "red") +
                   theme_void())
        }
        
        lncRNAs <- unique(protein_data$lncRNA_Name)
        lncRNAs <- lncRNAs[!is.na(lncRNAs) & nchar(lncRNAs) > 0]
        original_count <- length(lncRNAs)
        
        if(length(lncRNAs) > 300) {
            set.seed(123)
            lncRNAs <- sample(lncRNAs, 300)
            sampled_flag <- TRUE
        } else {
            sampled_flag <- FALSE
        }
        
        if(length(lncRNAs) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                           label = "No lncRNA interactions found.",
                           size = 6, color = "red") +
                   theme_void())
        }
        
        network_stats$original_count <- original_count
        network_stats$displayed_count <- length(lncRNAs)
        network_stats$sampled <- sampled_flag
        network_stats$protein_name <- protein_name
        
        edges <- data.frame(
            from = rep(protein_name, length(lncRNAs)),
            to = lncRNAs,
            stringsAsFactors = FALSE
        )
        
        nodes <- data.frame(
            name = c(protein_name, lncRNAs),
            type = c("Protein", rep("lncRNA", length(lncRNAs))),
            stringsAsFactors = FALSE
        )
        
        g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
        set.seed(123)
        layout <- layout_with_fr(g)
        coords <- data.frame(
            x = layout[, 1],
            y = layout[, 2],
            name = V(g)$name,
            type = V(g)$type,
            stringsAsFactors = FALSE
        )
        
        edge_coords <- data.frame()
        for(i in 1:nrow(edges)) {
            from_idx <- which(coords$name == edges$from[i])
            to_idx <- which(coords$name == edges$to[i])
            
            if(length(from_idx) > 0 && length(to_idx) > 0) {
                edge_coords <- rbind(edge_coords, data.frame(
                    x = coords$x[from_idx],
                    y = coords$y[from_idx],
                    xend = coords$x[to_idx],
                    yend = coords$y[to_idx]
                ))
            }
        }
        
        coords$color <- ifelse(coords$type == "Protein", "#E64B35", "#4DBBD5")
        coords$size <- ifelse(coords$type == "Protein", 15, 6)
        coords$label_color <- ifelse(coords$type == "Protein", "black", "gray30")
        
        p <- ggplot() +
            geom_segment(data = edge_coords, 
                        aes(x = x, y = y, xend = xend, yend = yend),
                        color = "gray70", linewidth = 0.5, alpha = 0.6) +
            geom_point(data = coords, 
                      aes(x = x, y = y, color = color, size = size),
                      alpha = 0.9) +
            geom_text(data = coords,
                     aes(x = x, y = y, label = name, color = label_color),
                     family = "DejaVu Sans",
                     size = 3.5,
                     vjust = -0.8,
                     hjust = 0.5) +
            scale_color_identity() +
            scale_size_identity() +
            labs(title = paste("LncRNA Interaction Network for Protein:", protein_name),
                 subtitle = paste0("Total interacting lncRNAs: ", original_count,
                                  ifelse(sampled_flag, paste0(" (Displayed: ", length(lncRNAs), " - random sample)"), ""))) +
            theme_void() +
            theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(family = "DejaVu Sans", face = "bold", 
                                          size = 14, hjust = 0.5),
                  plot.subtitle = element_text(family = "DejaVu Sans", size = 10, 
                                             color = "gray50", hjust = 0.5))
        
        return(p)
    })
            
    # 渲染网络图
    output$protein_lncRNA_network <- renderPlot({
        req(protein_lncRNA_network())
        protein_lncRNA_network()
    })
    
    # 下载网络图
    output$download_protein_network <- downloadHandler(
        filename = function() {
            paste0("protein_network_", selected_protein_network(), "_", Sys.Date(), ".png")
        },
        content = function(file) {
            req(protein_lncRNA_network())
            ggsave(file, plot = protein_lncRNA_network(), 
                   device = "png", width = 14, height = 12, dpi = 300)
        }
    )



##---- RPI Analysis ----
    analysis_result <- eventReactive(input$analyze_btn, {
        tryCatch({
            rna_name <- trimws(input$lncrna_input)
            cell_f <- trimws(input$cellline_input)
            meth_f <- input$method_input
            
            req(data)
            
            empty_df <- data.frame(
                lncRNA_Name = character(),
                Protein_name = character(),
                Cell_Line = character(),
                Method = character(),
                Data = character(),
                Entry = character(),
                KEGG = character()
            )
            
            df <- data
            
            if (rna_name != "") {
                if (rna_name %in% df$lncRNA_Name) {
                    df <- filter(df, lncRNA_Name == rna_name)
                } else {
                    showNotification(paste("Not found lncRNA:", rna_name), type = "warning")
                    return(NULL)
                }
            } else {
                showNotification("Please select an lncRNA", type = "warning")
                return(NULL)
            }
            
            if (cell_f != "") {
                if (cell_f %in% df$Cell_Line) {
                    df <- filter(df, Cell_Line == cell_f)
                } else {
                    showNotification(paste("Not found Cell Line:", cell_f), type = "warning")
                    return(NULL)
                }
            }
            
            if (meth_f != "" && meth_f != "All") {
                if (meth_f %in% df$Method) {
                    df <- df[df$Method == meth_f, ]
                } else {
                    showNotification(paste("Not found Method:", meth_f), type = "warning")
                    return(NULL)
                }
            }
            
            if (nrow(df) == 0) {
                showNotification("No matching data found", type = "warning")
                return(NULL)
            }
            
            cols_to_select <- c("lncRNA_Name", "Protein_name", "Cell_Line", "Method", "Data", "Entry", "KEGG")
            cols_exist <- cols_to_select[cols_to_select %in% names(df)]
            res <- df[, cols_exist, drop = FALSE]
            
            missing_cols <- setdiff(cols_to_select, names(res))
            for (col in missing_cols) {
                res[[col]] <- NA
            }
            
            lncRNA_info <- lncRNA_bed_data %>% filter(lncRNA_name == rna_name)
            
            if(nrow(lncRNA_info) == 0) {
                showNotification(paste("No lncRNA location information found:", rna_name), type = "warning")
                return(list(
                    res = res,
                    prot_bindings = data.frame(),
                    exons = data.frame(),
                    pos_min = NA,
                    pos_max = NA,
                    timestamp = Sys.time()
                ))
            }
            
            exons <- lncRNA_info %>% filter(region == "exon")
            chr_target <- as.character(lncRNA_info$chr[1])
            pos_min <- min(lncRNA_info$start)
            pos_max <- max(lncRNA_info$end)
            
            prot_bindings <- protein_binding_data %>%
                filter(chr == chr_target, between(position, pos_min, pos_max))
            
            list(
                res = res,
                prot_bindings = prot_bindings,
                exons = exons,
                pos_min = pos_min,
                pos_max = pos_max,
                timestamp = Sys.time()
            )
            
        }, error = function(e) {
            print(paste("Error:", e$message))
            showNotification(paste("Error:", e$message), type = "error")
            return(NULL)
        })
    })

    # Generate wordcloud from Protein_name
    wordcloud_plot <- reactive({
        req(analysis_result())
        
        res <- analysis_result()$res
        if(is.null(res) || nrow(res) == 0) {
            return(wordcloud2(data.frame(word = "No data", freq = 1), 
                              size = 1, color = "red"))
        }
        
        protein_freq <- as.data.frame(table(res$Protein_name))
        colnames(protein_freq) <- c("word", "freq")
        
        if (nrow(protein_freq) > 100) {
            protein_freq <- protein_freq[order(-protein_freq$freq), ][1:100, ]
        }
        
        wordcloud2(protein_freq, size = 0.8,
                   color = "random-light", backgroundColor = "white")
    })

    output$wordcloud <- renderWordcloud2({
        req(wordcloud_plot())
        wordcloud_plot()
    })

    output$download_wordcloud <- downloadHandler(
        filename = function() {
            paste0("wordcloud_", input$lncrna_input, ".html")
        },
        content = function(file) {
            req(wordcloud_plot())
            saveWidget(wordcloud_plot(), file = file, selfcontained = TRUE)
        }
    )

    # Generate lolliplot
    RBP_plot <- reactive({
        req(input$lncrna_input)
        req(analysis_result())
        
        res <- analysis_result()$res
        if(is.null(res) || nrow(res) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "Non-interacting RBP.", size = 6, color = "red") +
                   theme_void())
        }
        
        prot_bindings <- analysis_result()$prot_bindings
        exons <- analysis_result()$exons
        pos_min <- analysis_result()$pos_min
        pos_max <- analysis_result()$pos_max
        
        if(nrow(prot_bindings) == 0 || is.na(pos_min) || is.na(pos_max)) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No RBP binding data available.", size = 6, color = "red") +
                   theme_void())
        }
        
        rna_name <- trimws(input$lncrna_input)
        
        p <- ggplot() +
            geom_hline(yintercept = 0.1, linewidth = 1.5, color = "black", linetype = "solid") +
            geom_rect(data = exons, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.2), 
                      fill = rep("#7EB7DC", nrow(exons))) +
            geom_segment(data = prot_bindings, aes(x = position, xend = position, y = 0.2, yend = 1),
                        linewidth = 0.3, colour = "black") +
            geom_point(data = prot_bindings, aes(x = position, y = 1), 
                       size = 3, alpha = 0.7, fill = "#E64B35", shape = 21) +
            geom_text(data = prot_bindings, aes(x = position, y = 1.08, label = protein_name),
                      family = "DejaVu Sans", size = 6, angle = 90, hjust = 0, color = "black") +
            xlim(pos_min, pos_max) +
            ylim(-0.2, 1.8) +
            theme_void() +
            labs(title = rna_name) +
            theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(hjust = 0.5, size = 16, family = "DejaVu Sans", face = "bold"))
        
        return(p)
    })

    output$lolliplot <- renderPlot({
        req(RBP_plot())
        RBP_plot()
    })

    output$download_lolliplot <- downloadHandler(
        filename = function() {
            paste0("ENCODE_RBP_binding_plot_", input$lncrna_input, ".png")
        },
        content = function(file) {
            req(RBP_plot())
            ggsave(file, plot = RBP_plot(), 
                   device = "png", width = 10, height = 6, dpi = 600)
        }
    )

    # Generate PPI from Protein_name
    PPI_plot <- reactive({
        req(analysis_result())
        res <- analysis_result()$res
        if(is.null(res) || nrow(res) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "Non-interacting proteins.", size = 6, color = "red") +
                   theme_void())
        }
        
        proteins_freq <- as.data.frame(table(res$Protein_name))
        colnames(proteins_freq) <- c("protein_symbol","freq")
        proteins_freq <- proteins_freq[which(proteins_freq$freq > 0), ]
        proteins_freq$protein_symbol <- as.character(proteins_freq$protein_symbol)

        if(nrow(proteins_freq) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "Non-interacting proteins.", size = 6, color = "red") +
                   theme_void())
        }

        g <- getPPI(proteins_freq$protein_symbol, taxID="9606")

        if(is.null(g) || length(V(g)) == 0 || length(E(g)) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "Non-interacting proteins.", size = 6, color = "red") +
                   theme_void())
        }
        
        # ==================== Calculate various centrality metrics ====================
        degree_centrality <- degree(g, mode = "all")
        betweenness_centrality <- betweenness(g, directed = FALSE)
        closeness_centrality <- closeness(g, mode = "all")
        eigen_centrality <- eigen_centrality(g, directed = FALSE)$vector
        degree_norm <- degree_centrality / max(degree_centrality)
        betweenness_norm <- betweenness_centrality / max(betweenness_centrality)
        composite_score <- degree_norm * 0.4 + betweenness_norm * 0.3 + 
                        (closeness_centrality / max(closeness_centrality)) * 0.3
        hub_threshold <- quantile(composite_score, 0.8)
        is_hub <- composite_score >= hub_threshold
        hub_level <- ifelse(composite_score >= quantile(composite_score, 0.95), "Level 3 (Key Hub)",
                    ifelse(composite_score >= quantile(composite_score, 0.8), "Level 2 (Moderate Hub)",
                        "Level 1 (Peripheral)"))
        node_attributes <- data.frame(
            name = names(degree_centrality),
            degree = degree_centrality,
            betweenness = betweenness_centrality,
            closeness = closeness_centrality,
            composite_score = composite_score,
            is_hub = is_hub,
            hub_level = hub_level,
            freq = proteins_freq$freq[match(names(degree_centrality), proteins_freq$protein_symbol)]
        )
        
        node_attributes$node_color <- ifelse(node_attributes$is_hub, "#E64B35", "#4DBBD5")
        node_attributes$node_size <- ifelse(node_attributes$is_hub, 10, 5)
        node_attributes$label_color <- ifelse(node_attributes$is_hub, "black", "gray40")
        hub_proteins <- node_attributes[node_attributes$is_hub, ]
        hub_proteins <- hub_proteins[order(hub_proteins$composite_score, decreasing = TRUE), ]

        # ==================== Plotting a PPI network ====================
        p <- ggplot(g, layout='circle') %<+% node_attributes + 
            geom_edge(aes(alpha = 0.3), color = "gray70") + 
            geom_point(aes(color = node_color, size = node_size), alpha = 0.9) + 
            shadowtext::geom_shadowtext(aes(label = name, color = label_color), 
                                        family = "DejaVu Sans", 
                                        bg.color = "white",
                                        size = 3.5) +
            scale_color_identity() +
            scale_size_identity() +
            annotate("text", x = -Inf, y = Inf, 
                    label = paste0("★ Hub Protein (N = ", sum(is_hub), ")"), 
                    hjust = -0.1, vjust = 1.5, 
                    color = "#E64B35", fontface = "bold", size = 5,
                    family = "DejaVu Sans") +
            theme_void() +
            theme(text = element_text(family = "DejaVu Sans", size = 10),
                plot.background = element_rect(fill = "white", color = NA),
                plot.title = element_text(family = "DejaVu Sans", face = "bold", size = 13, hjust = 0.5),
                plot.subtitle = element_text(family = "DejaVu Sans", size = 9, color = "gray50", hjust = 0.5))
        
        return(p)
    })

    output$network <- renderPlot({
        req(PPI_plot())
        PPI_plot()
    })
    
    output$download_network <- downloadHandler(
        filename = function() {
            paste0("PPI_network_plot_", input$lncrna_input, ".png")
        },
        content = function(file) {
            req(PPI_plot())
            ggsave(file, plot = PPI_plot(), 
                   device = "png", width = 12, height = 10, dpi = 600)
        }
    )

    # Generate GO from Protein_name
    GO_plot <- reactive({
        req(analysis_result())
        res <- analysis_result()$res
        if(is.null(res) || nrow(res) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No enriched pathways", size = 6, color = "red") +
                   theme_void())
        }
        
        proteins_list <- unique(res$Protein_name)
        proteins_list <- proteins_list[nchar(proteins_list) > 0]
        
        if(length(proteins_list) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No enriched pathways", size = 6, color = "red") +
                   theme_void())
        }

        go_result <- enrichGO(gene = proteins_list,
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.99,
                              qvalueCutoff = 0.99)
        
        if(is.null(go_result) || nrow(go_result) == 0) {
            return(ggplot() + 
                   annotate("text", x = 0.5, y = 0.5, 
                            label = "No enriched pathways", size = 6, color = "red") +
                   theme_void())
        }

        p <- dotplot(go_result, showCategory = 15, 
                    color = "pvalue",
                    label_format = 80,
                    font.size = 11,
                    title = "GO Pathway Enrichment") +
             theme(text = element_text(size = 12),
                   axis.text.y = element_text(family = "DejaVu Sans", size = 10),
                   axis.text.x = element_text(family = "DejaVu Sans", size = 9),
                   axis.title = element_text(family = "DejaVu Sans", size = 12),
                   legend.text = element_text(family = "DejaVu Sans"),
                   legend.title = element_text(family = "DejaVu Sans"),
                   plot.background = element_rect(fill = "white", color = NA),
                   plot.title = element_text(family = "DejaVu Sans", face = "bold"))
        return(p)
    })

    output$godotplot <- renderPlot({
        req(GO_plot())
        GO_plot()
    })

    output$download_godotplot <- downloadHandler(
        filename = function() {
            paste0("GO_enrichment_Plot_", input$lncrna_input, ".png")
        },
        content = function(file) {
            req(GO_plot())
            ggsave(file, plot = GO_plot(), 
                   device = "png", width = 12, height = 8, dpi = 300)
        }
    )

    # Display detailed table
    output$analysis_table <- DT::renderDataTable({
        req(analysis_result())
        res <- analysis_result()$res
        if (is.null(res) || nrow(res) == 0) return(NULL)
        
        datatable(res[, c("lncRNA_Name", "Protein_name", "Cell_Line", "Method", "Data")], 
                  escape = FALSE, options = list(pageLength = 10), rownames = FALSE)
    })
    
    ##---- Download table ----
    output$download_table_csv <- downloadHandler(
        filename = function() {
            paste0(input$lncrna_input, "_binding_protein.csv")
        },
        content = function(file) {
            res <- analysis_result()$res
            if(!is.null(res) && nrow(res) > 0) {
                write.csv(res[, c("lncRNA_Name", "Protein_name", "Cell_Line", "Method", "Data")], 
                         file, row.names = FALSE)
            } else {
                write.csv(data.frame(), file, row.names = FALSE)
            }
        }
    )

    # Display RBP binding table
    output$analysis_table2 <- DT::renderDataTable({
        req(analysis_result())
        prot_bindings <- analysis_result()$prot_bindings
        if (is.null(prot_bindings) || nrow(prot_bindings) == 0) return(NULL)
        
        datatable(prot_bindings, escape = FALSE, options = list(pageLength = 10), rownames = FALSE)
    })
    
    ##---- Download table ----
    output$download_table_csv2 <- downloadHandler(
        filename = function() {
            paste0(input$lncrna_input, "_RBP_binding_data.csv")
        },
        content = function(file) {
            prot_bindings <- analysis_result()$prot_bindings
            if(!is.null(prot_bindings) && nrow(prot_bindings) > 0) {
                write.csv(prot_bindings, file, row.names = FALSE)
            } else {
                write.csv(data.frame(), file, row.names = FALSE)
            }
        }
    )

##---- DOWNLOAD ----
    ##---- Download ALL ----
    output$downloas_all_data <- downloadHandler(
        filename = function() {"HuRInterDB_all_RIPs_data.csv"},
        content = function(file) {
            cols <- c("lncRNA_Name","Protein_name","Cell_Line","Method")
            existing_cols <- cols[cols %in% colnames(data)]
            write.csv(data[, existing_cols], file, row.names = FALSE)
        }
    )
    
    ##---- Download TF ----
    output$downloas_TF_data <- downloadHandler(
        filename = function() {"HuRInterDB_TF_only_RIPs_data.csv"},
        content = function(file) {
            if("TF" %in% colnames(data)) {
                TF_data <- data %>%
                    filter(TF == "TF") %>%
                    select(any_of(c("lncRNA_Name", "Protein_name", "Cell_Line", "Method")))
                write.csv(TF_data, file, row.names = FALSE)
            } else {
                write.csv(data.frame(), file, row.names = FALSE)
            }
        }
    )
    
    ##---- Download RBP ----
    output$downloas_RBP_data <- downloadHandler(
        filename = function() {"HuRInterDB_RBP_only_RIPs_data.csv"},
        content = function(file) {
            if("RBP" %in% colnames(data)) {
                RBP_data <- data %>%
                    filter(RBP == "RBP") %>%
                    select(any_of(c("lncRNA_Name", "Protein_name", "Cell_Line", "Method")))
                write.csv(RBP_data, file, row.names = FALSE)
            } else {
                write.csv(data.frame(), file, row.names = FALSE)
            }
        }
    )
})

shinyApp(ui, server)