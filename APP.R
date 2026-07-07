library(shiny)
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
library(arrow)

# Global plot configuration
theme_set(theme_bw(base_family = "DejaVu Sans"))
options(ggplot2.device = "png")
# Disable multi-session to simplify error debugging
# plan(multisession, workers = 2)
# File path constants
PATH_DATA_PARQUET <- "data/data.parquet"
PATH_BED_PARQUET <- "data/lncRNA_bed_data.parquet"
PATH_PROT_BIND_PARQUET <- "data/protein_binding_data.partition/"
# Cache environment for large protein binding dataset only
data_cache_env <- new.env(hash = TRUE, parent = emptyenv(), size = 3L)
# Load lightweight lncRNA coordinate data (small volume, no startup block)
lncRNA_bed_data <- read_parquet(PATH_BED_PARQUET)

fetch_filtered_data <- function(filter_rna = NULL, filter_prot = NULL, filter_cell = NULL, filter_method = NULL) {
  ds <- open_dataset(PATH_DATA_PARQUET)
  df <- ds %>% collect()

  if (!is.null(filter_rna) && trimws(filter_rna) != "") {
    kw <- trimws(filter_rna)
    mask <- grepl(paste0("^", kw, "$"), df$lncRNA_Name, ignore.case = TRUE)
    df <- df[mask, , drop = FALSE]
  }
  if (!is.null(filter_prot) && trimws(filter_prot) != "") {
    kw <- trimws(filter_prot)
    mask <- grepl(paste0("^", kw, "$"), df$Protein_name, ignore.case = TRUE)
    df <- df[mask, , drop = FALSE]
  }
  if (!is.null(filter_cell) && trimws(filter_cell) != "") {
    kw <- trimws(filter_cell)
    mask <- grepl(paste0("^", kw, "$"), df$Cell_Line, ignore.case = TRUE)
    df <- df[mask, , drop = FALSE]
  }
  if (!is.null(filter_method) && trimws(filter_method) != "") {
    df <- df[df$Method == filter_method, , drop = FALSE]
  }
  return(df)
}

# Fuzzy match cell line names by input keyword
match_cell <- function(keyword, limit = 30) {
  kw <- trimws(keyword)
  if (nchar(kw) < 2) return(character(0))
  ds <- open_dataset(PATH_DATA_PARQUET)
  dt <- ds %>% filter(str_detect(tolower(Cell_Line), tolower(kw))) %>% select(Cell_Line) %>% collect()
  res <- dt %>% distinct() %>% slice_head(n = limit) %>% pull(Cell_Line)
  return(res)
}

# Pre-generate hyperlink HTML columns for data table rendering
prepare_data_with_links <- function(df) {
  if (!"Entry" %in% colnames(df)) df$Entry <- ""
  if (!"KEGG" %in% colnames(df)) df$KEGG <- ""

  
  df$lncRNA_Name_link <- with(df, sprintf(
    '<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank">%s</a>',
    lncRNA_Name, lncRNA_Name
  ))
  
  # Fix protein link: href="#" removed, add onclick preventDefault
  df$Protein_name_link <- with(df, sprintf(
    '<a class="protein-network-link" data-protein="%s" style="color: #0066cc; text-decoration: underline; cursor: pointer;" onclick="event.preventDefault();">%s</a>',
    Protein_name, Protein_name
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

  df$KEGG_link <- with(df, ifelse(KEGG != "", 
    sprintf('<a href="https://www.kegg.jp/entry/%s" target="_blank">%s</a>', KEGG, KEGG),
    ""
  ))
  
  generate_single_link <- function(data_id) {
    if (data_id == "" || is.na(data_id)) return("")
    if (grepl("^GSE", data_id)) {
      return(sprintf('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" target="_blank">%s</a>', data_id, data_id))
    } else if (grepl("^ENCSR", data_id)) {
      return(sprintf('<a href="https://www.encodeproject.org/%s/" target="_blank">%s', data_id, data_id))
    } else if (grepl("^[0-9]+$", data_id)) {
      return(sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s" target="_blank">%s</a>', data_id, data_id))
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
  return(df)
}



# Lazy load chromosome partitioned protein binding data with cache
get_protein_binding_data <- function(chr_filter = NULL, pos_min = NULL, pos_max = NULL) {
  if (exists("prot_bind_full", envir = data_cache_env)) {
    full_dt <- get("prot_bind_full", envir = data_cache_env)
  } else {
    ds <- open_dataset(PATH_PROT_BIND_PARQUET)
    full_dt <- ds %>% collect()
    assign("prot_bind_full", full_dt, envir = data_cache_env)
  }
  if (!is.null(chr_filter) && !is.null(pos_min) && !is.null(pos_max)) {
    full_dt <- full_dt %>% filter(chr == chr_filter, between(position, pos_min, pos_max))
  }
  return(full_dt)
}

# UI Definition: Replace undefined autocompleteInput with selectize textInput
ui <- shinyUI(
  fluidPage(
    tags$head(
      tags$link(rel = "shortcut icon", href = "img/logo.ico"),
      tags$link(rel="stylesheet", type = "text/css",
                href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"),
      tags$link(rel="stylesheet", type = "text/css",
                href = "https://fonts.googleapis.com/css?family=Open+Sans|Source+Sans+Pro")
    ),
    list(tags$head(HTML('<link rel="icon", href="img/logo.png",type="image/png" />'))),
    div(style="padding: 1px 0px; width: '50%'",
        titlePanel(title ="", windowTitle = "HuRInterDB")
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
    navbarPage(
      title = div(img(src = "img/logo.png", height = "50px"), style = "padding-left:40px;"),
      id = "navbar",
      selected = "home",
      theme = "styles.css",
      fluid = TRUE,
      # Home Tab
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
      # Search Tab: Replace autocompleteInput with textInput + selectize
      tabPanel(title = "Search", value = "searchs", icon = icon("magnifying-glass"),
               br(), hr(),
               column(width = 10,style = "padding-top: 0px;",
                      column(2, textInput("filter_g", label = "lncRNA", value = "")),
                      column(2, textInput("filter_p", label = "Protein", value = "")),
                      column(2, textInput("filter_c", label = "Cell Line", value = "")),
                      column(2, selectInput("filter_m", label = "Method",
                                            c("","RNA Pulldown","ChIRP-MS","RAP-MS","HyPR-MS","HPLC-MS","CARPID","SILAC-MS","TREX","RIP-seq","eCLIP-seq","CLIP-seq","PAR-CLIP","HITS-CLIP","LACE-seq","ARTR-seq","PRIM-seq")))
               ),
               column(width = 2, style = "padding-top:55px;",
                      actionBttn(inputId = "apply_filter", label = "Select", style = "fill", color = "success", icon = icon("check"), size = "sm")),
               column(12,
                      h3("Search results:", style = "color: #0277bd; text-align: left;"),
                      withSpinner(DTOutput("result_table"),type=6,color="#0277bd"),
                      downloadButton("download_data","Download Results"),
                      br(),br(),
                      conditionalPanel(condition = "output.protein_network_visible == true",
                                       column(width = 12,
                                       h3("LncRNA Interaction Network for Selected Protein:", style = "color: #0277bd;"),
                                       withSpinner(plotOutput("protein_lncRNA_network",height = "400px", width = "100%"),type = 6, color = "#0277bd"),
                                       downloadButton("download_protein_network","Download Network (PNG)",class="btn-sm btn-primary")))
               )
      ),

      # Analysis Tab: Replace autocompleteInput with textInput
      tabPanel(title = "RPIanalysis", value = "RPIanalysis", icon = icon("atom"),
               br(), hr(),
               tags$style(".inputBox{background:#E1F5FE;padding:10px;margin-bottom:10px}.resultBox{background:#E1F5FE;padding:10px;border-radius:10px;margin:10px 0}"),
               h2("Online lncRNA-protein Interaction Analysis", style="color:#0277bd;text-align:center"),
               column(width = 10,style="padding-top:0px;",
                      column(2, textInput("ana_rna", label = "lncRNA", value = "")),
                      column(2, textInput("ana_cell", label = "Cell Line", value = "")),
                      column(2, selectInput("method_input","Method",
                                            c("","RNA Pulldown","ChIRP-MS","RAP-MS","HyPR-MS","HPLC-MS","CARPID","SILAC-MS","TREX","RIP-seq","eCLIP-seq","CLIP-seq","PAR-CLIP","HITS-CLIP","LACE-seq","ARTR-seq","PRIM-seq")))
               ),
              column(width = 2, style = "padding-top: 55px;",
                     actionBttn(inputId = "analyze_btn", label = "Continue", style = "fill", color = "success", icon = icon("arrow-right"), size = "sm")),

               column(12, div(class="resultBox",
                              h3("Interaction Analysis Results"),
                              tabsetPanel(
                                tabPanel("Protein Wordcloud", withSpinner(wordcloud2Output("wordcloud"),type=6),downloadButton("download_wordcloud")),
                                tabPanel("RBP binding", withSpinner(plotOutput("lolliplot"),type=6),downloadButton("download_lolliplot")),
                                tabPanel("PPI Network", withSpinner(plotOutput("network"),type=6),downloadButton("download_network")),
                                tabPanel("GO Enrichment", withSpinner(plotOutput("godotplot"),type=6),downloadButton("download_godotplot"))
                              )),
               h3("List of Interacting Proteins"),br(),
               column(7,withSpinner(DTOutput("analysis_table"),type=6),downloadLink("download_table_csv","CSV")),
               column(5,withSpinner(DTOutput("analysis_table2"),type=6),downloadLink("download_table_csv2","CSV"))
               )
      ),

      # Download Tab
      tabPanel(title = "Download", value = "download", icon = icon("download"),
               br(),hr(),
               h2("Bulk Interaction Data Download", style="color:#0277bd"),
               wellPanel(style="background:#E1F5FE",
                         fluidRow(column(8,p("Download all lncRNA-protein interaction records")),column(4,downloadButton("downloas_all_data","Download"))),hr(),
                         fluidRow(column(8,p("Download TF subset interaction records")),column(4,downloadButton("downloas_TF_data","Download"))),hr(),
                         fluidRow(column(8,p("Download RBP subset interaction records")),column(4,downloadButton("downloas_RBP_data","Download")))
               )
      ),
      # About Tab
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
    div(class="footer", includeHTML("www/footer.html"))
  )
)

# Server Backend Logic
server <- shinyServer(function(input, output, session){
  current_results <- reactiveVal(NULL)
  selected_protein_network <- reactiveVal(NULL)
  network_stats <- reactiveValues(original_count=0, displayed_count=0, sampled=F, protein_name="")

  # Close loading spinner after UI render
  session$onFlushed(function(){
    waiter_hide()
  }, once = TRUE)

  # Homepage lncRNA search
  observeEvent(input$RNA_search, {
    query <- trimws(input$RNA_search)
    if (query == "") return()
    df <- fetch_filtered_data(filter_rna = query)
    df_linked <- prepare_data_with_links(df)
    if(nrow(df_linked) == 0){
      showNotification(paste("No record matched lncRNA:", query), type="warning", duration=4)
      return()
    }
    current_results(df_linked)
    updateNavbarPage(session, "navbar", selected = "searchs")
  })

  # Homepage protein search
  observeEvent(input$protein_search, {
    query <- trimws(input$protein_search)
    if (query == "") return()
    df <- fetch_filtered_data(filter_prot = query)
    df_linked <- prepare_data_with_links(df)
    if(nrow(df_linked) == 0){
      showNotification(paste("No record matched protein:", query), type="warning", duration=4)
      return()
    }
    current_results(df_linked)
    updateNavbarPage(session, "navbar", selected = "searchs")
  })

  # Search filter button
  search_data <- eventReactive(input$apply_filter, {
    gr <- trimws(input$filter_g)
    pr <- trimws(input$filter_p)
    cl <- trimws(input$filter_c)
    mt <- input$filter_m
    fetch_filtered_data(filter_rna = gr, filter_prot = pr, filter_cell = cl, filter_method = mt)
  })
  observeEvent(input$apply_filter, {
    df_raw <- search_data()
    df_linked <- prepare_data_with_links(df_raw)
    current_results(df_linked)
    if(nrow(df_linked) == 0) showNotification("No matched interaction records", type="warning")
  })

  # Render search table
  output$result_table <- renderDT({
    req(current_results())
    res <- current_results()
    if(nrow(res) == 0){
      return(datatable(data.frame(NOTE = "No matched data"), rownames = FALSE))
    }
    target_cols <- c("lncRNA_Name_link","lncRNA_RNALocate_link","Protein_name_link","Protein_Domains_link","KEGG_link","Cell_Line","Method","Score","Data_link")
    valid_cols <- target_cols[target_cols %in% colnames(res)]
    disp <- unique(res[, valid_cols, drop = FALSE])

    if ("Score" %in% colnames(disp)) {
      disp$Score <- round(disp$Score, 2)
    }

    # Dynamic JS callback for paginated table rows
    table_callback <- JS('
      table.on("click", ".protein-network-link", function(e) {
        e.preventDefault();
        var prot = $(this).data("protein");
        Shiny.setInputValue("selected_protein_network", prot, {priority: "event"});
      });
    ')

    datatable(disp,escape = FALSE,filter = "top",rownames = FALSE,colnames = c("lncRNA","lncRNA Locate","Protein","Protein Domain","Protein KEGG","Cell Line","Method","Score","Reference"),
      callback = table_callback, # Critical dynamic event binding
      options = list(pageLength = 10,lengthMenu = c(10, 25, 50),deferRender = TRUE # Speed up large dataset rendering
      )
    )
  })

  # Search result download
  output$download_data <- downloadHandler(
    "search_results.csv",
    function(file){
      res <- current_results()
      if(!is.null(res) && nrow(res)>0){
        exp <- c("lncRNA_Name","Protein_name","Entry","Cell_Line","Method","Data")
        exp <- exp[exp %in% colnames(res)]
        write.csv(res[,exp], file, row.names=F)
      }else write.csv(data.frame(), file, row.names=F)
    }
  )

  # Protein click trigger
  observeEvent(input$selected_protein_network, {
    selected_protein_network(input$selected_protein_network)
    showNotification(paste("Selected protein:", input$selected_protein_network), type = "message",duration=3)
  })
  output$protein_network_visible <- reactive(!is.null(selected_protein_network()) && selected_protein_network() != "")
  outputOptions(output, "protein_network_visible", suspendWhenHidden=F)

  # Protein-lncRNA network plot
  protein_lncRNA_network <- reactive({
    req(selected_protein_network())
    pn <- selected_protein_network()
    dt <- fetch_filtered_data(filter_prot = pn)
    if(nrow(dt)==0) return(ggplot()+annotate("text",0.5,0.5,paste("No interaction for",pn),size=6,color="red")+theme_void())
    lnc_vec <- unique(dt$lncRNA_Name[!is.na(dt) & nchar(dt)>0])
    sample_flag <- F
    if(length(lnc_vec)>300){
      set.seed(123)
      lnc_vec <- sample(lnc_vec,300)
      sample_flag <- T
    }
    edges <- data.frame(from=rep(pn,length(lnc_vec)), to=lnc_vec)
    nodes <- data.frame(name=c(pn,lnc_vec), type=c("Protein",rep("lncRNA",length(lnc_vec))))
    g <- graph_from_data_frame(edges, F, nodes)
    set.seed(123)
    lay <- layout_with_fr(g)
    coords <- data.frame(x=lay[,1],y=lay[,2],name=V(g)$name,type=V(g)$type)
    coords$color <- ifelse(coords$type=="Protein","#E64B35","#4DBBD5")
    coords$size <- ifelse(coords$type=="Protein",15,6)
    ggplot() +
      geom_segment(data=inner_join(edges,coords,by=c("from"="name")) %>% inner_join(coords,by=c("to"="name")), aes(x=x.x,y=y.x,xend=x.y,yend=y.y), color="gray70") +
      geom_point(data=coords,aes(x=x,y,color=color,size=size),alpha=0.9) +
      geom_text(data=coords,aes(x=x,y,label=name), family = "DejaVu Sans",vjust=-0.8,size=3.5) +
      scale_color_identity() + scale_size_identity() + theme_void() +
      labs(title=paste("Interaction Network of",pn)) +
      theme(plot.title = element_text(family = "DejaVu Sans"))
  })

  output$protein_lncRNA_network <- renderPlot(protein_lncRNA_network())
  output$download_protein_network <- downloadHandler(
    paste0("network_",selected_protein_network(),".png"),
    function(file) 
    ggsave(file, protein_lncRNA_network(), w=14,h=12,dpi=300)
  )

  # Main analysis pipeline
  analysis_result <- eventReactive(input$analyze_btn, {
    rna <- trimws(toupper(input$ana_rna))
    cell <- trimws(input$ana_cell)
    meth <- input$method_input
    if(rna==""){
      showNotification("Please input target lncRNA", type="warning")
      return(NULL)
    }
    tryCatch({
      df <- fetch_filtered_data(filter_rna=rna,filter_cell=cell,filter_method=meth)
      if(nrow(df)==0){
        showNotification("No matched data", type="warning")
        return(NULL)
      }
      lnc_info <- lncRNA_bed_data %>% filter(lncRNA_name == rna)
      if(nrow(lnc_info)==0) return(list(res=df,prot_bind=data.frame(),exons=data.frame(),pos_min=NA,pos_max=NA))
      chr_t <- as.character(lnc_info$chr[1])
      pmin <- min(lnc_info$start)
      pmax <- max(lnc_info$end)
      bind_dt <- get_protein_binding_data(chr_filter=chr_t,pos_min=pmin,pos_max=pmax)
      list(res=df,prot_bind=bind_dt,exons=lnc_info %>% filter(region=="exon"),pos_min=pmin,pos_max=pmax)
    }, error=function(e){
      showNotification(paste("Analysis failed:",e$message),type="error")
      NULL
    })
  })

  # Wordcloud
  wordcloud_plot <- reactive({
    req(analysis_result())
    res <- analysis_result()$res
    freq <- as.data.frame(table(res$Protein_name))
    colnames(freq) <- c("word","freq")
    if(nrow(freq)>100) freq <- head(freq[order(-freq$freq),],100)
    wordcloud2(freq, size=0.8, backgroundColor="white")
  })
  output$wordcloud <- renderWordcloud2(wordcloud_plot())
  output$download_wordcloud <- downloadHandler("wordcloud.html", function(file) saveWidget(wordcloud_plot(), file))

  # RBP Lollipop plot
  RBP_plot <- reactive({
    req(analysis_result())
    dat <- analysis_result()
    bind <- dat$prot_bind
    ex <- dat$ex
    if(nrow(bind)==0 || is.na(dat$pos_min)) return(ggplot()+annotate("text",0.5,0.5,"No RBP binding data",size=6,color="red")+theme_void())
    ggplot() + 
      geom_hline(yintercept = 0.1, linewidth = 1.5, color = "black", linetype = "solid") +
      geom_rect(data=ex,aes(xmin=start,xmax=end,ymin=0,ymax=0.2),fill=rep("#7EB7DC", nrow(ex))) +
      geom_segment(data=bind,aes(x=position,xend=position,y = 0.2, yend = 1),linewidth = 0.3, colour = "black") +
      geom_point(data=bind,aes(x=position, y = 1),size = 3, alpha = 0.7,fill="#E64B35",shape=21) +
      geom_text(data=bind,aes(x=position,y=1.08,label=protein_name), family = "DejaVu Sans",size = 6,angle=90,hjust=0, color = "black") +
      xlim(dat$pos_min,dat$pos_max) + ylim(-0.2, 1.8) + theme_void() + labs(title=input$ana_rna)+
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5, size = 16, family = "DejaVu Sans", face = "bold"))
  })
  output$lolliplot <- renderPlot(RBP_plot())
  output$download_lolliplot <- downloadHandler("rbp_plot.png",function(file) 
  ggsave(file,RBP_plot(),w=10,h=6,dpi=600))

  # PPI Plot with error capture
  PPI_plot <- reactive({
    req(analysis_result())   
    proteins_freq <- as.data.frame(table(analysis_result()$res$Protein_name))
    colnames(proteins_freq) <- c("protein_symbol","freq")
    proteins_freq <- proteins_freq[which(proteins_freq$freq > 0), ]
    proteins_freq$protein_symbol <- as.character(proteins_freq$protein_symbol)
    if(nrow(proteins_freq) == 0) {
        return(ggplot() + 
                annotate("text", x = 0.5, y = 0.5, 
                        label = "Non-interacting proteins.", size = 6, color = "red") +
                theme_void())
    }

    safe_ppi <- tryCatch(getPPI(proteins_freq$protein_symbol, taxID="9606"), error=function(e) NULL)
    if(is.null(safe_ppi) || length(V(safe_ppi)) == 0 || length(E(safe_ppi)) == 0) {
        return(ggplot() + 
                annotate("text", x = 0.5, y = 0.5, 
                        label = "Non-interacting proteins.", size = 6, color = "red") +
                theme_void())
    }
    degree_centrality <- degree(safe_ppi, mode = "all")
    betweenness_centrality <- betweenness(safe_ppi, directed = FALSE)
    closeness_centrality <- closeness(safe_ppi, mode = "all")
    eigen_centrality <- eigen_centrality(safe_ppi, directed = FALSE)$vector
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
    ggplot(safe_ppi, layout='circle') %<+% node_attributes + 
    geom_edge(aes(alpha = 0.3), color = "gray70") + 
    geom_point(aes(color = node_color, size = node_size), alpha = 0.9) + 
    shadowtext::geom_shadowtext(aes(label = name, color = label_color),  family = "DejaVu Sans",  bg.color = "white", size = 3.5) +
    scale_color_identity() + scale_size_identity() +
    annotate("text", x = -Inf, y = Inf, label = paste0("★ Hub Protein (N = ", sum(is_hub), ")"), 
            hjust = -0.1, vjust = 1.5,  color = "#E64B35", fontface = "bold", size = 5, family = "DejaVu Sans") +
    theme_void() +
    theme(text = element_text(family = "DejaVu Sans", size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(family = "DejaVu Sans", face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(family = "DejaVu Sans", size = 9, color = "gray50", hjust = 0.5))
  })
  output$network <- renderPlot(PPI_plot())
  output$download_network <- downloadHandler("ppi.png",
    function(file)
    ggsave(file,PPI_plot(),w=12,h=10)
  )

  # GO Dotplot
  GO_plot <- reactive({
    req(analysis_result())
    prots <- unique(analysis_result()$res$Protein_name)
    prots <- prots[nchar(prots)>0]
    go_obj <- enrichGO(prots, keyType="SYMBOL", OrgDb=org.Hs.eg.db, ont="ALL")
    dotplot(go_obj, showCategory = 15, 
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

  })
  output$godotplot <- renderPlot(GO_plot())
  output$download_godotplot <- downloadHandler("go.png",function(file) 
  ggsave(file,GO_plot(),w=12,h=8))

  # Analysis tables
  output$analysis_table <- renderDT({
    req(analysis_result())
    if (is.null(analysis_result()$res) || nrow(analysis_result()$res) == 0) return(NULL)
    datatable(analysis_result()$res[,c("lncRNA_Name","Protein_name","Cell_Line","Method","Data")], rownames=F,
              colnames=c("lncRNA","Protein","Cell Line","Method","Reference"))
  })
  output$download_table_csv <- downloadHandler("analysis.csv",function(file) write.csv(analysis_result()$res,file,row.names=F))
  output$analysis_table2 <- renderDT({
    req(analysis_result())
    datatable(analysis_result()$prot_bind[,c("chr","position","protein_name")], rownames=F,
              colnames=c("Chr","Binding Position","Protein"))
  })
  output$download_table_csv2 <- downloadHandler("rbp_bind.csv",function(file) write.csv(analysis_result()$prot_bind,file,row.names=F))


  
  # Bulk download handlers
  output$downloas_all_data <- downloadHandler(
    filename = "HuRInterDB_all.csv",
    content = function(file){
      all_dt <- fetch_filtered_data()
      cols <- c("lncRNA_Name","Protein_name","Cell_Line","Method")
      cols <- cols[cols %in% colnames(all_dt)]
      write.csv(all_dt[, cols], file, row.names = FALSE)
    }
  )
  output$downloas_TF_data <- downloadHandler(
    "HuRInterDB_TF.csv",
    function(file) write.csv(read_parquet("data/tf_data.parquet"), file, row.names=F)
  )
  output$downloas_RBP_data <- downloadHandler(
    "HuRInterDB_RBP.csv",
    function(file) write.csv(read_parquet("data/rbp_data.parquet"), file, row.names=F)
  )
})

shinyApp(ui, server)