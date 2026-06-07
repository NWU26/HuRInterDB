shinyUI(
  fluidPage(
    ##-- Favicon ----
    tags$head(
      tags$link(rel = "shortcut icon", href = "img/logo.ico"),
      #-- biblio js ----
      tags$link(rel="stylesheet", type = "text/css",
                href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"),
      tags$link(rel="stylesheet", type = "text/css",
                href = "https://fonts.googleapis.com/css?family=Open+Sans|Source+Sans+Pro")
    ),
    ##-- Logo ----
    list(tags$head(HTML('<link rel="icon", href="img/logo.png",type="image/png" />'))),
    div(style="padding: 1px 0px; width: '50%'",
        titlePanel(
          title ="HuRInterDB", windowTitle = "HuRInterDB"
        )
    ),
    ##-- Header ----
    navbarPage(title = div(img(src="img/logo.png", height = "50px"), style = "padding-left:40px;"),
               id = "navbar", selected = "Home", theme = "styles.css", fluid = T,
###------ Home  ------###
              home <- tabPanel(title = "Home", value = "home", icon = icon("house"),
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
                            )),

                            ##-- Statistics ----
                            fluidRow(
                                   column(4, div(class = "stat-card",
                                                 tags$i(class = "fas fa-dna fa-2x", style = "margin-bottom: 1px;"),
                                                 h2("55925 lncRNA", style = "margin: 0; font-size: 24px;"))
                                          ),
                                   column(4, div(class = "stat-card",
                                                 tags$i(class = "fas fa-project-diagram fa-2x", style = "margin-bottom: 1px;"),
                                                 h2("7864 Protein", style = "margin: 0; font-size: 24px;"))
                                          ),
                                   column(4, div(class = "stat-card",
                                                 tags$i(class = "fas fa-atom fa-2x", style = "margin-bottom: 1px;"),
                                                 h2("2195828 Interaction", style = "margin: 0; font-size: 24px;"))
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
                                          list(name = "miRBase", url = "http://mirbase.org"),
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
              search <- tabPanel(title = "Search", value = "searchs", icon = icon("magnifying-glass"),
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
                                   ##-- Outputs ----
                                   titlePanel(h3("Search results:", style = "color: #0277bd; text-align: left;")) ,
                                   withSpinner(dataTableOutput("result_table"), type = 6, color = "#0277bd"),
                                   ##-- Download ----
                                   downloadButton("download_data", "Download Results")
                            )
              ),
###------ Analysis  ------###
              RPIanalysis <- tabPanel(title = "RPIanalysis", value = "RPIanalysis", icon = icon("atom"),
                                   br(), hr(), 
                                   tags$head(
                                   tags$style("
                                   .inputBox {background-color: #E1F5FE; padding: 10px; margin-bottom: 10px;
                                              box-shadow: 0 2px 4px rgba(0,0,0,0.1);}
                                   .resultBox {background-color: #E1F5FE; padding: 10px; border-radius: 10px; 
                                          margin-top: 10px; margin-bottom: 10px;}
                                   ")
                                   ),
                                   ##---- Head ----
                                   titlePanel(h2("🔍 On-line lncRNA-protein interaction analysis", style = "color: #0277bd; text-align: center;")),
                                   ##---- Input ----
                                   div(id = "inputBox", h3("📌 Input your lncRNA list:"),
                                       textAreaInput(inputId= "lncrna_input", label = NULL, 
                                                     width = '100%', rows = 3, resize = "horizontal"),
                                       actionButton(inputId = "analyze_btn", label = "Continue", icon = icon("arrow-right"),class = "btn-primary btn-lg")
                                   ),
                                   ##-- Outputs ----
                                   # Fig result
                                   div(class = "resultBox", h3("📊 Interaction Analysis Results"),
                                   tabsetPanel(
                                          tabPanel("Wordcloud", withSpinner(wordcloud2Output("wordcloud", width = "100%", height = "400px"), type = 6, color = "#0277bd")),
                                          tabPanel("lolliplot", withSpinner(plotOutput("lolliplot", width = "100%", height = "400px"), type = 6, color = "#0277bd")),
                                          tabPanel("PPI", withSpinner(plotOutput('network', width = "100%", height = "400px"), type = 6, color = "#0277bd")),
                                          tabPanel("GO", withSpinner(plotOutput("godotplot", width = "100%", height = "400px"), type = 6, color = "#0277bd"))
                                          )
                                   ),
                                   # Table result
                                   div(class = "resultBox",
                                          h3("📋 List of Interacting Proteins"),
                                          withSpinner(DT::dataTableOutput("analysis_table"), type = 6, color = "#0277bd"),
                                          downloadLink(outputId = "download_table_csv", icon = icon("download"), label = "Download (CSV)")
                                   )
              ),
              
###------ Download  ------###
              download <- tabPanel(title = "Download", value = "download", icon = icon("download"),
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
              about <- tabPanel(title = "About", value = "about", icon = icon("ghost"),
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


