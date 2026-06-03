shinyServer(function(input, output, session){
##---- HOME ----
            current_results <- reactiveVal(NULL)
            ##-- RNA search  ----
            observeEvent(input$RNA_search, {
                        query <- input$RNA_search
                        filtered <- filter(data, lncRNA_Name == query)
                        filtered <- prepare_data_with_links(filtered)
                        current_results(filtered)
                        updateTabsetPanel(session = session, inputId = "navbar", selected = "searchs")
            })

            ##-- RNA search  ----
            observeEvent(input$protein_search, {
                        query <- input$protein_search
                        filtered <- filter(data, Protein_name == query)
                        filtered <- prepare_data_with_links(filtered)
                        current_results(filtered)
                        updateTabsetPanel(session = session, inputId = "navbar", selected = "searchs")
            })

##---- SEARCH ----
            # === Search Page: Apply filters ===
            search_data <- eventReactive(input$apply_filter, {
            req(data)
            
            df <- data
            
            gene_f <- trimws(input$filter_gene)
            prot_f <- trimws(input$filter_protein)
            cell_f <- trimws(input$filter_cellline)
            meth_f <- input$filter_method
                        
            if (gene_f != "") 
                df <- filter(df, lncRNA_Name == gene_f)
            if (prot_f != "") 
                df <- filter(df, Protein_name == prot_f)
            if (cell_f != "") 
                df <- filter(df, Cell_Line == cell_f)
            if (meth_f != "") 
                df <- df[df$Method == meth_f, ]
            
            df
            df <- prepare_data_with_links(df)
            }, ignoreNULL = FALSE)

            observeEvent(input$apply_filter, {
                        current_results(search_data())
            })

            # Render table in Search page
            output$result_table <- DT::renderDataTable({
                                                        res <- current_results()
                                                        if (is.null(res) || nrow(res) == 0) return(NULL)
                                                        
                                                        display_df <- unique(res[, c("lncRNA_Name_link", "lncRNA_RNALocate_link", "Protein_name_link", "Protein_Domains_link", "AlphaFoldDB_link", 
                                                                               "KEGG_link", "Cell_Line", "Method", "Data_link")])
                                                        display_df <- subset(display_df, !is.na(Method))
                                                        colnames(display_df) <- c("lncRNA", "lncRNA Localization", "Protein", "Protein Domains", "Protein Structure", "Protein KEGG", "Cell Line", "Method", "Data")
                                                        datatable(display_df, escape = FALSE,
                                                                  options = list(pageLength = 10, lengthMenu = c(10, 25, 50)),
                                                                  rownames = FALSE,filter = "top"
                                                        )
            })

            # Download filtered data (plain text, no HTML links)
            output$download_data <- downloadHandler(
                                    filename = function() "search_results.csv",
                                    content = function(file) {
                                        res <- current_results()
                                        if (!is.null(res)) {
                                        write.csv(res[c("lncRNA_Name", "Protein_name", "AlphaFoldDB", "KEGG", 
                                                        "Cell_Line", "Method", "Data")], 
                                                    file, row.names = FALSE)
                                        } else {
                                        write.csv(data.frame(), file, row.names = FALSE)
                                        }
                                    }
            )

##---- RPI Analysis ----
                # === Analysis Page ===
                analysis_result <- eventReactive(input$analyze_btn, {
                                    rna_name <- trimws(input$lncrna_input)
                                    if (rna_name == "") return(NULL)
                                    res <- filter(data, lncRNA_Name == rna_name)
                                    res <- res[,c("lncRNA_Name", "Protein_name", "Cell_Line", "Method", "Data")]
                                    res
                })
                

                # Generate wordcloud from Protein_name
                output$wordcloud <- renderWordcloud2({
                                    res <- analysis_result()
                                    if (is.null(res) || nrow(res) == 0) return(NULL)
                                    
                                    proteins <- res$Protein_name
                                    freq_table <- head(sort(table(proteins), decreasing = TRUE), n = 20)
                                    wordcloud2(data = data.frame(word = names(freq_table), freq = as.numeric(freq_table)),
                                            color = "random-light", backgroundColor = "white")
                })

                # Generate lolliplot
                output$lolliplot <- renderPlot({
                                                req(input$lncrna_input)
                                                res <- analysis_result()

                                                load(file = "data/lncRNA_bed_data.RData")
                                                exons <- lncRNA_bed_data %>%
                                                        filter(lncRNA_name == input$lncrna_input, region == "exon")
                                                if (nrow(exons) == 0) return(NULL)
                                                exons_gr <- makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE,
                                                            seqnames.field = "chr", start.field = "start", end.field = "end"
                                                )
                                                num_exons <- nrow(exons)
                                                exons_gr$height <- rep(list(c(0.1)), length(num_exons))
                                                exon_colors <- sample(colors(), num_exons)
                                                exons_gr$fill <- list(exon_colors)

                                                lncRNA_info <- lncRNA_bed_data %>%
                                                            filter(lncRNA_name == input$lncrna_input)
                                                if (nrow(lncRNA_info) == 0) return(NULL)
                                                chr_target <- as.character(lncRNA_info$chr[1])
                                                range <- c(min(lncRNA_info$start), max(lncRNA_info$end))
                                                rm(lncRNA_bed_data)
                                                
                                                load(file = "data/protein_binding_data.RData")
                                                prot_bindings <- protein_binding_data %>%
                                                                    filter(chr == chr_target, between(position, range[1], range[2]))
                                                if (nrow(prot_bindings) == 0) return(NULL)
                                                rm(protein_binding_data)
                                                proteins <- unique(prot_bindings$protein_name)
                                                protein_gr <- GRanges(chr_target, IRanges(prot_bindings$position, width=1, names=c(prot_bindings$protein_name)))

                                                dandelion.plot(protein_gr, exons_gr,type="circle",ylab="",yaxis=FALSE, xaxis = F)
                                                grid.text(input$lncrna_input, x=.5, y=.98, just="top",  gp=gpar(cex=1.5, fontface="bold"))
                                            })
                
                # Generate PPI from Protein_name
                  output$network <- renderPlot({res <- analysis_result()
                                               proteins_freq <- as.data.frame(table(res$Protein_name))
                                                colnames(proteins_freq) <- c("protein_symbol","freq")
                                                proteins_freq <- proteins_freq[which(proteins_freq$freq > 1), ]
                                                proteins_freq$protein_symbol <- as.character(proteins_freq$protein_symbol)

                                                g <- getPPI(proteins_freq$protein_symbol, taxID="9606")
                                                ggplot(g, layout='circular') %<+% proteins_freq + 
                                                geom_edge() + 
                                                geom_point(aes(color = freq), size = 8) + 
                                                shadowtext::geom_shadowtext(aes(label = name), color="black", bg.color="white") +
                                                enrichplot::set_enrichplot_color(reverse=F) +
                                                theme_void()
                                                })

                # Generate GO from Protein_name
                output$godotplot <- renderPlot({
                    res <- analysis_result()
                    proteins_list <- unique(res$Protein_name)
                    proteins_list <- proteins_list[nchar(proteins_list) > 0]

                    go_result <- enrichGO(gene = proteins_list,
                                          OrgDb = org.Hs.eg.db,
                                          keyType = "SYMBOL",
                                          ont = "ALL",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.99,
                                          qvalueCutoff = 0.99)                    
                    dotplot(go_result, showCategory = 15, 
                            color = "pvalue",
                            label_format = 50,
                            font.size = 11,
                            title = "GO Pathway Enrichment") +
                           theme(text = element_text(size = 12))
                })

                # Display detailed table (same format as Search page)
                output$analysis_table <- DT::renderDataTable({
                    res <- analysis_result()
                    if (is.null(res) || nrow(res) == 0) return(NULL)
                    
                    datatable(res, escape = FALSE, options = list(pageLength = 10), rownames = FALSE)
                })

                ##---- Download table ----
                output$download_table_csv <- downloadHandler(
                                            filename = function() {
                                                    paste0(input$lncrna_input, "_binding_protein.csv")
                                            },
                                            content = function(file) {
                                               res <- analysis_result()

                                                if(!is.null(res) && nrow(res) > 0) {
                                                write.csv(res, file, row.names = FALSE)
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
                                        write.csv(data[,c("lncRNA_Name","Protein_name","Cell_Line","Method")], file)})
            ##---- Download TF ----
            output$downloas_TF_data <- downloadHandler(
                                        filename = function() {"HuRInterDB_TF_only_RIPs_data.csv"},
                                        content = function(file) {
                                        TF_data <- data %>%
                                        filter(TF == "TF") %>%
                                        select(lncRNA_Name, Protein_name, Cell_Line, Method)
                                        write.csv(TF_data, file)})
            ##---- Download RBP ----
            output$downloas_RBP_data <- downloadHandler(
                                        filename = function() {"HuRInterDB_RBP_only_RIPs_data.csv"},
                                        content = function(file) {
                                        RBP_data <- data %>%
                                        filter(RBP == "RBP") %>%
                                        select(lncRNA_Name, Protein_name, Cell_Line, Method)
                                        write.csv(RBP_data, file)})
})


