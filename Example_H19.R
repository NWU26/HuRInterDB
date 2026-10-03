library(tidyverse)
library(wordcloud2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggtangle)
library(ggplot2)
library(igraph)


# Wordcloud
freq <- as.data.frame(table(analysis_result$Protein_name))
colnames(freq) <- c("word","freq")
wordcloud2(freq, size=0.8, backgroundColor="white")

# RBP Lollipop plot
bind <- transcript_binding_data$prot_bind
ex   <- transcript_binding_data$exons
ggplot() +
  geom_hline(yintercept = 0.1, linewidth = 1.5, color = "black", linetype = "solid") +
  geom_rect(data = ex,
            aes(xmin = start, xmax = end, ymin = 0, ymax = 0.2),
            fill = "#7EB7DC") +
  geom_segment(data = bind,
              aes(x = Position, xend = Position, y = 0.2, yend = 1),
              linewidth = 0.3, colour = "black") +
  geom_point(data = bind,
            aes(x = Position, y = 1),
            size = 3, alpha = 0.7, fill = "#E64B35", shape = 21) +
  geom_text(data = bind,
            aes(x = Position, y = 1.08, label = Protein_name),
            family = "DejaVu Sans", size = 5, angle = 90,
            hjust = 0, color = "black") +
  xlim(transcript_binding_data$pos_min, transcript_binding_data$pos_max) +
  ylim(-0.2, 1.8) +
  theme_void() +
  labs(title = paste0(analysis_result()$rna, "  (", transcript_binding_data$transcript_id, ")")) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, size = 14, family = "DejaVu Sans", face = "bold"))

# PPI Plot with error capture
proteins_freq <- as.data.frame(table(analysis_result$Protein_name))
colnames(proteins_freq) <- c("protein_symbol","freq")
proteins_freq$protein_symbol <- as.character(proteins_freq$protein_symbol)
result <- getPPI(batch, taxID="9606", add_nodes = 0, required_score = 700)
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

# GO Dotplot
prots <- unique(analysis_result$Protein_name)
go_obj <- enrichGO(prots, keyType="SYMBOL", 
                   OrgDb=org.Hs.eg.db, 
                   universe = all_protein, 
                   ont="ALL",
                   pAdjustMethod = "BH", 
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)
dotplot(go_obj, showCategory = 15, color = "pvalue", label_format = 80, font.size = 11, title = "GO Pathway Enrichment") +
        theme(text = element_text(size = 12),
              axis.text.y = element_text(family = "DejaVu Sans", size = 10),
              axis.text.x = element_text(family = "DejaVu Sans", size = 9),
              axis.title = element_text(family = "DejaVu Sans", size = 12),
              legend.text = element_text(family = "DejaVu Sans"),
              legend.title = element_text(family = "DejaVu Sans"),
              plot.background = element_rect(fill = "white", color = NA),
              plot.title = element_text(family = "DejaVu Sans", face = "bold"))
