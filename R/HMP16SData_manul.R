## --------------------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE)

## --------------------------------------------------------------------------
library(HMP16SData)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(tibble)
library(dplyr)
library(dendextend)
library(circlize)
library(curatedMetagenomicData)
library(gridExtra)
library(cowplot)
library(readr)
library(haven)

## ---- eval=FALSE-----------------------------------------------------------
#  BiocManager::install("HMP16SData")

## --------------------------------------------------------------------------
V13()
V35()

## --------------------------------------------------------------------------
V13() %>%
    table_one() %>%
    head()

## --------------------------------------------------------------------------
list(V13 = V13(), V35 = V35()) %>%
    table_one() %>%
    kable_one()

## --------------------------------------------------------------------------
V35_stool <-
    V35() %>%
    subset(select = HMP_BODY_SUBSITE == "Stool")

V35_stool

## ---- eval=FALSE-----------------------------------------------------------
#  V35_stool_protected <-
#      V35_stool %>%
#      attach_dbGaP("~/prj_12146.ngc")

## ---- eval=FALSE-----------------------------------------------------------
#  colData(V35_stool_protected)

## --------------------------------------------------------------------------
V13_tonsils <-
    V13() %>%
    subset(select = HMP_BODY_SUBSITE == "Palatine Tonsils")

V13_stool <-
    V13() %>%
    subset(select = HMP_BODY_SUBSITE == "Stool")

## --------------------------------------------------------------------------
V13_tonsils_phyloseq <-
    as_phyloseq(V13_tonsils)

V13_stool_phyloseq <-
    as_phyloseq(V13_stool)

## --------------------------------------------------------------------------
sample_samples <- function(x, size) {
    sampled_names <-
        sample_names(x) %>%
        sample(size)

    prune_samples(sampled_names, x)
}

## --------------------------------------------------------------------------
V13_tonsils_phyloseq %<>%
    sample_samples(25)

V13_stool_phyloseq %<>%
    sample_samples(25)

## --------------------------------------------------------------------------
sample_data(V13_tonsils_phyloseq)$Study <- "Tonsils"

sample_data(V13_stool_phyloseq)$Study <- "Stool"

## --------------------------------------------------------------------------
V13_phyloseq <-
    merge_phyloseq(V13_tonsils_phyloseq, V13_stool_phyloseq)

## --------------------------------------------------------------------------
V13_phyloseq %<>%
    taxa_sums() %>%
    is_greater_than(0) %>%
    prune_taxa(V13_phyloseq)

## --------------------------------------------------------------------------
richness_measures <-
    c("Observed", "Shannon", "Simpson")

## ---- fig.height=5, fig.width=8--------------------------------------------
V13_phyloseq %>%
    plot_richness(x = "Study", color = "Study", measures = richness_measures) +
    stat_boxplot(geom ="errorbar") +
    geom_boxplot() +
    theme_bw() +
    theme(axis.title.x = element_blank(), legend.position = "none")

## --------------------------------------------------------------------------
V13_dendrogram <-
    distance(V13_phyloseq, method = "bray") %>%
    hclust() %>%
    as.dendrogram()

## --------------------------------------------------------------------------
V13_sample_data <-
    sample_data(V13_phyloseq) %>%
    data.frame()

## --------------------------------------------------------------------------
V13_sample_data %<>%
    rownames_to_column(var = "PSN") %>%
    mutate(labels_col = if_else(Study == "Stool", "#F8766D", "#00BFC4")) %>%
    mutate(leaves_col = if_else(Study == "Stool", "#F8766D", "#00BFC4")) %>%
    mutate(leaves_pch = if_else(Study == "Stool", 16, 17))

## --------------------------------------------------------------------------
V13_sample_order <-
    labels(V13_dendrogram) %>%
    match(V13_sample_data$PSN)

## --------------------------------------------------------------------------
labels_col <- V13_sample_data$labels_col[V13_sample_order]
leaves_col <- V13_sample_data$leaves_col[V13_sample_order]
leaves_pch <- V13_sample_data$leaves_pch[V13_sample_order]

## --------------------------------------------------------------------------
V13_dendrogram %<>%
    set("labels_col", labels_col) %>%
    set("leaves_col", leaves_col) %>%
    set("leaves_pch", leaves_pch)

## ---- fig.height=8, fig.width=8--------------------------------------------
V13_dendrogram %>%
    circlize_dendrogram(labels_track_height = 0.2)

## --------------------------------------------------------------------------
V13_ordination <-
    ordinate(V13_phyloseq, method = "PCoA", distance = "bray")

## ---- fig.height=5, fig.width=8--------------------------------------------
V13_phyloseq %>%
    plot_ordination(V13_ordination, color = "Study", shape = "Study") +
    theme_bw() +
    theme(legend.position = "bottom")

## ---- fig.height=5, fig.width=8--------------------------------------------
V13_ordination %>%
    plot_scree() +
    theme_bw()

## --------------------------------------------------------------------------
V35_stool_phyloseq <-
    V35_stool %>%
    as_phyloseq()

MGX_stool_phyloseq <-
    HMP_2012.metaphlan_bugs_list.stool(cmdversion = "20170526") %>%
    ExpressionSet2phyloseq()

## --------------------------------------------------------------------------
MGX_stool_melted <-
    MGX_stool_phyloseq %>%
    subset_taxa(is.na(Class) & !is.na(Phylum)) %>%
    psmelt()

## --------------------------------------------------------------------------
V35_stool_melted <-
    V35_stool_phyloseq %>%
    tax_glom(taxrank = "PHYLUM") %>%
    psmelt()

## --------------------------------------------------------------------------
SRS_SAMPLE_IDS <-
    intersect(V35_stool_melted$SRS_SAMPLE_ID, MGX_stool_melted$NCBI_accession)

V35_stool_melted %<>%
    filter(SRS_SAMPLE_ID %in% SRS_SAMPLE_IDS) %>%
    rename(Phylum = PHYLUM) %>%
    mutate(Sample = SRS_SAMPLE_ID) %>%
    group_by(Sample) %>%
    mutate(`Relative Abundance` = Abundance / sum(Abundance)) %>%
    arrange(Phylum, -`Relative Abundance`) %>%
    group_by(Phylum) %>%
    mutate(phylum_rank = sum(`Relative Abundance`)) %>%
    arrange(-phylum_rank) %>%
    select(Sample, Phylum, `Relative Abundance`) %>%
    ungroup()

MGX_stool_melted %<>%
    filter(NCBI_accession %in% SRS_SAMPLE_IDS) %>%
    group_by(Sample) %>%
    mutate(`Relative Abundance` = Abundance / sum(Abundance)) %>%
    arrange(Phylum, -`Relative Abundance`) %>%
    group_by(Phylum) %>%
    mutate(phylum_rank = sum(`Relative Abundance`)) %>%
    arrange(-phylum_rank) %>%
    select(Sample, Phylum, `Relative Abundance`) %>%
    ungroup()

## --------------------------------------------------------------------------
V35_top_phyla <-
    V35_stool_melted %$%
    unique(Phylum) %>%
    as.character()

MGX_top_phyla <-
    MGX_stool_melted %$%
    unique(Phylum) %>%
    as.character()

top_eight_phyla <-
    intersect(V35_top_phyla, MGX_top_phyla) %>%
    extract(1:8)

V35_stool_melted %<>%
    filter(Phylum %in% top_eight_phyla)

MGX_stool_melted %<>%
    filter(Phylum %in% top_eight_phyla)

## --------------------------------------------------------------------------
sample_order <-
    V35_stool_melted %$%
    unique(Sample)

## --------------------------------------------------------------------------
phylum_order <-
    V35_stool_melted %$%
    unique(Phylum)

## --------------------------------------------------------------------------
bang_wong_colors <-
    c("#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73", "#56B4E9",
      "#E69F00", "#000000")

## --------------------------------------------------------------------------
V35_plot <-
    V35_stool_melted %>%
    mutate(Sample = factor(Sample, sample_order)) %>%
    mutate(Phylum = factor(Phylum, phylum_order)) %>%
    ggplot(aes(Sample, `Relative Abundance`, fill = Phylum)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    scale_fill_manual(values = bang_wong_colors) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none",
          legend.title = element_blank()) +
    ggtitle("Phylum-Level Relative Abundance", "16S Stool Samples")

MGX_plot <-
    MGX_stool_melted %>%
    mutate(Sample = factor(Sample, sample_order)) %>%
    mutate(Phylum = factor(Phylum, phylum_order)) %>%
    ggplot(aes(Sample, `Relative Abundance`, fill = Phylum)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    scale_fill_manual(values = bang_wong_colors) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none",
          legend.title = element_blank()) +
    ggtitle("", "MGX Stool Samples")

## --------------------------------------------------------------------------
plot_legend <- {
        MGX_plot +
            theme(legend.position = "bottom")
    } %>%
    get_legend()

## ---- fig.height=8, fig.width=8--------------------------------------------
grid.arrange(V35_plot, MGX_plot, plot_legend, ncol = 1, heights = c(3, 3, 1))

## ---- eval=FALSE-----------------------------------------------------------
#  V35_plot +
#      theme(text = element_text(size = 19)) +
#      labs(title = NULL, subtitle = NULL, tag = "A")
#  
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2A.eps", device = "eps")
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2A.pdf", device = "pdf")
#  
#  MGX_plot +
#      theme(text = element_text(size = 19)) +
#      labs(title = NULL, subtitle = NULL, tag = "B")
#  
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2B.eps", device = "eps")
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2B.pdf", device = "pdf")
#  
#  plot_legend <- {
#          MGX_plot +
#              theme(legend.position = "bottom", text = element_text(size = 19)) +
#              guides(fill = guide_legend(byrow = TRUE))
#      } %>%
#      get_legend()
#  
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2 Legend 1.eps", plot = plot_legend,
#         device = "eps")
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2 Legend 1.pdf", plot = plot_legend,
#         device = "pdf")
#  
#  plot_legend <- {
#          MGX_plot +
#              theme(legend.position = "right", text = element_text(size = 19))
#      } %>%
#      get_legend()
#  
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2 Legend 2.eps", plot = plot_legend,
#         device = "eps")
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2 Legend 2.pdf", plot = plot_legend,
#         device = "pdf")
#  
#  plot_legend <- {
#          MGX_plot +
#              theme(legend.position = "right", text = element_text(size = 19)) +
#              guides(fill = guide_legend(ncol = 2, byrow = TRUE))
#      } %>%
#      get_legend()
#  
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2 Legend 3.eps", plot = plot_legend,
#         device = "eps")
#  ggsave("~/AJE-00611-2018 Schiffer Figure 2 Legend 3.pdf", plot = plot_legend,
#         device = "pdf")
#  
#  V35_plot <-
#      V35_plot +
#      theme(text = element_text(size = 19)) +
#      labs(title = NULL, subtitle = NULL, tag = "A")
#  
#  MGX_plot <-
#      MGX_plot +
#      theme(text = element_text(size = 19)) +
#      labs(title = NULL, subtitle = NULL, tag = "B")
#  
#  plot_legend <- {
#          MGX_plot +
#              theme(legend.position = "bottom", text = element_text(size = 19)) +
#              guides(fill = guide_legend(byrow = TRUE))
#      } %>%
#      get_legend()
#  
#  grid_object <-
#      grid.arrange(V35_plot, MGX_plot, plot_legend, ncol = 1,
#                   heights = c(3, 3, 1))
#  
#  ggsave("~/AJE-00611-2018 Schiffer Figure 3.eps", plot = grid_object,
#         device = "eps", width = 8, height = 8)
#  ggsave("~/AJE-00611-2018 Schiffer Figure 3.pdf", plot = grid_object,
#         device = "pdf", width = 8, height = 8)

## ---- eval=FALSE-----------------------------------------------------------
#  V13_participant_data <-
#    V13() %>%
#    colData() %>%
#    as.data.frame() %>%
#    rownames_to_column(var = "PSN")

## ---- eval=FALSE-----------------------------------------------------------
#  V13_OTU_counts <-
#    V13() %>%
#    assay() %>%
#    t() %>%
#    as.data.frame() %>%
#    rownames_to_column(var = "PSN")

## ---- eval=FALSE-----------------------------------------------------------
#  V13_data <-
#    merge.data.frame(V13_participant_data, V13_OTU_counts, by = "PSN")

## ---- eval=FALSE-----------------------------------------------------------
#  V13_dictionary <-
#    V13() %>%
#    rowData() %>%
#    t.data.frame() %>%
#    as.data.frame()

## ---- eval=FALSE-----------------------------------------------------------
#  colnames(V13_data) <-
#      colnames(V13_data) %>%
#      gsub(pattern = "\\.", replacement = "_", x = .)

## ---- eval=FALSE-----------------------------------------------------------
#  colnames(V13_dictionary) <-
#      colnames(V13_dictionary) %>%
#      gsub(pattern = "\\.", replacement = "_", x = .)

## ---- eval=FALSE-----------------------------------------------------------
#  write_csv(V13_data, "~/V13_data.csv")
#  write_csv(V13_dictionary, "~/V13_dictionary.csv")

## ---- eval=FALSE-----------------------------------------------------------
#  write_sas(V13_data, "~/V13_data.sas7bdat")
#  write_sas(V13_dictionary, "~/V13_dictionary.sas7bdat")

## ---- eval=FALSE-----------------------------------------------------------
#  write_sav(V13_data, "~/V13_data.sav")
#  write_sav(V13_dictionary, "~/V13_dictionary.sav")

## ---- eval=FALSE-----------------------------------------------------------
#  write_dta(V13_data, "~/V13_data.dta", version = 13)
#  write_dta(V13_dictionary, "~/V13_dictionary.dta", version = 13)

## --------------------------------------------------------------------------
sessionInfo()

