# _targets.R

#Load libraries
library(targets)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(matrixStats)
library(pheatmap)
library(glue)

#targets options
tar_option_set(
  packages = c("tidyverse","DESeq2","edgeR","matrixStats","pheatmap","glue"),
  format   = "rds"
)

#pick top n variable genes
select_top_variable <- function(vsd, n = 50) {
  rv <- rowVars(assay(vsd))
  head(order(rv, decreasing = TRUE), n)
}

#Pipeline
pipeline <- list(
  
  # 3.1) Data files
  tar_target(counts_file, "data/GSE60424_raw_counts_GRCh38.p13_NCBI.tsv", format = "file"),
  tar_target(meta_file,   "data/SraRunTable-3.csv",                    format = "file"),
  
  # Read & clean counts
  tar_target(
    counts_df,
    {
      df  <- read_tsv(counts_file, col_types = cols())
      df2 <- df %>%
        column_to_rownames("gene_id") %>%
        select(where(~ !any(is.na(.))))
      df2
    }
  ),
  
  #Raw metadata
  tar_target(
    meta_raw,
    read_csv(meta_file, show_col_types = FALSE)
  ),
  
  #Build col_data from Sample Name & celltype
  tar_target(
    col_data,
    {
      cd <- meta_raw %>%
        transmute(
          sample    = `Sample Name`,
          condition = celltype %>% as.character() %>% str_trim()
        ) %>%
        filter(
          sample    %in% colnames(counts_df),
          !is.na(condition),
          condition != ""
        )
      if (nrow(cd) == 0) {
        stop("No samples survived filtering on Sample Name & celltype.")
      }
      cd %>% column_to_rownames("sample")
    }
  ),
  
  #  samples per condition (group_by + summarise)
  tar_target(
    cond_table,
    {
      tibble(
        sample    = rownames(col_data),
        condition = col_data$condition
      ) %>%
        group_by(condition) %>%
        summarise(n = n(), .groups = "drop")
    }
  ),
  
  # Build & run DESeq2
  tar_target(
    dds,
    {
      ct  <- counts_df[, rownames(col_data), drop = FALSE]
      dds <- DESeqDataSetFromMatrix(
        countData = ct,
        colData   = col_data,
        design    = ~ condition
      )
      DESeq(dds)
    }
  ),
  
  #Median‐of‐ratios normalization
  tar_target(
    norm_median,
    as_tibble(counts(dds, normalized = TRUE), rownames = "gene")
  ),
  
  #  TMM normalization (edgeR)
  tar_target(
    norm_tmm,
    {
      dge <- DGEList(counts = counts_df[, rownames(col_data), drop = FALSE])
      dge <- calcNormFactors(dge, method = "TMM")
      as_tibble(cpm(dge, normalized.lib.sizes = TRUE), rownames = "gene")
    }
  ),
  
  #  PCA plots for both methods (zero‐variance genes filtered out)
  tar_target(
    pca_plots,
    {
      make_pca <- function(df, label) {
        # Build matrix: rows = samples, cols = genes
        mat0 <- df %>%
          column_to_rownames("gene") %>%
          t()
        # Remove any zero‐variance genes (columns)
        vars <- apply(mat0, 2, var, na.rm = TRUE)
        mat  <- mat0[, vars > 0, drop = FALSE]
        
        # Run PCA with scaling
        pca    <- prcomp(mat, scale. = TRUE)
        scores <- as_tibble(pca$x, rownames = "sample")
        scores <- left_join(
          scores,
          col_data %>% rownames_to_column("sample"),
          by = "sample"
        )
        
        # Plot
        p <- ggplot(scores, aes(PC1, PC2, color = condition)) +
          geom_point(size = 3) +
          ggtitle(glue("PCA: {label}")) +
          theme_minimal()
        ggsave(glue("results/PCA_{label}.png"), plot = p)
        p
      }
      
      list(
        Median = make_pca(norm_median, "Median"),
        TMM    = make_pca(norm_tmm,    "TMM")
      )
    }
  ),
  
  #  Combined boxplots of distributions
  tar_target(
    boxplot_norm,
    {
      df_all <- bind_rows(
        norm_median %>% mutate(method = "Median"),
        norm_tmm    %>% mutate(method = "TMM")
      ) %>%
        pivot_longer(
          cols      = -c(gene, method),
          names_to  = "sample",
          values_to = "value"
        )
      p <- ggplot(df_all, aes(x = sample, y = value)) +
        geom_boxplot() +
        facet_wrap(~ method, scales = "free") +
        theme(
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        ggtitle("Expression Distributions by Normalization")
      ggsave("results/boxplot_norm.png", plot = p, width = 12, height = 5)
      p
    }
  )
  
) 

# Execute it
pipeline
