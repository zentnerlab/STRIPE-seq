
library("tidyverse")
library("cowplot")

#################################
## STRIPE-seq vs. RNA-seq DEGS ##
#################################

## Load and prepare data.

degs <- list.files("references", pattern = ".*\\.tsv$", full.names = TRUE) %>%
	map(
		~ read.delim(., sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
			as_tibble(.name_repair = "unique") %>%
			mutate(change = case_when(
				log2FC <= -1 & FDR <= 0.05 ~ "downregulated",
				log2FC >= 1 & FDR <= 0.05 ~ "upregulated",
				TRUE ~ "unchaged"
			))
	)

## Merge data.

degs <- purrr::reduce(degs, full_join, by = "geneId", suffix = c("_RNAseq", "_STRIPEseq")) %>%
	replace_na(list(
		FDR_RNAseq = 1, FDR_STRIPEseq = 1,
		log2FC_RNAseq = 0, log2FC_STRIPEseq = 0,
		change_RNAseq = "unchaged", change_STRIPEseq = "unchaged"
	))

## Calculate cumulative average.

merged_de <- merged_de %>%
	dplyr::arrange(desc(abs(log2FC_RNAseq))) %>%
	mutate(
		match = ifelse(change_RNAseq == change_STRIPE, 1, 0),
		cumfrac = cummean(match)
	)

## Prepare the data for plotting.

plot_data <- merged_de %>%
	dplyr::select(geneId, log2FC_RNAseq, cumfrac) %>%
	rowid_to_column %>%
	gather(key = "metric", value = "value", -rowid, -geneId)

## Plot the data.

p <- ggplot(plot_data, aes(x = rowid, y = abs(value))) +
	theme_bw() +
	geom_line() +
	theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
		text = element_text(size = 12)
        ) +
	facet_grid(metric ~ ., scales = "free")

ggsave("cumulative_plot.pdf", plot = p, device = cairo_pdf, height = 3, width = 5)
