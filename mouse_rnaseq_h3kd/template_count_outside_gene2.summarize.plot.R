library(data.table)
library(tidyverse)

counts <- fread("template_count_outside_gene2.summarize.csv")
output <- "template_count_outside_gene2.summarize.plot.pdf"

sample_names <- c("MII #1", "MII #2", "2C 18hpi #1", "2C 18hpi #2", "2C No Injection #1", "2C No Injection #2",
    "2C Control #1", "2C Control #2", "2C H3.1/3.2 KD #1", "2C H3.1/3.2 KD #2", "2C DNA(-) #1", "2C DNA(-) #2")
sample_names_level <- factor(sample_names, levels = sample_names)

sample_categories <- c("MII", "2C 18hpi", "2C No Injection",
    "2C Control", "2C H3.1/3.2 KD", "2C DNA(-)")
sample_categories_level <- factor(sample_categories, levels = sample_categories)
counts$sample_category <- factor(counts$sample_category, levels = sample_categories)
counts$sample_name <- factor(counts$sample_name, levels = sample_names)

counts <- counts[chr_set == "main"]

mean_counts <- counts[, .(mean_intergenic_read_ratio = mean(intergenic_read_count / total_read_count)), by=c("sample_category", "chr_set", "gene_extension_width")]
mean_counts <- mean_counts[chr_set == "main"]

common_settings <- list(scale_x_continuous(trans = "log10", labels = scales::label_number_si()),
    coord_cartesian(xlim = c(1, NA)),
    annotation_logticks(sides = "b"),
    xlab("Gene extension length"), ylab("Intergenic read ratio"), labs(color = "Category", shape = "Category"))

p1 <- ggplot() +
    geom_point(data = counts, aes(gene_extension_width, intergenic_read_count / total_read_count, color = sample_category, shape = sample_category), alpha = 0.3) +
    geom_point(data = mean_counts, aes(gene_extension_width, mean_intergenic_read_ratio, color = sample_category, shape = sample_category)) +
    geom_line(data = mean_counts, aes(gene_extension_width, mean_intergenic_read_ratio, color = sample_category)) +
    common_settings

p1_each <- ggplot() +
    geom_point(data = counts, aes(gene_extension_width, intergenic_read_count / total_read_count, color = sample_category, shape = sample_category)) +
    common_settings

p1_mean <- ggplot() +
    geom_point(data = mean_counts, aes(gene_extension_width, mean_intergenic_read_ratio, color = sample_category, shape = sample_category)) +
    geom_line(data = mean_counts, aes(gene_extension_width, mean_intergenic_read_ratio, color = sample_category)) +
    common_settings

widths <- counts[, sort(unique(gene_extension_width))]
common_settings2 <- list(scale_y_continuous(labels = scales::comma),
    theme(axis.text.x = element_text(angle = 90)),
    coord_cartesian(ylim = c(0, NA)),
    ylab("Intergenic read ratio"),
    theme(legend.justification = c(0.5, 1), panel.grid.minor = element_blank()))
p_each_width <- list()
for (i in widths) {
    p_each_width1 <- ggplot(counts[gene_extension_width == i]) +
        geom_col(aes(sample_name, intergenic_read_count / total_read_count, fill = sample_category)) +
        common_settings2 + xlab("Sample") + labs(fill = "Category") + ggtitle(sprintf("Gene extension length: %d", i))
    p_each_width2 <- ggplot() + geom_point(data = counts[gene_extension_width == i], aes(sample_category, intergenic_read_count / total_read_count)) +
        geom_crossbar(data = mean_counts[gene_extension_width == i], aes(sample_category, mean_intergenic_read_ratio), ymax = NA, ymin = NA) +
        common_settings2 + xlab("Category") + ggtitle(sprintf("Gene extension length: %d", i))
    p_each_width <- c(p_each_width, list(p_each_width1, p_each_width2))
}

cairo_pdf(output, width = 4, height = 5, family = "Arial", onefile = TRUE)
print(p1)
print(p1_each)
print(p1_mean)
print(p_each_width)
invisible(dev.off())
