library(data.table)
library(tidyverse)
library(cowplot)
m2_rep1_name <- "MII #1"
m2_rep1_path <- "./NoRibs_MII201209/Aligned.sortedByCoord.out.uniq_collated.template.bed"
m2_rep2_name <- "MII #2"
m2_rep2_path <- "./NoRibs_MII201218/Aligned.sortedByCoord.out.uniq_collated.template.bed"
m2_mean_name <- "MII"

twocell18h_rep1_name <- "2C 18hpi #1"
twocell18h_rep1_path <- "./NoRibs_2c_18hpi201209/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell18h_rep2_name <- "2C 18hpi #2"
twocell18h_rep2_path <- "./NoRibs_2c_18hpi201218/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell18h_mean_name <- "2C 18hpi"

twocell_noinj_rep1_name <- "2C NoInjection #1"
twocell_noinj_rep1_path <- "./NoRibs_2c_Noinjection201209/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_noinj_rep2_name <- "2C NoInjection #2"
twocell_noinj_rep2_path <- "./NoRibs_2c_Noinjection201218/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_noinj_mean_name <- "2C NoInjection"

twocell_control_rep1_name <- "2C Control #1"
twocell_control_rep1_path <- "./NoRibs_2c_Control201209/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_control_rep2_name <- "2C Control #2"
twocell_control_rep2_path <- "./NoRibs_2c_Control201218/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_control_mean_name <- "2C Control"

twocell_h3kd_rep1_name <- "2C H3.1/3.2 KD #1"
twocell_h3kd_rep1_path <- "./NoRibs_2c_H3_1_3_2KD201209/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_h3kd_rep2_name <- "2C H3.1/3.2 KD #2"
twocell_h3kd_rep2_path <- "./NoRibs_2c_H3_1_3_2KD201218/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_h3kd_mean_name <- "2C H3.1/3.2 KD"

twocell_dnaminus_rep1_name <- "2C DNA(-) #1"
twocell_dnaminus_rep1_path <- "./NoRibs_2cAphidicolin201209/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_dnaminus_rep2_name <- "2C DNA(-) #2"
twocell_dnaminus_rep2_path <- "./NoRibs_2c_Aphidicoilin201218/Aligned.sortedByCoord.out.uniq_collated.template.bed"
twocell_dnaminus_mean_name <- "2C DNA(-)"

header <- c("chr", "start", "end")
m2_rep1_data <- fread(m2_rep1_path, col.names = header)[, .(chr, tlen = end - start)]
m2_rep2_data <- fread(m2_rep2_path, col.names = header)[, .(chr, tlen = end - start)]
twocell18h_rep1_data <- fread(twocell18h_rep1_path, col.names = header)[, .(chr, tlen = end - start)]
twocell18h_rep2_data <- fread(twocell18h_rep2_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_noinj_rep1_data <- fread(twocell_noinj_rep1_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_noinj_rep2_data <- fread(twocell_noinj_rep2_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_control_rep1_data <- fread(twocell_control_rep1_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_control_rep2_data <- fread(twocell_control_rep2_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_h3kd_rep1_data <- fread(twocell_h3kd_rep1_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_h3kd_rep2_data <- fread(twocell_h3kd_rep2_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_dnaminus_rep1_data <- fread(twocell_dnaminus_rep1_path, col.names = header)[, .(chr, tlen = end - start)]
twocell_dnaminus_rep2_data <- fread(twocell_dnaminus_rep2_path, col.names = header)[, .(chr, tlen = end - start)]

level_order <- c(m2_rep1_name, m2_rep2_name, twocell18h_rep1_name, twocell18h_rep2_name, twocell_noinj_rep1_name, twocell_noinj_rep2_name, twocell_control_rep1_name, twocell_control_rep2_name, twocell_h3kd_rep1_name, twocell_h3kd_rep2_name, twocell_dnaminus_rep1_name, twocell_dnaminus_rep2_name)
level_order <- factor(level_order, levels = level_order)
#level_order2 <- c(m2_mean_name, twocell18h_mean_name, twocell_noinj_mean_name, twocell_control_mean_name, twocell_h3kd_mean_name, twocell_dnaminus_mean_name)
#level_order2 <- factor(level_order2, levels = level_order2)
#level_order <- c(m2_rep1_name, m2_rep2_name, twocell18h_rep1_name, twocell18h_rep2_name)
#level_order <- factor(level_order, levels = level_order)

# check if the filter is correct
major_chr_pattern <- "^(?!JH|GL).*$"
#message("[INFO] filtered chr:")
#message(paste0(m2_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), unique(chr)], collapse = "\t"))

m2_rep1_data_main <- m2_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
m2_rep2_data_main <- m2_rep2_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell18h_rep1_data_main <- twocell18h_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell18h_rep2_data_main <- twocell18h_rep2_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_noinj_rep1_data_main <- twocell_noinj_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_noinj_rep2_data_main <- twocell_noinj_rep2_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_control_rep1_data_main <- twocell_control_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_control_rep2_data_main <- twocell_control_rep2_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_h3kd_rep1_data_main <- twocell_h3kd_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_h3kd_rep2_data_main <- twocell_h3kd_rep2_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_dnaminus_rep1_data_main <- twocell_dnaminus_rep1_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]
twocell_dnaminus_rep2_data_main <- twocell_dnaminus_rep2_data[grepl(major_chr_pattern, chr, perl = TRUE), .(tlen)]

data_main_list <- list(m2_rep1_data_main, m2_rep2_data_main, twocell18h_rep1_data_main, twocell18h_rep2_data_main, twocell_noinj_rep1_data_main, twocell_noinj_rep2_data_main, twocell_control_rep1_data_main, twocell_control_rep2_data_main, twocell_h3kd_rep1_data_main, twocell_h3kd_rep2_data_main, twocell_dnaminus_rep1_data_main, twocell_dnaminus_rep2_data_main)
#data_main_list <- list(m2_rep1_data_main, m2_rep2_data_main, twocell18h_rep1_data_main, twocell18h_rep2_data_main)
names(data_main_list) <- level_order
rm(m2_rep1_data, m2_rep2_data, twocell18h_rep1_data, twocell18h_rep2_data, twocell_noinj_rep1_data, twocell_noinj_rep2_data, twocell_control_rep1_data, twocell_control_rep2_data, twocell_h3kd_rep1_data, twocell_h3kd_rep2_data, twocell_dnaminus_rep1_data, twocell_dnaminus_rep2_data)
gc()

xmax <- max(sapply(data_main_list, function(d){ d[, max(tlen)] }))
binwidth1 <- 10000
xmax2 <- 10000
binwidth2 <- 100
make_plots <- function(data, title) {
    p1 <- ggplot(data) + geom_histogram(aes(tlen), binwidth = binwidth1, boundary = 0) +
        ggtitle(title) + xlab("Template length") +
        coord_cartesian(xlim = c(0, xmax + binwidth1)) +
        scale_y_continuous(labels = scales::comma)
    p1_yzoom <- p1 + coord_cartesian(ylim = c(0, 100))
    p1_ylog <- p1 + scale_y_continuous(labels = scales::comma, trans = "log10")
    p2 <- ggplot(data) + geom_histogram(aes(tlen), binwidth = binwidth2, boundary = 0) +
        ggtitle(title) + xlab("Template length") +
        coord_cartesian(xlim = c(0, xmax2 + binwidth2)) +
        scale_y_continuous(labels = scales::comma)
    p2_ylog <- p2 + scale_y_continuous(labels = scales::comma, trans = "log10")
    return(list("base" = p1, "base_yzoom" = p1_yzoom, "base_ylog" = p1_ylog, "xzoom" = p2, "xzoom_ylog" = p2_ylog))
}

data_main_plotlist <- mapply(make_plots, data_main_list, level_order, SIMPLIFY = FALSE)

pdf_path <- "plot_template_length.pdf"
pdf(pdf_path, height = 1.5 * length(data_main_plotlist), width = 3)
plot_grid(plotlist = lapply(data_main_plotlist, function(l){l[["base"]]}), ncol = 1, align = "v")
plot_grid(plotlist = lapply(data_main_plotlist, function(l){l[["base_yzoom"]]}), ncol = 1, align = "v")
plot_grid(plotlist = lapply(data_main_plotlist, function(l){l[["base_ylog"]]}), ncol = 1, align = "v")
plot_grid(plotlist = lapply(data_main_plotlist, function(l){l[["xzoom"]]}), ncol = 1, align = "v")
plot_grid(plotlist = lapply(data_main_plotlist, function(l){l[["xzoom_ylog"]]}), ncol = 1, align = "v")
invisible(dev.off())
