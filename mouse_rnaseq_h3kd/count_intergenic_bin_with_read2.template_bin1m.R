library(data.table)
library(tidyverse)
bin_width <- "1m"

m2_rep1_name <- "MII #1"
m2_rep1_path <- sprintf("./NoRibs_MII201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
m2_rep1_total_read_count_main <- fread("./NoRibs_MII201209/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
m2_rep1_intergenic_read_count_main <- fread(sprintf("./NoRibs_MII201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
m2_rep2_name <- "MII #2"
m2_rep2_path <- sprintf("./NoRibs_MII201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
m2_rep2_total_read_count_main <- fread("./NoRibs_MII201218/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
m2_rep2_intergenic_read_count_main <- fread(sprintf("./NoRibs_MII201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
m2_mean_name <- "MII"

twocell18h_rep1_name <- "2C 18hpi #1"
twocell18h_rep1_path <- sprintf("./NoRibs_2c_18hpi201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell18h_rep1_total_read_count_main <- fread("./NoRibs_2c_18hpi201209/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell18h_rep1_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_18hpi201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell18h_rep2_name <- "2C 18hpi #2"
twocell18h_rep2_path <- sprintf("./NoRibs_2c_18hpi201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell18h_rep2_total_read_count_main <- fread("./NoRibs_2c_18hpi201218/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell18h_rep2_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_18hpi201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell18h_mean_name <- "2C 18hpi"

twocell_noinj_rep1_name <- "2C NoInjection #1"
twocell_noinj_rep1_path <- sprintf("./NoRibs_2c_Noinjection201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_noinj_rep1_total_read_count_main <- fread("./NoRibs_2c_Noinjection201209/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_noinj_rep1_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_Noinjection201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_noinj_rep2_name <- "2C NoInjection #2"
twocell_noinj_rep2_path <- sprintf("./NoRibs_2c_Noinjection201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_noinj_rep2_total_read_count_main <- fread("./NoRibs_2c_Noinjection201218/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_noinj_rep2_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_Noinjection201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_noinj_mean_name <- "2C NoInjection"

twocell_control_rep1_name <- "2C Control #1"
twocell_control_rep1_path <- sprintf("./NoRibs_2c_Control201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_control_rep1_total_read_count_main <- fread("./NoRibs_2c_Control201209/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_control_rep1_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_Control201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_control_rep2_name <- "2C Control #2"
twocell_control_rep2_path <- sprintf("./NoRibs_2c_Control201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_control_rep2_total_read_count_main <- fread("./NoRibs_2c_Control201218/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_control_rep2_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_Control201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_control_mean_name <- "2C Control"

twocell_h3kd_rep1_name <- "2C H3.1/3.2 KD #1"
twocell_h3kd_rep1_path <- sprintf("./NoRibs_2c_H3_1_3_2KD201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_h3kd_rep1_total_read_count_main <- fread("./NoRibs_2c_H3_1_3_2KD201209/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_h3kd_rep1_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_H3_1_3_2KD201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_h3kd_rep2_name <- "2C H3.1/3.2 KD #2"
twocell_h3kd_rep2_path <- sprintf("./NoRibs_2c_H3_1_3_2KD201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_h3kd_rep2_total_read_count_main <- fread("./NoRibs_2c_H3_1_3_2KD201218/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_h3kd_rep2_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_H3_1_3_2KD201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_h3kd_mean_name <- "2C H3.1/3.2 KD"

twocell_dnaminus_rep1_name <- "2C DNA(-) #1"
twocell_dnaminus_rep1_path <- sprintf("./NoRibs_2cAphidicolin201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_dnaminus_rep1_total_read_count_main <- fread("./NoRibs_2cAphidicolin201209/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_dnaminus_rep1_intergenic_read_count_main <- fread(sprintf("./NoRibs_2cAphidicolin201209/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_dnaminus_rep2_name <- "2C DNA(-) #2"
twocell_dnaminus_rep2_path <- sprintf("./NoRibs_2c_Aphidicoilin201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin%s.bed", bin_width)
twocell_dnaminus_rep2_total_read_count_main <- fread("./NoRibs_2c_Aphidicoilin201218/Aligned.sortedByCoord.out.uniq_collated.template.bed.count_main", header = FALSE)[[1,1]]
twocell_dnaminus_rep2_intergenic_read_count_main <- fread(sprintf("./NoRibs_2c_Aphidicoilin201218/Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin%s.bed.count_main", bin_width), header = FALSE)[[1,1]]
twocell_dnaminus_mean_name <- "2C DNA(-)"

header <- c("chr", "start", "end", "read_count_per_bin")
m2_rep1_data <- fread(m2_rep1_path, col.names = header)
m2_rep2_data <- fread(m2_rep2_path, col.names = header)
twocell18h_rep1_data <- fread(twocell18h_rep1_path, col.names = header)
twocell18h_rep2_data <- fread(twocell18h_rep2_path, col.names = header)
twocell_noinj_rep1_data <- fread(twocell_noinj_rep1_path, col.names = header)
twocell_noinj_rep2_data <- fread(twocell_noinj_rep2_path, col.names = header)
twocell_control_rep1_data <- fread(twocell_control_rep1_path, col.names = header)
twocell_control_rep2_data <- fread(twocell_control_rep2_path, col.names = header)
twocell_h3kd_rep1_data <- fread(twocell_h3kd_rep1_path, col.names = header)
twocell_h3kd_rep2_data <- fread(twocell_h3kd_rep2_path, col.names = header)
twocell_dnaminus_rep1_data <- fread(twocell_dnaminus_rep1_path, col.names = header)
twocell_dnaminus_rep2_data <- fread(twocell_dnaminus_rep2_path, col.names = header)

level_order <- c(m2_rep1_name, m2_rep2_name, twocell18h_rep1_name, twocell18h_rep2_name, twocell_noinj_rep1_name, twocell_noinj_rep2_name, twocell_control_rep1_name, twocell_control_rep2_name, twocell_h3kd_rep1_name, twocell_h3kd_rep2_name, twocell_dnaminus_rep1_name, twocell_dnaminus_rep2_name)
level_order <- factor(level_order, levels = level_order)
level_order2 <- c(m2_mean_name, twocell18h_mean_name, twocell_noinj_mean_name, twocell_control_mean_name, twocell_h3kd_mean_name, twocell_dnaminus_mean_name)
level_order2 <- factor(level_order2, levels = level_order2)

#total_count <- 1473376
#l1_count <- l1_data[, .N]
#l2_count <- l2_data[, .N]
#l3_count <- l3_data[, .N]
#l5_count <- l5_data[, .N]
#g1_data <- data.table(name = factor(level_order, levels = level_order),
#    count = c(l2_count, l1_count, l3_count, l5_count))
#g1_data$name <- factor(g1_data$name, levels = level_order)
#g1 <- ggplot(g1_data) + geom_col(aes(name, count))
#print(g1)

# check if the filter is correct
message("[INFO] filtered chr:")
message(paste0(m2_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), unique(chr)], collapse = "\t"))

# total count of the intergenic bins in the main chromosomes (1, 2, ..., X, Y, and MT; no PATCH chromosomes)
total_count_main <- 1467996

# main chromosome data
m2_rep1_data_main <- m2_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
m2_rep2_data_main <- m2_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell18h_rep1_data_main <- twocell18h_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell18h_rep2_data_main <- twocell18h_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_noinj_rep1_data_main <- twocell_noinj_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_noinj_rep2_data_main <- twocell_noinj_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_control_rep1_data_main <- twocell_control_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_control_rep2_data_main <- twocell_control_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_h3kd_rep1_data_main <- twocell_h3kd_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_h3kd_rep2_data_main <- twocell_h3kd_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_dnaminus_rep1_data_main <- twocell_dnaminus_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]
twocell_dnaminus_rep2_data_main <- twocell_dnaminus_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE)]

# intergenic bin count
m2_rep1_count_main <- m2_rep1_data_main[read_count_per_bin > 0, .N]
m2_rep2_count_main <- m2_rep2_data_main[read_count_per_bin > 0, .N]
twocell18h_rep1_count_main <- twocell18h_rep1_data_main[read_count_per_bin > 0, .N]
twocell18h_rep2_count_main <- twocell18h_rep2_data_main[read_count_per_bin > 0, .N]
twocell_noinj_rep1_count_main <- twocell_noinj_rep1_data_main[read_count_per_bin > 0, .N]
twocell_noinj_rep2_count_main <- twocell_noinj_rep2_data_main[read_count_per_bin > 0, .N]
twocell_control_rep1_count_main <- twocell_control_rep1_data_main[read_count_per_bin > 0, .N]
twocell_control_rep2_count_main <- twocell_control_rep2_data_main[read_count_per_bin > 0, .N]
twocell_h3kd_rep1_count_main <- twocell_h3kd_rep1_data_main[read_count_per_bin > 0, .N]
twocell_h3kd_rep2_count_main <- twocell_h3kd_rep2_data_main[read_count_per_bin > 0, .N]
twocell_dnaminus_rep1_count_main <- twocell_dnaminus_rep1_data_main[read_count_per_bin > 0, .N]
twocell_dnaminus_rep2_count_main <- twocell_dnaminus_rep2_data_main[read_count_per_bin > 0, .N]

g2_data <- data.table(sample = level_order,
    bin_count = c(m2_rep1_count_main, m2_rep2_count_main,
        twocell18h_rep1_count_main, twocell18h_rep2_count_main,
        twocell_noinj_rep1_count_main, twocell_noinj_rep2_count_main,
        twocell_control_rep1_count_main, twocell_control_rep2_count_main,
        twocell_h3kd_rep1_count_main, twocell_h3kd_rep2_count_main,
        twocell_dnaminus_rep1_count_main, twocell_dnaminus_rep2_count_main),
    read_count = c(m2_rep1_intergenic_read_count_main, m2_rep2_intergenic_read_count_main,
        twocell18h_rep1_intergenic_read_count_main, twocell18h_rep2_intergenic_read_count_main,
        twocell_noinj_rep1_intergenic_read_count_main, twocell_noinj_rep2_intergenic_read_count_main,
        twocell_control_rep1_intergenic_read_count_main, twocell_control_rep2_intergenic_read_count_main,
        twocell_h3kd_rep1_intergenic_read_count_main, twocell_h3kd_rep2_intergenic_read_count_main,
        twocell_dnaminus_rep1_intergenic_read_count_main, twocell_dnaminus_rep2_intergenic_read_count_main),
    total_read_count = c(m2_rep1_total_read_count_main, m2_rep2_total_read_count_main,
        twocell18h_rep1_total_read_count_main, twocell18h_rep2_total_read_count_main,
        twocell_noinj_rep1_total_read_count_main, twocell_noinj_rep2_total_read_count_main,
        twocell_control_rep1_total_read_count_main, twocell_control_rep2_total_read_count_main,
        twocell_h3kd_rep1_total_read_count_main, twocell_h3kd_rep2_total_read_count_main,
        twocell_dnaminus_rep1_total_read_count_main, twocell_dnaminus_rep2_total_read_count_main),
    category = rep(level_order2, each = 2))

# Mean count data
g2_mean_data <- g2_data[, .(bin_count = mean(bin_count), read_count = mean(read_count), total_read_count = mean(total_read_count)), by = .(category)]

make_plots <- function(yfunc, ylab_text) {
    yfunc <- enquo(yfunc)
    data_separate <- g2_data
    data_mean <- g2_mean_data
    g_common <- list(scale_y_continuous(labels = scales::comma),
        theme(axis.text.x = element_text(angle = 90)),
        coord_cartesian(ylim = c(0, NA)),
        ylab(ylab_text))
    g_separate <- ggplot(data_separate) + geom_col(aes(sample, !!yfunc, fill = category)) + g_common
    g_mean <- ggplot(data_mean) + geom_col(aes(category, !!yfunc)) + g_common
    g_mean2 <- ggplot(data_separate, aes(category, !!yfunc)) + geom_point() + geom_crossbar(data = data_mean, ymax = NA, ymin = NA) + g_common
    return(list(g_separate, g_mean, g_mean2))
}

g2_list <- make_plots(bin_count, "Intergenic bin count")
g3_list <- make_plots(total_read_count, "Total mapped read count")
g4_list <- make_plots(bin_count / total_read_count, "Intergenic bin count / Total mapped read count")
g5_list <- make_plots(read_count, "Intergenic read count")
g6_list <- make_plots(read_count / total_read_count, "Intergenic read count / Total mapped read count")
g7_list <- make_plots(read_count / bin_count, "Intergenic read count / Intergenic bin count")

# read count per bin
make_plots_per_bin <- function(plotdata) {
    g <- ggplot(plotdata) + geom_histogram(aes(y = read_count_per_bin / total_read_count), bins = 50) +
        facet_grid(. ~ sample) +
        ylab("Intergenic read count in a bin / Total read count") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90), panel.grid.minor = element_blank()) +
        scale_x_continuous(labels = scales::comma)
    g_xlog <- g + scale_x_continuous(trans = "log1p", breaks = 10 ** c(-Inf, 1:9), labels = scales::label_number_si())
    g_xzoom <- g + coord_cartesian(xlim = c(0, 1000))
    g_xzoom2 <- g + coord_cartesian(xlim = c(0, 100))
    return(list(g, g_xlog, g_xzoom, g_xzoom2))
}

data_main_list <- list(m2_rep1_data_main, m2_rep2_data_main,
        twocell18h_rep1_data_main, twocell18h_rep2_data_main,
        twocell_noinj_rep1_data_main, twocell_noinj_rep2_data_main,
        twocell_control_rep1_data_main, twocell_control_rep2_data_main,
        twocell_h3kd_rep1_data_main, twocell_h3kd_rep2_data_main,
        twocell_dnaminus_rep1_data_main, twocell_dnaminus_rep2_data_main)
mapply(function(data1, name1){data1[, "total_read_count" := .(g2_data[sample == name1, total_read_count])]}, data_main_list, level_order)
names(data_main_list) <- level_order
data_main_merged <- rbindlist(data_main_list, idcol = "sample")
data_main_merged$sample <- factor(data_main_merged$sample, levels = level_order)
g8_list <- make_plots_per_bin(data_main_merged)

pdf_path <- sprintf("count_intergenic_bin_with_read2.template_bin%s.pdf", bin_width)
cairo_pdf(pdf_path, height = 5, width = 5, onefile = TRUE)
print(g2_list[[1]])
print(g2_list[[2]])
print(g2_list[[3]])
print(g3_list[[1]])
print(g3_list[[2]])
print(g3_list[[3]])
print(g4_list[[1]])
print(g4_list[[2]])
print(g4_list[[3]])
print(g5_list[[1]])
print(g5_list[[2]])
print(g5_list[[3]])
print(g6_list[[1]])
print(g6_list[[2]])
print(g6_list[[3]])
print(g7_list[[1]])
print(g7_list[[2]])
print(g7_list[[3]])
invisible(dev.off())

pdf_path2 <- sprintf("count_intergenic_bin_with_read2.template_bin%s.per_bin.pdf", bin_width)
cairo_pdf(pdf_path2, height = 4, width = 14, onefile = TRUE)
print(g8_list[[1]])
print(g8_list[[2]])
print(g8_list[[3]])
print(g8_list[[4]])
invisible(dev.off())

quit()
library(RVAideMemoire)
table_main <- matrix(c(g2_data$bin_count, total_count_main - g2_data$bin_count), nrow = 2, byrow = TRUE,
    dimnames = list(c("WithRead", "WithoutRead"), level_order))
fisher.multcomp(table_main, p.method = "bonferroni")
