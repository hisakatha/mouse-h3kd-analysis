library(data.table)
library(ggplot2)
m2_rep1_name <- "MII #1"
m2_rep1_path <- "./NoRibs_MII201209/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
m2_rep2_name <- "MII #2"
m2_rep2_path <- "./NoRibs_MII201218/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
m2_mean_name <- "MII"

twocell18h_rep1_name <- "2C 18hpi #1"
twocell18h_rep1_path <- "./NoRibs_2c_18hpi201209/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell18h_rep2_name <- "2C 18hpi #2"
twocell18h_rep2_path <- "./NoRibs_2c_18hpi201218/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell18h_mean_name <- "2C 18hpi"

twocell_noinj_rep1_name <- "2C NoInjection #1"
twocell_noinj_rep1_path <- "./NoRibs_2c_Noinjection201209/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_noinj_rep2_name <- "2C NoInjection #2"
twocell_noinj_rep2_path <- "./NoRibs_2c_Noinjection201218/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_noinj_mean_name <- "2C NoInjection"

twocell_control_rep1_name <- "2C Control #1"
twocell_control_rep1_path <- "./NoRibs_2c_Control201209/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_control_rep2_name <- "2C Control #2"
twocell_control_rep2_path <- "./NoRibs_2c_Control201218/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_control_mean_name <- "2C Control"

twocell_h3kd_rep1_name <- "2C H3.1/3.2 KD #1"
twocell_h3kd_rep1_path <- "./NoRibs_2c_H3_1_3_2KD201209/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_h3kd_rep2_name <- "2C H3.1/3.2 KD #2"
twocell_h3kd_rep2_path <- "./NoRibs_2c_H3_1_3_2KD201218/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_h3kd_mean_name <- "2C H3.1/3.2 KD"

twocell_dnaminus_rep1_name <- "2C DNA(-) #1"
twocell_dnaminus_rep1_path <- "./NoRibs_2cAphidicolin201209/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_dnaminus_rep2_name <- "2C DNA(-) #2"
twocell_dnaminus_rep2_path <- "./NoRibs_2c_Aphidicoilin201218/Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"
twocell_dnaminus_mean_name <- "2C DNA(-)"

header <- c("chr", "start", "end")
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

# total count of the intergenic bins in the main chromosomes (no PATCH chromosomes)
total_count_main <- 1467996

m2_rep1_count_main <- m2_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
m2_rep2_count_main <- m2_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell18h_rep1_count_main <- twocell18h_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell18h_rep2_count_main <- twocell18h_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_noinj_rep1_count_main <- twocell_noinj_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_noinj_rep2_count_main <- twocell_noinj_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_control_rep1_count_main <- twocell_control_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_control_rep2_count_main <- twocell_control_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_h3kd_rep1_count_main <- twocell_h3kd_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_h3kd_rep2_count_main <- twocell_h3kd_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_dnaminus_rep1_count_main <- twocell_dnaminus_rep1_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]
twocell_dnaminus_rep2_count_main <- twocell_dnaminus_rep2_data[grepl("^(?!JH|GL).*$", chr, perl = TRUE), .N]

g2_data <- data.table(name = level_order,
    count = c(m2_rep1_count_main, m2_rep2_count_main,
        twocell18h_rep1_count_main, twocell18h_rep2_count_main,
        twocell_noinj_rep1_count_main, twocell_noinj_rep2_count_main,
        twocell_control_rep1_count_main, twocell_control_rep2_count_main,
        twocell_h3kd_rep1_count_main, twocell_h3kd_rep2_count_main,
        twocell_dnaminus_rep1_count_main, twocell_dnaminus_rep2_count_main),
    category = rep(level_order2, each = 2))
g2 <- ggplot(g2_data) + geom_col(aes(name, count, fill = category)) + scale_y_continuous(labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 90))

# Mean count data
g3_data <- data.table(name = level_order2,
    count = c((m2_rep1_count_main + m2_rep2_count_main) / 2,
        (twocell18h_rep1_count_main + twocell18h_rep2_count_main) / 2,
        (twocell_noinj_rep1_count_main + twocell_noinj_rep2_count_main) / 2,
        (twocell_control_rep1_count_main + twocell_control_rep2_count_main) / 2,
        (twocell_h3kd_rep1_count_main + twocell_h3kd_rep2_count_main) / 2,
        (twocell_dnaminus_rep1_count_main + twocell_dnaminus_rep2_count_main) / 2))
g3 <- ggplot(g3_data) + geom_col(aes(name, count)) + scale_y_continuous(labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 90)) + ylab("Mean count")

g4 <- ggplot(g2_data, aes(category, count)) + scale_y_continuous(labels = scales::comma) +
    geom_point() + theme(axis.text.x = element_text(angle = 90)) +
    stat_summary(fun = mean, geom = "crossbar")

pdf_path <- "count_intergenic_bin_with_read.raw.pdf"
pdf(pdf_path, height = 5, width = 5)
print(g2)
print(g3)
print(g4)
invisible(dev.off())

quit()
library(RVAideMemoire)
table_main <- matrix(c(g2_data$count, total_count_main - g2_data$count), nrow = 2, byrow = TRUE,
    dimnames = list(c("WithRead", "WithoutRead"), level_order))
fisher.multcomp(table_main, p.method = "bonferroni")
