#' rnbtn_intersect_effects combines results from the total and unique counts and produces intersected .csv files
#' with specified fdr cutoffs.It also produces intersection plots

#' @param TC_files   : A list of Total counts files from rnbtn_fdr_effects along with full path.
#' Each file should have columns of fdr,log2_coefficient,locus_tag,mean,controlmean
#' @param UC_files   : A list of Unique counts files from rnbtn_fdr_effects along with full path.
#' Each file should have fdr,log2_coefficient,locus_tag,mean,controlmean
#' @param path: A path to write output files.
#' @param fdrcutoff: A cutoff for significant locus_tags.Default is 0.05


#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import stringr


#' @export

rnbtn_intersect_effects <- function(TC_files, UC_files, path, fdrcutoff = 0.05) {

    # Required packages tidyverse,reshape and doParallel

    ## validate input TC_files and UC_files
    if (length(TC_files) != length(UC_files)) {
        stop("The Total counts and Unique counts files have different lengths.Make sure that number of effect files are same for both")
    } else {
        cat("Running the intersection")
    }

    for (i in seq_len(TC_files)) {
        tryCatch({

            condition <- tail(stringr::str_split(TC_files, "/|.csv")[[i]], n = 2)

            # Total data frame
            tdf <- read.table(TC_files[i], header = T, sep = ",")
            if (all(c("fdr", "log2_coefficient", "locus_tag", "mean", "controlmean") %in%
                names(tdf)) == TRUE) {
                tdf <- tdf %>%
                  dplyr::select(fdr, log2_coefficient, locus_tag, mean, controlmean)
            } else {
                cat("The fdr results input of TC does not have some of fdr,log2_coefficient,locus_tag,mean,controlmean columns ",
                  "\n")
            }
            # Unique data frame

            udf <- read.table(UC_files[i], header = T, sep = ",")
            if (all(c("fdr", "log2_coefficient", "locus_tag", "mean", "controlmean") %in%
                names(udf)) == TRUE) {
                udf <- udf %>%
                  dplyr::select(fdr, log2_coefficient, locus_tag, mean, controlmean)
            } else {
                cat("The fdr results input of UC does not have some of fdr,log2_coefficient,locus_tag,mean,controlmean columns ",
                  "\n")
            }

            # merge them by using locus_tag
            df <- na.omit(merge(tdf, udf, by = "locus_tag", all = TRUE))
            colnames(df) <- c("locus_tag", "Total.fdr", "Total.coefficient", "Total.mean",
                "Total.Control.Mean", "Unique.fdr", "Unique.coefficient", "Unique.mean",
                "Unique.Control.Mean")

            # seperate out uniquely significant
            dfplot1 <- df %>%
                dplyr::mutate(fdr = case_when(Unique.fdr < fdrcutoff ~ "Unique_sig"))
            # seperate out total significant
            dfplot2 <- df %>%
                dplyr::mutate(fdr = case_when(Total.fdr < fdrcutoff ~ "Total_sig"))
            # seperate out non significant
            dfplot3 <- df %>%
                dplyr::mutate(fdr = case_when(Total.fdr > fdrcutoff & Unique.fdr > fdrcutoff ~
                  "ns"))
            # combine
            dfplot <- rbind(dfplot1, dfplot2, dfplot3) %>%
                tidyr::drop_na()
            # plot data frame
            dfplot <- dfplot %>%
                dplyr::mutate(pointsize = case_when(fdr == "ns" | fdr == "ns" ~ "1", fdr ==
                  "Total_sig" | fdr == "Unique_sig" ~ "1"))
            # ggplot
            suppressWarnings(p <- ggplot2::ggplot(dfplot, aes(Total.coefficient, Unique.coefficient)) +
                geom_point(aes(shape = fdr, color = fdr, size = pointsize)) + xlab("Total Counts Coefficient") +
                ylab("Unique Counts Coefficient") + theme_bw() + theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank()) +
                ggtitle(condition) + scale_shape_manual(labels = c("Not Significant",
                "Unique Significant", "Total Significant"), values = c(ns = 16, Unique_sig = 6,
                Total_sig = 2)) + scale_color_manual(labels = c("Not Significant",
                "Unique Significant", "Total Significant"), values = c(ns = "azure3",
                Unique_sig = "darkblue", Total_sig = "darkred")) + guides(size = FALSE) +
                xlim(-5, 5) + ylim(-5, 5) + theme(legend.title = element_blank()) +
                theme(legend.justification = "top") + theme(plot.title = element_text(size = 24)) +
                theme(legend.text = element_text(size = 18)) + theme(axis.text = element_text(size = 18),
                axis.title.y = element_text(size = 1, angle = 90)) + theme(aspect.ratio = 1) +
                theme(panel.border = element_blank(), axis.line = element_line(color = "black")))
            # Write csv files of total sig,unique sig and both to the paths
            write_csv(dfplot2 %>%
                dplyr::filter(fdr == "Total_sig"), (file.path(path, paste0("TC_", condition,
                "_sig.csv"))[1]))
            write_csv(dfplot1 %>%
                dplyr::filter(fdr == "Unique_sig"), (file.path(path, paste0("UC_", condition,
                "_sig.csv"))[1]))
            dfb <- arrange(df %>%
                dplyr::filter(Unique.fdr < fdrcutoff & Total.fdr < fdrcutoff), Total.fdr,
                Unique.fdr)
            write_csv(dfb, (file.path(path, paste0("TC_", condition, "_sig.csv"))[1]))
            # write plots
            suppressWarnings(ggpubr::ggsave(filename = (file.path(path, paste0(condition,
                "_figure.pdf"))[1]), device = "pdf", plot = p, width = 7, height = 7,
                dpi = 300, units = "in"))
        }, error = function(e) {
        })
    }
}
