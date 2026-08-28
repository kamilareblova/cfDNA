library(tidyr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

input_file <- args[1] # merged soubor
list_file <- args[2] # seznam genu
output_file <- args[3]

# nacteni seznamu genu pro filtrovani
selected_genes <- readLines(list_file)
selected_genes <- trimws(selected_genes)
selected_genes <- unique(selected_genes[selected_genes != ""])

print("filter genes from merged files")

print(paste("Selected genes:", paste(selected_genes, collapse = ", ")))

# def fce filtr
filter_df <- function(input_file, selected_genes) {
    print(paste("Processing file:", input_file))

    # nacteni dat z merged souboru
    df_data <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", fill = TRUE)

    # prejmenovani duplicitnich sloupcu -> CHROM.1	POS.1	REF.1	ALT.1
    names(df_data) <- make.unique(names(df_data))

    # dup_names <- names(df_data)[duplicated(names(df_data))]
    # print(dup_names)

    # porovnani sloupcu s hledanymi geny
    if (!"Gene.refGeneWithVer" %in% names(df_data)) {stop("Column 'Gene.refGeneWithVer' was not found in input file.")}

    # prepis \x3b na ;
    df_data$Gene.refGeneWithVer <- str_replace_all(as.character(df_data$Gene.refGeneWithVer), "\\\\x3b", ";")

    # vytvoreni regex patternu pro hledani genu
    pattern <- paste0("\\b(", paste(selected_genes, collapse = "|"), ")\\b")

    # filtrovani dat podle genu - hledani v sloupci Gene.refGeneWithVer
    filtered_df <- df_data %>%
        filter(str_detect(Gene.refGeneWithVer, pattern))

    return(filtered_df)
}

# pouziti funkce
filtered_df <- filter_df(input_file, selected_genes)

print(paste("Number of filtered rows:", nrow(filtered_df)))

write.table(filtered_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

print(paste("Saved filtered file to:", output_file))
