library(tidyverse)
library(Homo.sapiens)
library(xlsx)
hs <- Homo.sapiens

inter <- read_tsv("Interactions.txt")
filtered <- filter(inter, Left.Component.size > 20)
splitted <- split(filtered, filtered$Left.Component.Id)
ensg <- map(splitted, ~ unique(.x$Right.ENSEMBL, .x$Left.ENSEMBL))

rna <- read_tsv("GSE67295_GeneExpressionPolyA.txt")
colnames(rna)[1] <- "id"

rna <- dplyr::select(rna, id,
                     "E2_1" = contains("Sample21"),
                     "E2_2" = contains("Sample22"),
                     "Veh_1" = contains("Sample29"),
                     "Veh_2" = contains("Sample30"))
                     
nm_2_ensg <- mapIds(hs, keys = rna$id, column = "ENSEMBL", keytype = "REFSEQ")
rna$ensg <- nm_2_ensg[rna$id]

# Output total
df <- imap(ensg, ~ tibble(community = .y, ensg = .x)) %>%
    map_df(~ .x) %>%
    left_join(rna, by = "ensg") %>%
    group_by(community) %>%
    summarize(fc_1 = sum(abs(E2_1 - Veh_1), na.rm = TRUE), fc_2 = sum(abs(E2_2 - Veh_2), na.rm = TRUE),
              sd_1 = sd(abs(E2_1 - Veh_1), na.rm = TRUE), sd_2 = sd(abs(E2_2 - Veh_2), na.rm = TRUE))

fc <- dplyr::select(df, community, fc_1, fc_2) %>% gather(rep, fc, - community) %>% mutate(rep = str_replace(rep, "fc_", ""))
sd <- dplyr::select(df, community, sd_1, sd_2) %>% gather(rep, sd, - community) %>% mutate(rep = str_replace(rep, "sd_", ""))
df <- left_join(fc, sd, by = c("community", "rep"))

pdf("GSE67295_output_total.pdf")
ggplot(df, aes(x = log2(fc), y = log2(sd), color = rep)) +
    geom_point() +
    theme_minimal() +
    xlab("log2(Total output FPKM)") +
    ylab("log2(Variance)") +
    ggtitle("Total output vs variance by community after E2 treatment")
dev.off()
write.xlsx2(as.data.frame(df), "GSE67295_output_total.xlsx", row.names = FALSE)

# New version
df <- imap(ensg, ~ tibble(community = .y, ensg = .x)) %>%
    map_df(~ .x) %>%
    filter(!is.na(ensg)) %>%
    left_join(rna, by = "ensg") %>%
    group_by(community) %>%
#    summarize(fc_1 = log2(mean(E2_1 / Veh_1, na.rm = TRUE)+1), fc_2 = log2(mean(E2_2 / Veh_2, na.rm = TRUE)+1),
#    summarize(fc_1 = log2(mean(E2_1 / (Veh_1+1), na.rm = TRUE)+1), fc_2 = log2(mean(E2_2 / (Veh_2+1), na.rm = TRUE)+1),
    summarize(fc_1 = log2(mean(E2_1 / (Veh_1+1), na.rm = TRUE)), fc_2 = log2(mean(E2_2 / (Veh_2+1), na.rm = TRUE)),
              sd_1 = log2(sd(E2_1 - Veh_1, na.rm = TRUE)), sd_2 = log2(sd(E2_2 - Veh_2, na.rm = TRUE)),
              tot_1 = log2(sum(E2_1, na.rm = TRUE)), tot_2 = log2(sum(E2_2, na.rm = TRUE)))

fc <- dplyr::select(df, community, fc_1, fc_2) %>% gather(rep, value, - community) %>% mutate(rep = str_replace(rep, "fc_", ""), type = "fold_change")
sd <- dplyr::select(df, community, sd_1, sd_2) %>% gather(rep, value, - community) %>% mutate(rep = str_replace(rep, "sd_", ""), type = "variance")
tot <- dplyr::select(df, community, tot_1, tot_2) %>% gather(rep, tot, - community) %>% mutate(rep = str_replace(rep, "tot_", ""))
fc <- fc %>% left_join(tot, by = c("community", "rep"))
sd <- sd %>% left_join(tot, by = c("community", "rep"))

pdf("GSE67295_output_vs_variance.pdf")
ggplot(fc, aes(x = tot, y = value, color = rep)) +
    geom_point() +
    theme_minimal() +
    geom_smooth() +
    xlab("log2(total output)") +
    ylab("log2(fold change)") +
    ggtitle("Total output vs fold change after E2 treatment")
ggplot(sd, aes(x = tot, y = value, color = rep)) +
    geom_point() +
    theme_minimal() +
    geom_smooth() +
    xlab("log2(total output)") +
    ylab("log2(variance)") +
    ggtitle("Total output vs variance after E2 treatment")
dev.off()

#df <- rbind(fc, sd) %>% left_join(tot, by = c("community", "rep"))

pdf("GSE67295_fold_change.pdf")
ggplot(df, aes(x = log2(fc), y = log2(sd), color = rep)) +
    geom_point() +
    theme_minimal() +
    xlab("log2(Fold Change)") +
    ylab("log2(Variance)") +
    ggtitle("Total output vs variance by community after E2 treatment")
dev.off()
write.xlsx2(as.data.frame(df), "GSE67295_fold_change.xlsx", row.names = FALSE)


