library(tidyverse)
annot <- read.table('~/Downloads/Mcavernosa.maker.coding.gff3', header = FALSE, sep = "\t", col.names = paste0("V",seq_len(9)), fill = TRUE, comment.char="", stringsAsFactors = FALSE)
index.rows <- as.numeric(rownames(annot[!startsWith(annot$V1, "S") & !startsWith(annot$V1, "x"),]))
gene.rows <- annot[startsWith(annot$V1, "S") | startsWith(annot$V1, "x"),]
genes <- unlist(lapply(strsplit(gene.rows$V9, ";"), `[[`, 1))
genes2 <- sub("ID=(Mcavernosa[0-9]+).*", "\\1", genes)
df.nospan <- mutate(gene.rows, genename = genes2) %>%
  group_by(genename) %>%
  mutate(keep = ifelse(length(unique(V1)) == 1, TRUE, FALSE)) %>%
  filter(keep == TRUE) %>%
  ungroup()

find.overlaps <- filter(df.nospan, V3 == 'gene') %>%
  group_by(V1, V7) %>%
  arrange(V4, .by_group = TRUE) %>%
  group_by(dups = cumsum(cummax(lag(V5, default = first(V5))) < V4))

df.whichoverlap <- find.overlaps %>%
  group_by(V1, V7, dups) %>%
  summarise(n = n()) %>%
  filter(n > 1)



df.overlaps <- right_join(find.overlaps, df.whichoverlap)

test <- group_by(df.overlaps, V1, V7, dups) %>%
  arrange(V4, .by_group = TRUE) %>%
  group_by(dups2 = cumsum(cummax(lag(V5, default = first(V5))) < V4))
  

group_by(V1, dups) %>%
  summarise(start = first(V4), stop = max(V5))

which.nospan <- as.numeric(rownames(gene.rows[gene.rows$V9 %in% df.nospan$V9,]))
allKeeps <- c(index.rows, which.nospan)
allKeeps.sort <- as.character(allKeeps[order(allKeeps)])

fileConn<-file("~/Downloads/rown")
writeLines(allKeeps.sort, fileConn)
close(fileConn)
