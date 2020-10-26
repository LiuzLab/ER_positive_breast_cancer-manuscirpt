library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(huge)
library(aricode)
library(superheat)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(TraceQC)

for (i in 1:20) {
  print(i)
  tmp <- readRDS(sprintf("./cache/traceQC_obj/LIBILCM_%s.rds",i))
  svg(sprintf("./figures/igor_MCF7/circular/LCM_%s_deletion.svg",i))
  plot_deletion_hotspot(tmp)
  dev.off()

  svg(sprintf("./figures/igor_MCF7/circular/LCM_%s_insertion.svg",i))
  plot_insertion_hotspot(tmp)
  dev.off()
  p <- mutation_type(tmp)
  ggsave(sprintf("./figures/igor_MCF7/circular/LCM_%s_mutation_type.svg",i))

}

igor_mcf7 <- readRDS("./cache/heatmap/igor_mcf7.rds")

df_all <- data.frame()
for (i in c(1:12,20)) {
  s1 <- igor_mcf7[[i]] %>%
    filter(type=="deletion") %>%
    # filter(type != "unmutated") %>%
    mutate(count=count*1e6/sum(count)) %>%
    select(cigar,type,count)
  if (nrow(df_all) == 0) {
    df_all <- rbind(df_all,s1)
  }
  else {
    df_all <- full_join(df_all,s1) %>%
      mutate(count=replace_na(count,0))
  }
  df_all[[sprintf("s%d",i)]] <- df_all$count
  df_all$count <- NULL
}
types <- df_all$type
df_all$type <- NULL
df_all <- column_to_rownames(df_all,"cigar")
m_all <- as.matrix(df_all)
m_all[is.na(m_all)] <- 0
m_all <- log10(m_all+1)

col_fun = colorRamp2(c(0, 6), c("white", "darkred"))
svg("./figures/igor_MCF7/allmutation_1-12.svg")
Heatmap(m_all,col = col_fun,
        show_row_names = FALSE,
        row_split = types,
        name="log10 count")
dev.off()

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "#EEEEEE", "red"))
ER_expression <- read_csv("./data/ER_quantificationtable.csv")
ER_expression$`20` <- NA
ER_expression <- ER_expression[c(1:12,20)]

ER_expression <- as.matrix(ER_expression)
er_mean <- mean(ER_expression,na.rm=TRUE)
ER_expression <- colMeans(ER_expression,na.rm=TRUE)
ER_expression[13] <- er_mean
ER_expression <- ER_expression/max(ER_expression)
names(ER_expression) <- c(1:12,20)


m_all <- huge.npn(m_all)
huge_all <- huge(m_all,method="glasso")
g <- graph_from_adjacency_matrix(huge_all$path[[2]],mode="undirected")
# ceb <- cluster_edge_betweenness(g)
ceb <- cluster_walktrap(g)


g <- set_vertex_attr(g,"label", value = c(1:12,20))

V(g)$size <- ER_expression * 20
svg("./figures/igor_MCF7/MRF_all_1-12.svg")
plot(ceb,g, vertex.label.color="black",
     edge.color="grey")
dev.off()
# 
# mst_tree <- mst(g)
# svg("./figures/igor_MCF7/mst_arr1.svg")
# plot(mst_tree)
# dev.off()
# 
# svg("./figures/igor_MCF7/mst_arr2.svg")
# plot(mst_tree,layout = layout.reingold.tilford(mst_tree, root=13))
# dev.off()
###############################

df_all <- data.frame()
for (i in 1:20) {
  s1 <- igor_mcf7[[i]] %>%
    # filter(type=="deletion") %>%
    filter(type != "unmutated") %>%
    mutate(count=count*1e6/sum(count)) %>%
    select(cigar,type,count)
  if (nrow(df_all) == 0) {
    df_all <- rbind(df_all,s1)
  }
  else {
    df_all <- full_join(df_all,s1) %>%
      mutate(count=replace_na(count,0))
  }
  df_all[[sprintf("s%d",i)]] <- df_all$count
  df_all$count <- NULL
}
types <- df_all$type
df_all$type <- NULL
df_all <- column_to_rownames(df_all,"cigar")
m_all <- as.matrix(df_all)
m_all[is.na(m_all)] <- 0
m_all <- log10(m_all+1)

col_fun = colorRamp2(c(0, 6), c("white", "darkred"))
svg("./figures/igor_MCF7/allmutation_1-12.svg")
Heatmap(m_all,col = col_fun,
        show_row_names = FALSE,
        row_split = types,
        name="log10 count")
dev.off()

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "#EEEEEE", "red"))
ER_expression <- read_csv("./data/ER_quantificationtable.csv")
ER_expression$`20` <- NA
# ER_expression <- ER_expression[c(1:12,20)]

ER_expression <- as.matrix(ER_expression)
er_mean <- mean(ER_expression,na.rm=TRUE)
ER_expression <- colMeans(ER_expression,na.rm=TRUE)
ER_expression[20] <- er_mean
ER_expression <- ER_expression/max(ER_expression)
# names(ER_expression) <- c(1:12,20)


m_all <- huge.npn(m_all)
huge_all <- huge(m_all,method="glasso")
g <- graph_from_adjacency_matrix(huge_all$path[[2]],mode="undirected")
# ceb <- cluster_edge_betweenness(g)
ceb <- cluster_walktrap(g)


g <- set_vertex_attr(g,"label", value = c(1:12,20))

V(g)$size <- ER_expression * 20
svg("./figures/igor_MCF7/MRF_all_1-12.svg")
plot(ceb,g, vertex.label.color="black",
     edge.color="grey")
dev.off()
