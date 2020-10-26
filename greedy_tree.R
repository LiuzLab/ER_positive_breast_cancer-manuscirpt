library(ggplot2)
library(ggtree)
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

ER_expression <- read_csv("./data/ER_quantificationtable.csv")
ER_expression$`20` <- NA
ER_expression <- ER_expression[c(1:12,20)]

ER_expression <- as.matrix(ER_expression)
er_mean <- mean(ER_expression,na.rm=TRUE)
ER_expression <- colMeans(ER_expression,na.rm=TRUE)
ER_expression[13] <- er_mean
ER_expression <- ER_expression/max(ER_expression)
names(ER_expression) <- c(1:12,20)

names(ER_expression) <- colnames(m_all)
cl1 <- c("s7","s1","s4","s5")
m_cl1 <- m_all[,cl1]
m_cl1 <- m_cl1 != 0
m_cl1 <- t(m_cl1[rowSums(m_cl1)>0,])
rownames(m_cl1) <- cl1

k <- nrow(m_cl1)
score_matrix <- matrix(nrow=k,ncol=k,data=0)
for (i in 1:k) {
  for (j in 1:k) {
    v <- sum(m_cl1[i,]*m_cl1[j,]) / sum((m_cl1[i,]+m_cl1[j,])>0)
    score_matrix[i,j] <- v}
}
colnames(score_matrix) <- cl1
rownames(score_matrix) <- cl1

colors <- c("white","#f4b5b4")
names(colors) <- c(0,1)
png("./figures/igor_MCF7/heatmap_cl1.png",width=1000,height=1000)
Heatmap(1*m_cl1, col = colors, cluster_rows = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        show_column_dend = FALSE, show_heatmap_legend = FALSE,
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()
ER <- ER_expression[cl1]
ggplot(data.frame(x=0,y=factor(names(ER),levels=rev(cl1)),er=ER)) + 
  geom_point(aes(x=x,y=y,size=ER))

cl2 <-  c("s2","s20","s12","s3","s9")
m_cl2 <- m_all[,cl2]
m_cl2 <- m_cl2 != 0
m_cl2 <- t(m_cl2[rowSums(m_cl2)>0,])
rownames(m_cl2) <- cl2

k <- nrow(m_cl2)
score_matrix <- matrix(nrow=k,ncol=k,data=0)
for (i in 1:k) {
  for (j in 1:k) {
    v <- sum(m_cl2[i,]*m_cl2[j,]) / sum((m_cl2[i,]+m_cl2[j,])>0)
    score_matrix[i,j] <- v}
}
colnames(score_matrix) <- cl2
rownames(score_matrix) <- cl2

colors <- c("white","#175110")
names(colors) <- c(0,1)
png("./figures/igor_MCF7/heatmap_cl2.png",width=1000,height=1000)
Heatmap(1*m_cl2, col = colors, cluster_rows = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        show_column_dend = FALSE, show_heatmap_legend = FALSE,
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()
ER <- ER_expression[cl2]

cl3 <-  c("s6","s11","s10","s8")
m_cl3 <- m_all[,cl3]
m_cl3 <- m_cl3 != 0
m_cl3 <- t(m_cl3[rowSums(m_cl3)>0,])
rownames(m_cl3) <- cl3

k <- nrow(m_cl3)
score_matrix <- matrix(nrow=k,ncol=k,data=0)
for (i in 1:k) {
  for (j in 1:k) {
    v <- sum(m_cl3[i,]*m_cl3[j,]) / sum((m_cl3[i,]+m_cl3[j,])>0)
    score_matrix[i,j] <- v}
}
colnames(score_matrix) <- cl3
rownames(score_matrix) <- cl3

colors <- c("white","#321051")
names(colors) <- c(0,1)
png("./figures/igor_MCF7/heatmap_cl3.png",width=1000,height=1000)
Heatmap(1*m_cl3, col = colors, cluster_rows = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        show_column_dend = FALSE, show_heatmap_legend = FALSE,
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()
ER <- ER_expression[cl3]
