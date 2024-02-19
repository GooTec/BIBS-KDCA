library(tidyverse)
library(glmnet)
library(ggplot2)
library(gridExtra)

# Working directory 설정
dir <- "C:/Users/KDCA_forR/"
setwd(dir)

# 데이터 및 결과 파일 경로 설정
data_root <- "Data/shannon_"
save_root <- "Result/"

# 데이터 읽어오기
BCR1H <- read.csv(paste0(data_root, "BCR_01_HC.csv"), header = TRUE)
BCR1L <- read.csv(paste0(data_root, "BCR_01_LC.csv"), header = TRUE)
BCR2H <- read.csv(paste0(data_root, "BCR_02_HC.csv"), header = TRUE)
BCR2L <- read.csv(paste0(data_root, "BCR_02_LC.csv"), header = TRUE)
TCR1 <- read.csv(paste0(data_root, "TCR_01.csv"), header = TRUE)
TCR2 <- read.csv(paste0(data_root, "TCR_02.csv"), header = TRUE)

col_title<-colnames(BCR1H)
# "Sample"            "Row_number"        "Total_readcount"   "Shannon_diversity" "Time_point"        "Type" 

##########################################
### Correct error 0 to 1

table(BCR1H$Time_point)
# 0   1   2   3 
# 1 429 186 224 
BCR1H$Time_point <- ifelse(BCR1H$Time_point == 0, 1, BCR1H$Time_point)
table(BCR1H$Time_point)
# 1   2   3 
# 430 186 224 

##########################################
### separate by Type

BCR1H_IGM <- BCR1H[BCR1H$Type == "IGM",] # 420
BCR1H_IGG <- BCR1H[BCR1H$Type == "IGG",] # 420
BCR1L_IGK <- BCR1L[BCR1L$Type == "IGK",] # 420
BCR1L_IGL <- BCR1L[BCR1L$Type == "IGL",] # 420

BCR2H_IGM <- BCR2H[BCR2H$Type == "IGM",] # 110
BCR2H_IGG <- BCR2H[BCR2H$Type == "IGG",] # 110
BCR2H_IGK <- BCR2H[BCR2H$Type == "IGK",] # 61
BCR2H_IGL <- BCR2H[BCR2H$Type == "IGL",] # 101

BCR2L_IGM <- BCR2L[BCR2L$Type == "IGM",] # 101
BCR2L_IGG <- BCR2L[BCR2L$Type == "IGG",] # 110
BCR2L_IGK <- BCR2L[BCR2L$Type == "IGK",] # 110 
BCR2L_IGL <- BCR2L[BCR2L$Type == "IGL",] # 110

TCR1_TRA <- TCR1[TCR1$Type == "TRA",] # 420
TCR1_TRB <- TCR1[TCR1$Type == "TRB",] # 420
TCR2_TRA <- TCR2[TCR2$Type == "TRA",] # 420
TCR2_TRB <- TCR2[TCR2$Type == "TRB",] # 420

all =c("BCR1H_IGM", "BCR1H_IGG", "BCR1L_IGK", "BCR1L_IGL", 
       "BCR2H_IGM", "BCR2H_IGG", "BCR2H_IGK", "BCR2H_IGL", 
       "BCR2L_IGM", "BCR2L_IGG", "BCR2L_IGK", "BCR2L_IGL", 
       "TCR1_TRA", "TCR1_TRB", "TCR2_TRA", "TCR2_TRB")
BCR_1 = c("BCR1H_IGM", "BCR1H_IGG", "BCR1L_IGK", "BCR1L_IGL")
BCR_2 = c("BCR2H_IGM", "BCR2H_IGG", "BCR2H_IGK", "BCR2H_IGL", 
          "BCR2L_IGM", "BCR2L_IGG", "BCR2L_IGK", "BCR2L_IGL")
BCR_H = c("BCR1H_IGM", "BCR1H_IGG",
          "BCR2H_IGM", "BCR2H_IGG", "BCR2H_IGK", "BCR2H_IGL")
BCR_L = c("BCR1L_IGK", "BCR1L_IGL",
          "BCR2L_IGM", "BCR2L_IGG", "BCR2L_IGK", "BCR2L_IGL")
TCR = c("TCR1_TRA", "TCR1_TRB", "TCR2_TRA", "TCR2_TRB")

##########################################
### in Sample ID, remove timepoint and type
### COV-CCO-0411-IGG into COV-CCO-041

# for (df_name in all) {
#   df <- get(df_name)
#   df$Sample <- substring(df$Sample, 1, nchar(df$Sample) - 5)
#   assign(df_name, df)
# }

for (df_name in all) {
  df <- get(df_name)
  df$Sample <- substr(df$Sample, 1, 11)
  assign(df_name, df)
}

#########################################
### Check if T1 - none - T3

plot_list <- list()

for (df_name in TCR) {
  df <- get(df_name)
  
  combi_counts <- df %>%
    group_by(Sample) %>%
    summarise(Combi = toString(unique(Time_point))) %>%
    count(Combi)
  
  plot <- ggplot(combi_counts, aes(x = Combi, y = n)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_text(aes(label = n), vjust = 1.5, size = 3) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8)) +
    labs(title = df_name, x = "Existing Time points", y = "# of samples")
  
  plot_list[[df_name]] <- plot
}

### save just manually (by Zoom)
# png(paste0(save_root, "Timpoint_combi.png"), width = 3000, height = 6800, res = 300)
grid.arrange(grobs = plot_list, ncol = 2)   # check plot, might take long to display on here.
# dev.off()

#########################################

for (df_name in all){
  df <- get(df_name) %>%
    group_by(Sample) %>%
    summarise(Combi = toString(unique(Time_point))) %>%
    left_join(get(df_name), by = "Sample")
  
  assign(df_name, df)
}

#########################################

plot_list02 <- list()

for (df_name02 in BCR_2) {
  df02 <- get(df_name02)

  plot02 <- ggplot(df02, aes(x = factor(Time_point))) +
    geom_bar(fill = "skyblue") +
    geom_text(stat = 'count', aes(label = ..count..), vjust = 1.5, size = 3) +
    scale_x_discrete(labels = c("1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5", "7" = "6", "9" = "7")) +
    theme_bw() +
    labs(title = df_name02, x = "Time point", y = "Cumulative # of samples")

  plot_list02[[df_name02]] <- plot02
}

grid.arrange(grobs = plot_list02, ncol = 2)

##########################################

BCR1_ID_list <- list() 
BCR2_ID_list <- list()

for (df_name in BCR_1) {
  df <- get(df_name)
  BCR1_ID_list <- union(BCR1_ID_list, unique(df$Sample))
  BCR1_ID_list <- unlist(BCR1_ID_list)
}

for (df_name in BCR_2) {
  df <- get(df_name)
  BCR2_ID_list <- union(BCR2_ID_list, unique(df$Sample))
  BCR2_ID_list <- unlist(BCR2_ID_list)
}

print(BCR1_ID_list)
print(BCR2_ID_list)
print(length(BCR1_ID_list))
print(length(BCR2_ID_list))
length(intersect(BCR1_ID_list, BCR2_ID_list))

BCR1_prefix <- substr(BCR1_ID_list, 1, 7)
BCR2_prefix <- substr(BCR2_ID_list, 1, 7)

table(BCR1_prefix)
table(BCR2_prefix)

write.table(BCR1_ID_list, file = paste0(save_root,"BCR1_ID_list.txt"), row.names = FALSE, col.names = FALSE)
write.table(BCR2_ID_list, paste0(save_root,"BCR2_ID_list.txt"), row.names = FALSE, col.names = FALSE)
