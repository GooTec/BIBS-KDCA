library(tidyverse)
library(glmnet)
library(ggplot2)
library(gridExtra)

# Working directory 설정
dir <- "C:/Users/"
setwd(dir)

# 데이터 및 결과 파일 경로 설정
data_root <- "KDCA_forR/Data/"
save_root <- "KDCA_forR/Result/"
data_root2 <- "정우빈/KDCA/Data/Clinical"

init <- read.csv(paste0(data_root, "COVID_initial_test.csv"), header = TRUE)
multiple <- read.csv(paste0(data_root, "COVID_multi_timepoints_test.csv"), header = TRUE)

cl1 <- read.csv(paste0(data_root, "clinical_data_1.csv"), header = TRUE)
cl2 <- read.csv(paste0(data_root, "clinical_data_2.csv"), header = TRUE)

######################################
### 5자리 연속 숫자 확인법

any(apply(init, 2, function(x) any(x == 66666, na.rm = TRUE)))
any(apply(init, 2, function(x) any(x == 77777, na.rm = TRUE)))
any(apply(init, 2, function(x) any(x == 99999, na.rm = TRUE)))

numbers <- c(66666, 77777, 99999)
result_df <- data.frame(Number = numbers)

for (col_name in names(init)) {
  freq <- table(init[[col_name]], exclude = NULL)
  result_df[[col_name]] <- sapply(numbers, function(x) freq[as.character(x)])
  result_final <- result_df[, colSums(is.na(result_df)) < nrow(result_df)]
}
print(result_df)

df66666 <- result_df[1, !is.na(result_df[1,])]
df77777 <- result_df[2, !is.na(result_df[2, ])]
df99999 <- result_df[3, !is.na(result_df[3, ])]

print(df66666)
print(df77777)

#####################################

start <- which(colnames(df66666) == 'FEV')
end <- which(colnames(df66666) == 'ETC')
target_col <- colnames(df66666)[start:end]
target_col2 <- paste("M1_", target_col, sep = "")
target_col <- paste("CMD_", target_col, sep = "")

target25_66666 <- init[init$FEV == 66666, ]
id25_66666 <- target25_66666$ID

target152_77777 <- init[init$FEV == 77777, ]
id152_77777 <- target152_77777$ID

intersect(id25_66666, id152_77777)    ### might be "0"

###

find25 <- cl1[cl1$ID %in% id25_66666, ]
find25 <- find25[, c("ID","PNTTM","CMD_ETCEXP",target_col)]
cat(paste(target_col, collapse = ", "))
print(length(target_col))

table_results <- lapply(find25, table)
print(table_results)

###

find152 <- cl2[cl2$ID %in% id152_77777, ]
find152 <- find152[, c("ID","PNTTM","M1_ETCEXP",target_col2)]
cat(paste(target_col2, collapse = ", "))
print(length(target_col2))

table_results <- lapply(find152, table)
print(table_results)

subset <- (find152[grep("CNC", find152$ID), ])

########################################

colname66666 <- colnames(df66666)
colname66666 <- colname66666[colname66666 %in% colnames(init)]
colname66666 <- paste("CMD_", colname66666, sep = "")
colname66666 <-intersect(colname66666, colnames(cl1))

target = data.frame()
for (col in colname66666) {
  trimmed <- cl1[, c("ID", "PNTTM", "CMD_ETCEXP", col)]
  trimmed <- trimmed[trimmed$PNTTM != 1, ]
  trimmed$Colinfo = col
  
    apply(trimmed, 1, function(row) {
    if (!is.na(row["CMD_ETCEXP"]) && (row["CMD_ETCEXP"] != 66666) || !is.na(row[col]) && any(row[col] != 66666)) {
      target <<- rbind(target, row)
    }
  })
}
target
table(target$NA_character_.)
yesclinic <- target[grep(1, target$NA_character_.), ]

#################################################
### have extra comorbidities that not includeded as column (1= yes, 2= No, 66666= Not account)

underexp <- target$X.COV.CCO.013.[target$X.CMD_CKIDNEYP.=="CMD_UNDEREXP"]

isthere <- data.frame(ID = numeric(0), PNTTM = character(0), CMD_UNDER = numeric(0))
for (each in underexp) {
  rows <- cl1[cl1$ID == each, c("ID", "PNTTM", "CMD_UNDER")]
  isthere <- rbind(isthere, rows)
}
table(isthere$CMD_UNDER)
length(unique(isthere$ID))

trimmed66666 <- isthere[duplicated(isthere$ID) & (isthere$CMD_UNDER == 1), ]
trimmed66666 <- na.omit(trimmed66666)
length(unique(trimmed66666$ID))

length(yesclinic[yesclinic$X.CMD_CKIDNEYP.== 'CMD_AVPU__1',])

##############################################
##############################################

colname77777 <- colnames(df77777)
colname77777 <- colname77777[colname77777 %in% colnames(init)]
colnameA <- paste("M1_", colname77777, sep = "")
colnameB <- paste("M2_", colname77777, sep = "")
colnameC <- paste("M3_", colname77777, sep = "")
colnameD <- paste("M4_", colname77777, sep = "")
colnameA <-intersect(colnameA, colnames(cl2))
colnameB <-intersect(colnameB, colnames(cl2))
colnameC <-intersect(colnameC, colnames(cl2))
colnameD <-intersect(colnameD, colnames(cl2))
totalcol = c(colnameA,colnameB,colnameC,colnameD)

target = data.frame()
for (col in totalcol) {
  trimmed02 <- cl2[, c("ID", "PNTTM", "M1_ETCEXP", col)]
  trimmed02 <- trimmed02[trimmed02$PNTTM != 1, ]
  trimmed02$Colinfo = col
  
  apply(trimmed, 1, function(row) {
    if (!is.na(row["CMD_ETCEXP"]) && (row["CMD_ETCEXP"] != 77777) || !is.na(row[col]) && any(row[col] != 77777)) {
      target <<- rbind(target, row)
    }
  })
}

target
colnames(target)
table(target$X.99999..1)
yesclinic <- target[grep(2, target$X.99999..1), ]

unique(target$X.M2_BDRUG.)

#################################################################

unique(cl1$CMD_UNDEREXP)
sum(grepl("lipidemia", cl1$CMD_UNDEREXP)) # 30
sum(grepl("dyslipidemia", cl1$CMD_UNDEREXP)) 
sum(grepl("Hyperlipidemia", cl1$CMD_UNDEREXP)) # 13
sum(grepl("hyperlipidemia", cl1$CMD_UNDEREXP)) # 3
