library(tidyverse)
library(ggplot2)
library(extrafont)
library(pals)
library(janitor)

font_import(prompt = FALSE)
loadfonts()

# Load data
clinData <- readRDS("/data/srlab/ssg34/SLE_kidney_v2/sugiarto_scripts/finalObjects/clinicalData.rds")
all_cell_meta <- read.csv('/data/srlab/ssg34/share/for_Roopa/QC_Annotated_Full_Kidney_Meta_sc_sn_03292024.csv')

# Take only patients whose single cell data is used
clinData <- clinData[clinData$sample %in% all_cell_meta$sample, ]
first_biop_pred <- readRDS("/data/srlab2/qxiao/AMP-SLE/data/clinical/df_pred_biop.rds") 
clinData <- merge(clinData, first_biop_pred, by.x = "AMP.Subject_ID", by.y = "Subject_ID")

# One-hot encode ISN
ISN <- as.data.frame(model.matrix(~Final_ISN-1, data=clinData))
ISN$sample <- clinData$sample[!(is.na(clinData$Final_ISN))]

clinData <- merge(clinData, ISN, by = "sample", , all.x = TRUE)

# One hot encode sex and race
clinData$Sex <- ifelse(clinData$Sex == "Male", 1, 0)
clinData$Hispanic <- ifelse(clinData$Ethnicity == "Hispanic or Latino", 1, 0)
clinData$Asian <- ifelse(grepl("\\[A", clinData$Race), 1, 0)
clinData$Black <- ifelse(grepl("\\[B", clinData$Race), 1, 0)
clinData$White <- ifelse(grepl("\\[W", clinData$Race), 1, 0)

# Discretize responder status
clinData <- clinData %>%
  mutate(Treatment_Response = case_when(
    Responder.Status == "NR" ~ 0,
    Responder.Status == "PR" ~ 1,
    Responder.Status == "CR" ~ 2,
    TRUE ~ NA_real_  # Handles unexpected values with NA
  ))

clinData <- clinData[, c(
    "sample",
    "Age",
    "Sex",
    "Site",
    "Final_Activity",
    "Final_Chronicity",
    "Pred_use",
    "First_biop",
    "Final_ISN[III]",
    "Final_ISN[III][V]",
    "Final_ISN[IV]",
    "Final_ISN[IV][V]",
    "Final_ISN[V]",
    "Hispanic",
    "Asian",
    "Black",
    "White",
    "Treatment_Response")
]

# Correlate two variables only for their non-NA values
correlate <- function(response, predictor) {
    temp <- clinData_clean[!(is.na(clinData_clean[[response]])) & !(is.na(clinData_clean[[predictor]])), ]
    return(cor(temp[[response]], temp[[predictor]])^2)
}

# Initialize dataframe
corrs <- data.frame(Var1 = character(), Var2 = character(), R2 = numeric())
addRow <- function(predictor, response){
    df <- data.frame(Var1 = c(predictor), Var2 = c(response), R2 = c(correlate(predictor, response)))
    return(df)
}

univariate <- c(
    "age", 
    "sex",  
    "final_activity", 
    "final_chronicity", 
    "pred_use", 
    "first_biop", 
    "final_isn_iii", 
    "final_isn_iii_v", 
    "final_isn_iv", 
    "final_isn_iv_v", 
    "final_isn_v", 
    "hispanic", 
    "asian", 
    "black", 
    "white", 
    "treatment_response"
)

clinData_clean <- clinData %>%
  clean_names()

siteCounts <- data.matrix(table(clinData_clean$site))

# Take only site with > 5 samples
clinData_clean <- clinData_clean[clinData_clean$site %in% row.names(siteCounts)[siteCounts[,1] > 5], ]

for (i in 1:length(univariate)) {
    for (j in i:length(univariate)) {
        corrs <- rbind(corrs, addRow(univariate[i], univariate[j]))
    }
}

# R for a categorical variable
categoricalR <- function(response) {
    formula <- as.formula(paste(response, "~ site"))
    model <- lm(formula, data = clinData_clean[!(is.na(clinData_clean[[response]])), ])
    return(summary(model)$r.squared)
}

for (var in univariate) {
    df <- data.frame(Var1 = c("site"), Var2 = c(var), R2 = c(categoricalR(var)))
    corrs <- rbind(corrs, df)
}
corrs <- rbind(corrs, data.frame(Var1 = c("site"), Var2 = c("site"), R2 = c(1)))

corrs_sym <- corrs %>%
  rename(var1_new = Var2, var2_new = Var1) %>%
  rename(Var1 = var1_new, Var2 = var2_new)

corrs <- rbind(corrs, corrs_sym)
corrs <- distinct(corrs)

order <- c(
    "age", 
    "sex", 
    "pred_use", 
    "first_biop",        
    "site", 
    "hispanic",          
    "asian", 
    "black", 
    "white",             
    "treatment_response",
    "final_isn_iii", 
    "final_isn_iii_v", 
    "final_isn_iv",      
    "final_isn_iv_v", 
    "final_isn_v",  
    "final_chronicity", 
    "final_activity"    
)

corrs$Var1 <- factor(corrs$Var1, levels = order, ordered = TRUE)
corrs$Var2 <- factor(corrs$Var2, levels = order, ordered = TRUE)

corrs_triangle <- corrs %>% filter(Var2 >= Var1)  # Keep only upper triangle

corrs_triangle$label <- round(corrs_triangle$R2, 2)

plot <- ggplot(data = corrs_triangle, aes(x = Var2, y = Var1, fill = R2)) + 
    geom_tile(height = 0.75, width = 0.75) + 
    geom_text(aes(label = label), size = 2) + 
    theme_classic() +
    labs(fill = "Rsquared") + 
    scale_fill_gradient(low = "white", high = "red", 
                      limits = c(0, 0.2), 
                      oob = scales::squish) + 
    theme(legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90),
          text=element_text(family="Arial Light"))

pdf("correlations.pdf", width = 7, height = 7)
    print(plot)
dev.off()
