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

plot_df <- clinData %>%
  count(`Final_ISN[III]`, `Final_ISN[III][V]`)

# Step 2: Plot with size = count
p <- ggplot(plot_df, aes(x = `Final_ISN[III]`, y = `Final_ISN[III][V]`, size = n)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Count") +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(
    title = "Scatter of One-Hot Encoded Variables",
    x = "Final_ISN[III]",
    y = "Final_ISN[III][V]"
  ) +
  theme_minimal() +
  coord_fixed()  # makes the plot square

pdf("iii_and_iii_v.pdf", width = 5, height = 5)
    print(p)
dev.off()