library(tidyverse)
library(concaveman)
library(ggplot2)
set.seed(1)

load("../cleaneddata.RData")
cannot <- read.csv("../extraAnalyses/processed_data/cell positions wrt gloms.csv") # loads "cannot"

cellNames <- readRDS("cellNames.Rds")
pred_labels <- readRDS("classifiedLabels.Rds")

macs <- data.frame(cell = cellNames, label = pred_labels)

viz_df <- data.frame(viz)
viz_df$id <- annot$cell_ID
fullSlides <- merge(cannot[, c("id", "position.vs.glom", "inside.glom", "closest.glom")], annot[, c("cell_ID", "tissuename", "fov")], by.x = "id", by.y = "cell_ID")
fullSlides <- merge(fullSlides, viz_df, by = "id")

myeloidExpr <- merge(fullSlides, macs, by.x = "id", by.y = "cell")

dir.create("spatialPlots", showWarnings = FALSE)

for (fov_id in unique(fullSlides$fov)) {
  tissue <- unique(fullSlides$tissuename[fullSlides$fov == fov_id])
  print(fov_id)
  dir.create(paste0("spatialPlots/", tissue), showWarnings = FALSE, recursive = TRUE)
  
  # Compute concave hulls for each position.vs.glom
  fovImage <- fullSlides %>% filter(fov == fov_id)

  coords <- as.matrix(fovImage[, c("sdimx", "sdimy")])

  hull <- concaveman(coords)

  hull_df <- as.data.frame(hull)
  colnames(hull_df) <- c("sdimx", "sdimy")

  # Bordering glomerulous → group by closest.glom
  bg_cells_list <- fovImage %>%
    filter(position.vs.glom == "bordering glomerulous") %>%
    group_by(closest.glom) %>%
    group_split()

  # Inside glomerulous → group by in.glom
  ig_cells_list <- fovImage %>%
    filter(position.vs.glom == "inside glomerulous") %>%
    group_by(inside.glom) %>%
    group_split()

  get_hull <- function(df, concavity = 2){
    coords <- as.matrix(df[, c("sdimx", "sdimy")])
    hull <- concaveman(coords, concavity = concavity, length_threshold = 0)
    hull_df <- as.data.frame(hull)
    colnames(hull_df) <- c("sdimx", "sdimy")
    return(hull_df)
  }

  get_hulls_per_group <- function(df, group_col){
  df %>%
    filter(!is.na(.data[[group_col]])) %>%   # drop rows with NA in the group
    group_by(across(all_of(group_col))) %>%
    group_split() %>%
    lapply(get_hull)                        # get_hull = concaveman wrapper
}

  # Inside glomerulous
  ig_hulls <- get_hulls_per_group(
    fovImage %>% filter(position.vs.glom == "inside glomerulous"),
    "inside.glom"
  )

  if(length(ig_hulls) == 0){
    next  # skip to the next iteration of the loop
  }

  # Bordering glomerulous
  bg_hulls <- get_hulls_per_group(
    fovImage %>% filter(position.vs.glom == "bordering glomerulous"),
    "closest.glom"
  )

  # Add type column to each hull and combine
  ig_hulls_df <- do.call(rbind, lapply(seq_along(ig_hulls), function(i){
    cbind(ig_hulls[[i]], hull_id = i, type = "inside")
  }))

  bg_hulls_df <- do.call(rbind, lapply(seq_along(bg_hulls), function(i){
    cbind(bg_hulls[[i]], hull_id = i, type = "bordering")
  }))


  myeloidExpr_plot <- myeloidExpr %>%
    filter(fov == fov_id) %>%             # keeps only 87 rows
    dplyr::mutate(color_sam = ifelse(label == "SAM", "green", "black"))

  p <- ggplot() +
    geom_polygon(
      data = hull_df,
      aes(x = sdimx, y = sdimy),
      fill = "lightgrey", alpha = 0.6
    ) +
    # Glomeruli hulls
    geom_polygon(data = bg_hulls_df, 
              aes(x = sdimx, y = sdimy, group = hull_id),
              fill = "steelblue1", alpha = 0.6) +
    geom_polygon(data = ig_hulls_df, 
              aes(x = sdimx, y = sdimy, group = hull_id),
              fill = "tomato", alpha = 0.6) + 
    # original points on top
    geom_point(
      data = myeloidExpr_plot %>% filter(label != "SAM"),
      aes(x = sdimx, y = sdimy, shape = label, color = label),
      size = 2,
      alpha = 0.7
    ) +
    geom_point(
      data = myeloidExpr_plot %>% filter(label == "SAM"),
      aes(x = sdimx, y = sdimy, shape = label, color = label),
      size = 2,
      alpha = 0.9
    ) +
    scale_color_manual(values = c(
      "SAM" = "green",
      "ambiguous" = "black",
      "nonSAMs" = "black"
    )) +
    scale_shape_manual(values = c(
      nonSAMs = 15,
      ambiguous = 16,
      SAM = 17
    )) + 
    theme_minimal(base_size = 14) +
    labs(title = paste("Tissue:", tissue),
    x = "width (mm)",
    y = "height (mm)")

  pdf(file = paste0("spatialPlots/", tissue, "/spatialPlots_", fov_id, ".pdf"), width = 5, height = 4)
    print(p)
  dev.off()
}