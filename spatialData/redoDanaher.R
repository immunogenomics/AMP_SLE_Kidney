library(spatstat.geom)
library(spatstat)

rm(list = ls())
gloms = readRDS("extraAnalyses/processed_data/gloms.rds")
db = readRDS("extraAnalyses/processed_data/db.rds")
cols <- readRDS("extraAnalyses/processed_data/cols.RDS")
load("cleaneddata.RData")
rownames(customlocs) = annot$cell_ID
rownames(viz) = annot$cell_ID

#### new attempt 5-12: annotate cells by position vis-a-vis glom --------------------------------------------

### for each glom, add neighboring Parietal.epithelium cells to the db clust:
## for each parietal epi cell, record nearest neighbor from !is.na(db) cells, and map to it
# point pattern of pe cells and db!=NA cells:
#use = is.element(names(clust), names(db)[!is.na(db)]) | (clust == "Parietal.epithelium")
use = rep(TRUE, nrow(customlocs)); names(use) = rownames(customlocs)
pp <- spatstat.geom::ppp(customlocs[use, 1], customlocs[use, 2],
                         xrange = range(customlocs[use, 1]), yrange = range(customlocs[use, 2]))
marks(pp) <- db[match(names(which(use)), names(db))]
#marks(pp)[is.na(marks(pp))] = "na"
marks(pp) = as.factor(marks(pp))
# count neighbors of each db cluster:
mt05 <- marktable(X = pp, R = 0.05, N = NULL, exclude=TRUE, collapse=FALSE)
mt01 <- marktable(X = pp, R = 0.01, N = NULL, exclude=TRUE, collapse=FALSE)
rownames(mt01) <- rownames(mt05) <- names(which(use))

## assign cells to their most neighboring glom:
# which glom wins the vote for each cell?
closestglom <- colnames(mt05)[apply(mt05,1,which.max)]
closestglom[rowSums(mt05) == 0] = NA
names(closestglom) = rownames(mt05)


## assign parietal epi cells to gloms with a more stringent radius:
closestglom.pe <- colnames(mt01)[apply(mt01,1,which.max)]
closestglom.pe[rowSums(mt01) == 0] = NA

is.glom.or.pe = is.element(names(clust), names(db)[!is.na(db)]) | (clust == "Parietal.epithelium")
closestglom.pe[!is.glom.or.pe] = NA
names(closestglom.pe) = rownames(mt01)

# made db object that includes pe cells along with the previous glom cells:
db.pe = db
db.pe[names(which(!is.glom.or.pe))] = closestglom.pe[names(which(!is.glom.or.pe))]

### start new cannot:
cannot = data.frame(id = annot$cell_ID)
rownames(cannot) = cannot$id

# variable for members of a glom (glomerular cells only):
cannot$in.glom = NA
cannot[names(db.pe), "in.glom"] = db.pe
# remove glom ids not present in "gloms" object:
removedgloms = setdiff(db.pe, c(rownames(gloms), NA))
print("removing these gloms from db.pe:")
print(removedgloms)
cannot$in.glom[is.element(cannot$in.glom, removedgloms)] = NA
table(clust[!is.na(cannot$in.glom)])

### for each glom, define a polygon:
gpolys = list()
for (gid in setdiff(unique(db.pe), NA)) {
  tempcells = names(db.pe)[(db.pe == gid) & !is.na(db.pe)]
  tempinds = chull(customlocs[tempcells, ])
  tempinds = c(tempinds, tempinds[1])
  gpolys[[gid]] = customlocs[tempcells[tempinds], ]
}

# test the polys:
png(paste0("extraAnalyses/results/glom polys test xy.png"), width = 17, height = 17, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, pch = 16, 
     cex = 0.01 + 0.05 * is.element(clust, c("Podocyte", "Mesangial.cell","Glomerular.endothelium","Parietal.epithelium")),
     col = cols[clust], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpolys)) {
  polygon(gpolys[[i]], border = "black", lwd = 0.3)
}
dev.off()

### record all cells within each polygon:
cannot$inside.glom = NA
for (i in 1:length(gpolys)) {
  if (i%%20==0){print(i)}
  insidethispoly = sp::point.in.polygon(point.x = customlocs[, 1], 
                                        point.y = customlocs[, 2], 
                                        pol.x = gpolys[[i]][, 1], 
                                        pol.y = gpolys[[i]][, 2], mode.checked=FALSE)
  cannot$inside.glom[insidethispoly > 0] = names(gpolys)[i]
}
#table(cannot$inside.glom)
#table(is.na(cannot$in.glom))
#table(is.na(cannot$inside.glom)) # great

### for each cell: what glom is it near?
cannot$closest.glom = closestglom

cannot$position.vs.glom = "tubulointerstitium"  
cannot$position.vs.glom[!is.na(cannot$closest.glom)] = "bordering glomerulous"
cannot$position.vs.glom[!is.na(cannot$inside.glom)] = "inside glomerulous"

png(paste0("extraAnalyses/results/glom position.vs.glom test xy.png"), width = 17, height = 17, res = 600, units = "in")
par(mar = c(0,0,0,0))
plot(customlocs, 
     pch = 16, #c(1, 16, 2)[(cannot$position.vs.glom == "inside glomerulous") + 2*(cannot$position.vs.glom == "bordering glomerulous") + 3*(cannot$position.vs.glom == "tubulointerstitium")], 
     cex = 0.02,
     col = c("red", "blue", "grey80")[
       (cannot$position.vs.glom == "inside glomerulous") + 2*(cannot$position.vs.glom == "bordering glomerulous") + 3*(cannot$position.vs.glom == "tubulointerstitium")
     ],
     #col = cols[clust], 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
for (i in 1:length(gpolys)) {
  polygon(gpolys[[i]], border = "black", lwd = 0.3)
}
dev.off()

# save
save(cannot, file = "extraAnalyses/processed_data/cannot with distances from gloms.RData")
write.csv(cannot, file = "extraAnalyses/processed_data/cell positions wrt gloms.csv", row.names = F)