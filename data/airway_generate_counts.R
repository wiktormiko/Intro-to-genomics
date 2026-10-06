library(airway)
data("airway")

# save raw counts
cnts <- assay(airway)
write.csv(cnts, "data/airway_raw_counts.csv")

# save metadata
mtd <- colData(airway)
write.csv(mtd, "data/airway_metadata.csv")
