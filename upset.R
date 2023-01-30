set1 <- read.csv(file = "Plants.txt", header = FALSE, sep = "\n")
set2 <- read.csv(file = "Soils.txt", header = FALSE, sep = "\n")
set3 <- read.csv(file = "Extreme.txt", header = FALSE, sep = "\n")
set4 <- read.csv(file = "Marine.txt", header = FALSE, sep = "\n")
set5 <- read.csv(file = "Air.txt", header = FALSE, sep = "\n")

set_1 <- as.vector(set1$V1)
set_2 <- as.vector(set2$V1)
set_3 <- as.vector(set3$V1)
set_4 <- as.vector(set4$V1)
set_5 <- as.vector(set5$V1)

Habitat_sets = list(Plants = set_1, Soils = set_2, Extreme = set_3, Marine = set_4, Air = set_5)
upset(fromList(Habitat_sets),
sets = c("Plants", "Soils", "Extreme", "Marine", "Air"),
number.angles = 20, point.size = 2.5, line.size = 1.5,
mainbar.y.label = "Habitat intersection", sets.x.label = "Habitat set size",
text.scale = c(1.5, 1.5, 1.25, 1.25, 1.5, 1.5), mb.ratio = c(0.65, 0.35),
group.by = "freq", keep.order = TRUE)
