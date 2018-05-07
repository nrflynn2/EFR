plant <- read.table('~/Documents/arabidopsis_thaliana_physical')

plant <- data.frame(plant$INTERACTOR_A, plant$INTERACTOR_B)

plant_uniq <- unique(plant)
library(plyr)
col1 <- count(plant_uniq, 'plant.INTERACTOR_A')
col2 <- count(plant_uniq, 'plant.INTERACTOR_B')

colnames(col1) <- c("protein", "freq")
colnames(col2) <- c("protein", "freq")
colBoth <- merge(col1, col2, by=c("protein"))
colBoth$count <- (colBoth$freq.x + colBoth$freq.y)

#Find Top 25 Physical
top25 <- data.frame(colBoth$protein, colBoth$count)
top25 <- top25[order(-top25$colBoth.count),]
list25 <- data.frame(top25$colBoth.protein[1:25], top25$colBoth.count[1:25])
write.table(list25, "~/Desktop/top_plant.txt", sep="\t")

#Log-log plotting
plant_interactions <- data.frame(colBoth$protein, colBoth$count)
total_count <- data.frame(plant_interactions$colBoth.count)

plant_dist <- count(total_count, 'plant_interactions.colBoth.count')

plot(plant_dist$plant_interactions.colBoth.count, (plant_dist$freq), log="xy", xlim=c(1,max(plant_dist$plant_interactions.colBoth.count)), ylim=c(1,max(plant_dist$freq)), main="plant Node Degree Distribution (Log Scale)", xlab="Number of Connections", ylab="Number of Nodes", pch=19)

#Linear Fitting
plot(log(plant_dist$freq) ~ log(plant_dist$plant_interactions.colBoth.count), main="plant Node Degree Distribution (Log Scale)", xlab="Number of Connections", ylab="Number of Nodes")
fit <- lm(log(plant_dist$freq) ~ log(plant_dist$plant_interactions.colBoth.count))
coef(fit)
abline(coef(fit)[1], coef(fit)[2])