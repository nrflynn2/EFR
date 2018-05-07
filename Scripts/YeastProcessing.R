yeast <- read.table('~/Downloads/Saccharomyces_cerevisiae_physical')

yeast <- data.frame(yeast$INTERACTOR_A, yeast$INTERACTOR_B)

yeast_uniq <- unique(yeast)
library(plyr)
col1 <- count(yeast_uniq, 'yeast.INTERACTOR_A')
col2 <- count(yeast_uniq, 'yeast.INTERACTOR_B')

colnames(col1) <- c("protein", "freq")
colnames(col2) <- c("protein", "freq")
colBoth <- merge(col1, col2, by=c("protein"))
colBoth$count <- (colBoth$freq.x + colBoth$freq.y)

#Find Top 25 Physical
top25 <- data.frame(colBoth$protein, colBoth$count)
top25 <- top25[order(-top25$colBoth.count),]
list25 <- data.frame(top25$colBoth.protein[1:25], top25$colBoth.count[1:25])
write.table(list25, "~/Desktop/top_fungi.txt", sep="\t")

#Log-log plotting
yeast_interactions <- data.frame(colBoth$protein, colBoth$count)
total_count <- data.frame(yeast_interactions$colBoth.count)

yeast_dist <- count(total_count, 'yeast_interactions.colBoth.count')

plot(yeast_dist$yeast_interactions.colBoth.count, (yeast_dist$freq), log="xy", xlim=c(1,max(yeast_dist$yeast_interactions.colBoth.count)), ylim=c(1,max(yeast_dist$freq)), main="Yeast Node Degree Distribution (Log Scale)", xlab="Number of Connections", ylab="Number of Nodes", pch=19)

#Linear Fitting
plot(log(yeast_dist$freq) ~ log(yeast_dist$yeast_interactions.colBoth.count), main="Yeast Node Degree Distribution (Log Scale)", xlab="Number of Connections", ylab="Number of Nodes")
fit <- lm(log(yeast_dist$freq) ~ log(yeast_dist$yeast_interactions.colBoth.count))
coef(fit)
abline(coef(fit)[1], coef(fit)[2])