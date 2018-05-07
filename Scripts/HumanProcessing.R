human <- read.table('~/Documents/Homo_sapiens_physical')

human <- data.frame(human$INTERACTOR_A, human$INTERACTOR_B)

human_uniq <- unique(human)
library(plyr)
col1 <- count(human_uniq, 'human.INTERACTOR_A')
col2 <- count(human_uniq, 'human.INTERACTOR_B')

colnames(col1) <- c("protein", "freq")
colnames(col2) <- c("protein", "freq")
colBoth <- merge(col1, col2, by=c("protein"))
colBoth$count <- (colBoth$freq.x + colBoth$freq.y)

#Find Top 25 Physical
top25 <- data.frame(colBoth$protein, colBoth$count)
top25 <- top25[order(-top25$colBoth.count),]
list25 <- data.frame(top25$colBoth.protein[1:25], top25$colBoth.count[1:25])
write.table(list25, "~/Desktop/top_human.txt", sep="\t")

#Log-log plotting
human_interactions <- data.frame(colBoth$protein, colBoth$count)
total_count <- data.frame(human_interactions$colBoth.count)

human_dist <- count(total_count, 'human_interactions.colBoth.count')

plot(human_dist$human_interactions.colBoth.count, (human_dist$freq), log="xy", xlim=c(1,max(human_dist$human_interactions.colBoth.count)), ylim=c(1,max(human_dist$freq)), main="human Node Degree Distribution (Log Scale)", xlab="Number of Connections", ylab="Number of Nodes", pch=19)

#Linear Fitting
plot(log(human_dist$freq) ~ log(human_dist$human_interactions.colBoth.count), main="human Node Degree Distribution (Log Scale)", xlab="Number of Connections", ylab="Number of Nodes")
fit <- lm(log(human_dist$freq) ~ log(human_dist$human_interactions.colBoth.count))
coef(fit)
abline(coef(fit)[1], coef(fit)[2])