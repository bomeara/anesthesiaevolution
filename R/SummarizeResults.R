library(adephylo)
library(RColorBrewer)

setwd("/Users/bomeara/Documents/MyDocuments/Active/WgPlay/March2015")
load("combinedData.RSave")
rm(list=ls()[!grepl("primate", ls())])
source("AnalysisFunctions.R")





#need to get dependent columns as logical


all.cols <- sequence(dim(primate.data)[2])
factor.cols <- sort(c(1, 26, 28, 38, 39:46, 122, 123, 127, 128, 134, 115, 117:120, 135:138, 147:150, 152:156, 33:36))
for (col.number in factor.cols) {
	primate.data[,col.number] <- as.factor(primate.data[, col.number])
}
non.factor.cols <- all.cols[-factor.cols]
for (col.number in non.factor.cols) {
	primate.data[,col.number] <- as.numeric(primate.data[, col.number])
}

primate.data.pruned <- primate.data



#primate.data.pruned <- primate.data[which(primate.data$overall.play.present==1),] #b/c we believe all primates play, so if this isn't the case, taxon might not be studied enough

primate.data.pruned$percent.leaves.in.diet.reconciled <- apply(cbind(as.numeric(primate.data.pruned$percent.leaves.in.diet), 100*as.numeric(primate.data.pruned$PropLeaves)), 1, mean, na.rm=TRUE)

primate.data.pruned$GroupSize.reconciled <- apply(cbind(as.numeric(primate.data.pruned$GS), as.numeric(primate.data.pruned$GroupSize)), 1, mean, na.rm=TRUE)


cols.to.delete <- c("Your.name", "Your.email", "species", "Genus", "species.1", "PropLeaves", "percent.leaves.in.diet", "GS", "GroupSize")
primate.data.pruned <- primate.data.pruned[, -which(colnames(primate.data.pruned) %in% cols.to.delete)]

all.independent<-c("BMR_mean", "percent.leaves.in.diet.reconciled", "percent.insects.and.animals.in.diet", "GroomingTimePercent", "log.weaning.months", "log.age.first.reproduction.months", "Habitat", "GroupSize.reconciled", "Substrate.y", "diet")
all.independent.nice.names <- c("BMR", "Leaf consumption", "Insect/animal consumption", "Grooming time", "Weaning age", "Reproduction age", "Habitat", "Group size", "Substrate", "Diet")


all.dependent <- c("overall.sex.play", "overall.social.play")
all.dependent.nice.names <- c("Sexual play", "Social play")

primate.data.pruned <- primate.data.pruned[, c(all.dependent, all.independent)]
primate.data.pruned$log.BMR_mean <- log(primate.data.pruned$BMR_mean)
all.independent[which(all.independent=="BMR_mean")]<-"log.BMR_mean"
primate.data.pruned$log.GroupSize.reconciled <- log(primate.data.pruned$GroupSize.reconciled)
all.independent[which(all.independent=="GroupSize.reconciled")]<-"log.GroupSize.reconciled"
all.independent.small <- all.independent



primate.data.to.plot <- primate.data.pruned[, c(all.dependent, all.independent)]
colnames(primate.data.to.plot) <- c(all.dependent.nice.names, all.independent.nice.names)
cols.to.factor <- c("Habitat", "Substrate", "Diet")
for (i in sequence(dim(primate.data.to.plot)[2])) {
	if(colnames(primate.data.to.plot)[i] %in% 	cols.to.factor ) {
			primate.data.to.plot[,i]<-as.numeric(as.factor(primate.data.to.plot[,i]))
	} else {
		primate.data.to.plot[,i]<-as.numeric(primate.data.to.plot[,i])

	}
}

pdf(file="RawDataPlot.pdf", width=7, height=7)
 #table.phylo4d(phylo4d(primate.phy, primate.data.to.plot), box=FALSE, cex.label=0.2, legend=FALSE, cex.symbol=0.2, grid=FALSE, col=brewer.pal(5,"RdYlGn"))
 par(mar=rep(10,4))
 bullseye(primate.phy, primate.data.to.plot, legend=FALSE, circ.n=1, show.tip.label=TRUE, cex=0.2, 
traits.cex=0.5, open.angle=45, rotate=90)
 dev.off()
 system("open RawDataPlot.pdf")
 cat(rev(c(all.dependent.nice.names, all.independent.nice.names)), sep="\n")

#USE phyloglm in package phylolm

#foreach (i=1:length(all.dependent)) %dopar% DredgeAnalyses(all.independent, all.dependent[i], primate.data.pruned, primate.phy, max.fields=3, intermediate.file=paste("Intermediate_", all.dependent[i], ".Rdata", sep=""), optimize.kappa=FALSE)


#all.independent.small<-c("SocialSystem", "Activity", "Substrate.y", "Habitat", "PropLeaves", "log.species.mass.grams", "percent.leaves.in.diet", "log.weaning.months", "log.age.first.reproduction.months", "percent.insects.and.animals.in.diet",  "GroupSize", "GroomingTimePercent", "BMR_mean")
all.independent.small <- all.independent




results <- system("ls -1 *Rdata | grep Intermediate_Small", intern=TRUE)
result.list <-NULL
results.short.names <- gsub("Intermediate_SmallRun_overall.", "", gsub(".Rdata", "", results))

for (result.index in sequence(length(results))) {
	print(result.index)
	load(results[result.index])
	source("AnalysisFunctions.R")
	all.results.filter <- FilterDredge(all.results)
	print(length(all.results.filter))
	if(length(all.results.filter)>0) {
		print(paste("saving to ", gsub("Rdata", "xlsx", results[result.index])))
		all.results.aggregate <- AggregateResults(all.results.filter)
		WriteAggregateResults(all.results.aggregate, file.name=gsub("Rdata", "xlsx", results[result.index]))
		summary.aggregate <- CountParameterInfluenceAllModelSets(all.results.aggregate)
#		
		if(result.index==1) {
			result.list <- list(summary.aggregate)
		} 
		if(result.index==2) {
			result.list <- list(result.list[[1]], summary.aggregate)
		} 
		if(result.index==3) {
			result.list <- list(result.list[[1]], result.list[[2]], summary.aggregate)
		} 
		write.xlsx(summary.aggregate, file=gsub("Rdata", "summary.xlsx", results[result.index]))
		#try(system(paste("open ", gsub("Rdata", "summary.xlsx", results[result.index]))))
	} else {
		print("all failed")
	}
	rm(all.results)
	rm(all.results.filter)
	rm(all.results.aggregate)
}

result.list.pruned <- result.list
for (i in sequence(length(result.list.pruned))) {
	result.list.pruned[[i]] <- result.list.pruned[[i]][result.list.pruned[[i]]$trait %in% all.independent,]
	colnames(result.list.pruned[[i]]) <- paste(results.short.names[[i]], colnames(result.list.pruned[[i]]), sep="_")
}

results.all <- result.list.pruned[[1]]
for (i in 2:length(result.list.pruned)) {
		results.all <- cbind(results.all, result.list.pruned[[i]])
}

results.all$Trait <- all.independent.nice.names[order(all.independent, decreasing=TRUE)]
write.xlsx(results.all, file="LogisticRegressionResults.xlsx")
