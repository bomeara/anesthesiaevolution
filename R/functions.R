GetData <- function() {
  raw_data <- read.csv("data/EffectsChart.csv", stringsAsFactors=FALSE, row.names=NULL)
  colnames(raw_data)[grepl("nitrous",colnames(raw_data))] <- "nitrous.oxide"
  raw_data <- raw_data[,-ncol(raw_data)] # deleting empty final column
  colnames(raw_data) <- tolower(colnames(raw_data))
  raw_data <- raw_data[nchar(raw_data[,1])>1,]
  rownames(raw_data) <- trimws(raw_data[,1], which="both")
  raw_data <- raw_data[,-1]
  return(raw_data)
}

AggregateData <- function(raw_data) {
  original_data <- raw_data
  local <- colnames(raw_data)[grepl("caine", colnames(raw_data))]
  general <- colnames(raw_data)[!grepl("caine", colnames(raw_data))]
  raw_data$local <- NA
  raw_data$general <- NA
  raw_data$intravenous <- NA
  raw_data$volatile <- NA
  raw_data$hydrcarbon <- NA
  raw_data$ethers <- NA
  raw_data$any <- 0
  for (species.index in sequence(nrow(raw_data))) {
    local.df <- raw_data[species.index, grepl("caine", colnames(raw_data))]
    if(any(local.df==1, na.rm=TRUE)) {
      raw_data$local[species.index] <- 1
    } else {
      if(any(local.df==0, na.rm=TRUE)) {
        raw_data$local[species.index] <- 0
      }
    }
    if(any(raw_data[species.index,]==1, na.rm=TRUE)) {
      raw_data$any[species.index] <- 1
    }
  }

  for (species.index in sequence(nrow(raw_data))) {
    general.df <- original_data [species.index, !grepl("caine", colnames(original_data))]
    if(any(general.df==1, na.rm=TRUE)) {
      raw_data$general[species.index] <- 1
    } else {
      if(any(general.df==0, na.rm=TRUE)) {
        raw_data$general[species.index] <- 0
      }
    }
  }
  for (species.index in sequence(nrow(raw_data))) {
    iv.df <- original_data[species.index, grepl("ketamine|propofol", colnames(original_data))]
    if(any(iv.df==1, na.rm=TRUE)) {
      raw_data$intravenous[species.index] <- 1
    } else {
      if(any(iv.df==0, na.rm=TRUE)) {
        raw_data$intravenous[species.index] <- 0
      }
    }
  }

  for (species.index in sequence(nrow(raw_data))) {
    volatile.df <- original_data[species.index, grepl("Halothane|chloroform|diethyl.ether|Isoflurane|sevoflurane", colnames(original_data), ignore.case=TRUE)]
    if(any(volatile.df==1, na.rm=TRUE)) {
      raw_data$volatile[species.index] <- 1
    } else {
      if(any(volatile.df==0, na.rm=TRUE)) {
        raw_data$volatile[species.index] <- 0
      }
    }
  }

  for (species.index in sequence(nrow(raw_data))) {
    hydrcarbon.df <- original_data[species.index, grepl("Halothane|chloroform", colnames(original_data), ignore.case=TRUE)]
    if(any(hydrcarbon.df==1, na.rm=TRUE)) {
      raw_data$hydrcarbon[species.index] <- 1
    } else {
      if(any(hydrcarbon.df==0, na.rm=TRUE)) {
        raw_data$hydrcarbon[species.index] <- 0
      }
    }
  }

  for (species.index in sequence(nrow(raw_data))) {
    ethers.df <- original_data[species.index, grepl("diethyl.ether|Isoflurane|sevoflurane", colnames(original_data), ignore.case=TRUE)]
    if(any(ethers.df==1, na.rm=TRUE)) {
      raw_data$ethers[species.index] <- 1
    } else {
      if(any(ethers.df==0, na.rm=TRUE)) {
        raw_data$ethers[species.index] <- 0
      }
    }
  }

  colnames(raw_data) <- tolower(gsub("hydrcarbon", "hydrocarbon", colnames(raw_data)))
  raw_data <- raw_data[,c("any", "local", "general", "intravenous", "volatile", "hydrocarbon", "ethers", "lidocaine", "procaine", "benzocaine", "ketamine", "propofol", "halothane", "chloroform", "diethyl.ether", "isoflurane", "sevoflurane", "nitrous.oxide")]
  return(raw_data)
}


FixCommonNames <- function(aggregate_data) {
  rownames(aggregate_data)[grepl("Rat", rownames(aggregate_data))] <- "Rattus novegicus"
  rownames(aggregate_data)[grepl("Centipede", rownames(aggregate_data))] <- "Scutigera coleoptrata" # 1936 study, so we're assuming "centipede" means house centipede, but we're not certain. It will not affect the phylogeny
  rownames(aggregate_data)[grepl("Katydid", rownames(aggregate_data))] <- "Pterophylla camellifolia" # Also no further detail, so we are picking a katydid that is commonly encountered
  rownames(aggregate_data)[grepl("Rose chafer", rownames(aggregate_data))] <- "Cetonia aurata"
  rownames(aggregate_data)[grepl("Canus familiaris", rownames(aggregate_data))] <- "Canis familiaris"
  rownames(aggregate_data)[grepl("Canis familiaris", rownames(aggregate_data))] <- "Canis lupus familiaris" # since OToL via rphylotastic is returning Canis famliaris
  rownames(aggregate_data)[grepl("Micrococcus lysodeikticus", rownames(aggregate_data))] <- "Micrococcus luteus SK58"

  aggregate_data <- aggregate_data[!grepl("coliform", rownames(aggregate_data), ignore.case=TRUE),] # too general
  aggregate_data <- aggregate_data[!grepl("Staphylococcus albus", rownames(aggregate_data), ignore.case=TRUE),] # synonym
  rownames(aggregate_data)[grepl("Staphylococcus aureus", rownames(aggregate_data))] <- "Staphylococcus aureus subsp. aureus C427"
  rownames(aggregate_data)[grepl("Streptococcus pneumoniae", rownames(aggregate_data))] <- "Streptococcus pneumoniae G54"
  rownames(aggregate_data)[grepl("Escherichia coli", rownames(aggregate_data))] <- "Escherichia coli MS 21-1"
  return(aggregate_data)
}

ResolveNames <- function(aggregate_data) {
  for(i in sequence(nrow(aggregate_data))) {
    try(rownames(aggregate_data)[i] <- rphylotastic::taxa_resolve_names_with_otol(rownames(aggregate_data)[i]))
  }
  return(aggregate_data)
}

GetTrees <- function(aggregate_data) {
  all_chronograms <- datelife::datelife_search(rownames(aggregate_data), use_tnrs=FALSE)
  biggest_chronogram <- datelife::datelife_search(rownames(aggregate_data), use_tnrs=FALSE, summary_format="phylo_biggest")
  otol_subtree <- get_dated_otol_induced_subtree(input= rownames(aggregate_data), use_tnrs=FALSE)

  return(list(all=all_chronograms, biggest=biggest_chronogram, otol=otol_subtree))
}

PlotTreeWithTraits <- function(phy, aggregate_data) {
  pdf(file="plots/bullseye.pdf", width=10, height=10)
  #adephylo::bullseye(phy, aggregate_data, legend=FALSE, traits.inset=1.1)
  bullseye(phy, aggregate_data, legend=FALSE, circ.n=1, show.tip.label=TRUE, cex=0.5,
 traits.cex=1, open.angle=45, rotate=90)
  dev.off()
  print("Made bullseye")
}

PlotPhyloHeatmap <- function(phy, aggregate_data) {
  pdf(file="plots/phyloheatmap.pdf", width=10, height=10)
  traits <- SanitizeData(aggregate_data)
  traits[is.na(traits)] <- 0.5
  phytools::phylo.heatmap(phy, traits, colors=c("blue", "white", "red"), fsize=0.6, legend=FALSE)
  dev.off()
}

PlotIndividualTraits <- function(phy, aggregate_data) {
  phy <- ape::collapse.singles(phy)
  for (trait_index in sequence(ncol(aggregate_data))) {
    focal_trait <- aggregate_data[,trait_index]
    names(focal_trait) <- rownames(aggregate_data)
    focal_trait <- focal_trait[!is.na(focal_trait)]
    pruned <- geiger::treedata(phy, focal_trait, warnings=FALSE, sort=TRUE)
    print(nrow(pruned$data))
    if(nrow(pruned$data)>2 & length(unique(pruned$data[,1]))>1 ) {
      pdf(file=paste0("plots/individual_", gsub(" ", "_", colnames(aggregate_data)[trait_index]), ".pdf"), width=10, height=10)
      simmap.result<-make.simmap(pruned$phy, pruned$data[,1],nsim=100)
      densityMap(simmap.result,states=c("no effect","effect"),plot=TRUE, main=colnames(aggregate_data)[trait_index])
      dev.off()
      print(paste0("Made ", colnames(aggregate_data)[trait_index]))

    }
  }
}

SanitizeNames <- function(x) {
  x <- gsub(" ", "_", x)
  for (i in seq_along(x)) {
    x[i] <- paste(strsplit(x[i], "_")[[1]][c(1,2)], collapse="_")
  }
  x[grepl("Canis", x)] <- "Canis_familiaris"
  return(x)
}

SanitizeData <- function(aggregate_data) {
  rownames(aggregate_data) <- SanitizeNames(rownames(aggregate_data))
  return(aggregate_data)
}

SanitizeTree <- function(phy) {
  phy$tip.label <- SanitizeNames(phy$tip.label)
  return(phy)
}

RunCorHMM <- function(phy, aggregate_data) {
  phy <- ape::collapse.singles(phy)
  trait_results <- list()
  for (trait_index in sequence(ncol(aggregate_data))) {
    print(paste0("trait is ", trait_index, " ", colnames(aggregate_data)[trait_index]))
    focal_trait <- data.frame(Taxon=rownames(aggregate_data), Trait=aggregate_data[,trait_index], stringsAsFactors=FALSE)
    focal_trait <- focal_trait[!is.na(focal_trait$Trait),]
    phy_focal <- phy
    to_cull <- phy$tip.label[!phy$tip.label %in% focal_trait$Taxon]
    if(length(to_cull)>0) {
      phy_focal <- ape::drop.tip(phy, to_cull)
    }
    rates <- c(1:2)
    corHMM_results <- list()
    for (i in seq_along(rates)) {
      local_result <- NA
      try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=rates[i], node.states="marginal", root.p="maddfitz"))
      corHMM_results[[i]] <- local_result
    }

    local_result <- NA
    try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=1, node.states="marginal", root.p=c(1,0)))
    corHMM_results[[length(corHMM_results)+1]] <- local_result

    local_result <- NA
    try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=1, node.states="marginal", root.p=c(0,1)))
    corHMM_results[[length(corHMM_results)+1]] <- local_result

    local_result <- NA
    rate.mat <- corHMM::rate.mat.maker(hrm=TRUE,rate.cat=1)
    rate.mat<-rate.par.drop(rate.mat,c(1))
    try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=1, rate.mat=rate.mat, node.states="marginal", root.p="maddfitz"))
    corHMM_results[[length(corHMM_results)+1]] <- local_result

    local_result <- NA
    rate.mat <- corHMM::rate.mat.maker(hrm=TRUE,rate.cat=1)
    rate.mat<-rate.par.drop(rate.mat,c(2))
    try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=1, rate.mat=rate.mat, node.states="marginal", root.p="maddfitz"))
    corHMM_results[[length(corHMM_results)+1]] <- local_result

    local_result <- NA
    rate.mat <- corHMM::rate.mat.maker(hrm=TRUE,rate.cat=1)
    rate.mat<-rate.par.drop(rate.mat,c(1))
    try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=1, rate.mat=rate.mat, node.states="marginal", root.p=c(1,0)))
    corHMM_results[[length(corHMM_results)+1]] <- local_result

    local_result <- NA
    rate.mat <- corHMM::rate.mat.maker(hrm=TRUE,rate.cat=1)
    rate.mat<-rate.par.drop(rate.mat,c(2))
    try(local_result <- corHMM::corHMM(phy_focal, focal_trait, rate.cat=1, rate.mat=rate.mat, node.states="marginal", root.p=c(0,1)))
    corHMM_results[[length(corHMM_results)+1]] <- local_result



    trait_results[[trait_index]] <- corHMM_results
    names(trait_results[[trait_index]]) <- colnames(aggregate_data)[trait_index]
  }
  return(trait_results)
}

PlotCorHMM <- function(corhmm_results) {
  for(trait in seq_along(corhmm_results)) {
    print(trait)
    trait_name <- names(corhmm_results[[trait]][1])
    best_model <- unname(which.min(unlist(lapply(corhmm_results[[trait]], "[", "AICc"))))[1]
    if(!any(is.na(corhmm_results[[trait]][[1]]))) {
      pdf(file=paste0("plots/corhmm_", gsub(" ", "_", trait_name), ".pdf"), width=10, height=10)
      try(corHMM::plotRECON(corhmm_results[[trait]][[best_model]]$phy, corhmm_results[[trait]][[best_model]]$states))
      dev.off()
    }
  }
}
