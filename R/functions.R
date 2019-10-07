GetData <- function() {
  raw_data <- read.csv("data/EffectsChart.csv", stringsAsFactors=FALSE, row.names=1)
  colnames(raw_data)[grepl("nitrous",colnames(raw_data))] <- "nitrous.oxide"
  raw_data <- raw_data[,-ncol(raw_data)] # deleting empty final column
  colnames(raw_data) <- tolower(colnames(raw_data))
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
  for (species.index in sequence(nrow(raw_data))) {
    local.df <- raw_data[species.index, grepl("caine", colnames(raw_data))]
    if(any(local.df==1, na.rm=TRUE)) {
      raw_data$local[species.index] <- 1
    } else {
      if(any(local.df==0, na.rm=TRUE)) {
        raw_data$local[species.index] <- 0
      }
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

  return(raw_data)
}


FixCommonNames <- function(aggregate_data) {
  rownames(aggregate_data)[grepl("Rat", rownames(aggregate_data))] <- "Rattus novegicus"
  rownames(aggregate_data)[grepl("Centipede", rownames(aggregate_data))] <- "Scutigera coleoptrata" # 1936 study, so we're assuming "centipede" means house centipede, but we're not certain. It will not affect the phylogeny
  rownames(aggregate_data)[grepl("Katydid", rownames(aggregate_data))] <- "Pterophylla camellifolia" # Also no further detail, so we are picking a katydid that is commonly encountered
  rownames(aggregate_data)[grepl("Rose chafer", rownames(aggregate_data))] <- "Cetonia aurata"
  rownames(aggregate_data)[grepl("Canus familiaris", rownames(aggregate_data))] <- "Canis familiaris"
  aggregate_data <- aggregate_data[!grepl("coliform", rownames(aggregate_data), ignore.case=TRUE),] # too general
  aggregate_data <- aggregate_data[!grepl("Staphylococcus albus", rownames(aggregate_data), ignore.case=TRUE),] # synonym
  return(aggregate_data)
}

ResolveNames <- function(aggregate_data) {
  for(i in sequence(nrow(aggregate_data))) {
    rownames(aggregate_data)[i] <- rphylotastic::taxa_resolve_names_with_otol(rownames(aggregate_data)[i])
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
  adephylo::bullseye(phy, aggregate_data, legend=FALSE, traits.inset=1.1)
  dev.off()
}
