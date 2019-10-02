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
    general.df <- original_data[species.index, !grepl("caine", colnames(original_data))]
    if(any(general.df==1, na.rm=TRUE)) {
      raw_data$general[species.index] <- 1
    } else {
      if(any(general.df==0, na.rm=TRUE)) {
        raw_data$general[species.index] <- 0
      }
    }
  }

  return(raw_data)
}
