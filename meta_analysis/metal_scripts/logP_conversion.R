#####
# Tanya Major
# June 2019
# University of Otago
#####

## Usage: Rscript logP_conversion.R filename
### p-value column MUST be named P-value / P.value

library(data.table)
library(dplyr)

args <- commandArgs(TRUE)

filename <- args[1]

# load data, specify P.value column must be a character to prevent 0e+0 replacing numbers < e-300
data <- as.data.frame(fread(file = filename))

# Function to make logP for e- numbers:
elog10 <- function(p) {
	if (is.character(p) & grepl('e-', p)) {
		split_p <- base::strsplit(p, split = "e-")
		tmp <- unlist(split_p)
		res <- as.numeric(tmp[2]) - log10(as.numeric(tmp[1]))
	} else {
		res <- -log10(as.numeric(p))
	}
	return(res)
}

data$log10P = lapply(data$`P-value`, elog10)

# save file
fwrite(data, file = paste0(filename, ".logP"), quote = F, sep = "\t", na = "NA", col.names = T, row.names = F)

