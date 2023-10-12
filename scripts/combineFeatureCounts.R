#!/usr/bin/env Rscript

library(data.table)
library(foreach)

# Get vector of all files within the FEATURE_COUNTS directory ending with .txt
files <- list.files('FEATURE_COUNTS', pattern='*.txt$', full.names=TRUE)

# Iterate over file names, rbind all into one long table
dat.long <- foreach(file=files, .combine='rbind') %do% {
	# Define sample_id by extracting it from the file name
	sample_id <- strsplit(file, split='/|\\.')[[1]][2]

	# Read in the feature counts file, starting with the header line
	tmp <- fread(file, skip='Geneid')

	# Rename the column (filename -> count)
	col_ids <- colnames(tmp)[1:6]
	setnames(tmp, c(col_ids,'count'))

	# Add new column identifying the sample
	tmp[, sample := sample_id]
}

# Convert to feature count matrix (rows = genes, columns = samples)
dat.wide <- dcast(dat.long, Geneid~sample, value.var='count')

# Save wide form matrix
fwrite(dat.wide, file='countsMatrix.tsv', sep='\t')



# Add treatment columns to long table
dat.long[sample %like% 'Foliate', treatment := 'Foliate']
dat.long[sample %like% 'Organoid', treatment := 'Organoid']


# Calculate each gene's mean expression across 5 replicates (per treatment)
dat.means <- dat.long[, mean(count), by=list(Geneid, treatment)]
dat.means.wide <- dcast(dat.means, Geneid~treatment, value.var='V1')

# Plot log2(mean+1)
ggplot(dat.means.wide, aes(x=log2(1+Foliate),y=log2(1+Organoid))) + geom_point() +
	labs(x='log2(mean foliate counts + 1)', y='log2(mean organoid counts + 1)') +
	geom_abline(slope=1, intercept=0, color='gray', linetype='dashed') +
	xlim(0,20) + ylim(0,20)

