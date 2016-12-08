# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 08/12/2016
# Desc: import and merge the 3m and 12m datasets


df3m = read.csv('Data_external/3m micro dataset with meta data for Umar 061216 (1).txt', sep='\t', stringsAsFactors = F)
dfLong = read.csv('Data_external/Umar longitudinal dataset IMs epi and bloods 071216.txt', sep='\t', stringsAsFactors = F)

rownames(df3m) = df3m$samples
rownames(dfLong) = dfLong$samples

## select the unique columns for merging
cn = unique(c(colnames(df3m), colnames(dfLong)))
# select columns from 3 month table
cn3m = cn[cn %in% colnames(df3m)]
# select remaining columns not present in 3m table
cnLong = cn[!(cn %in% colnames(df3m))]

# sanity check - if we cover all columns
table(cn %in% c(cn3m, cnLong))

# sample ids are unique so select these
samples = as.character(unique(c(df3m$samples, dfLong$samples)))

# fill in the 3month data first and select columns from long and 3m datasets
# match these samples with the 3m dataset as they should be present at least in 3m
i = match(samples, rownames(df3m))

lSamples.3m = lapply(samples[na.omit(i)], function(x){
  # copy columns from 3m data set and remaining extra columns, if any, from long data
  df = cbind(df3m[x,cn3m], dfLong[x,cnLong])
})

# convert to a data frame 
dfMerge.3m = do.call(rbind, lSamples.3m)

# repeat the process but this time we go through long data set
# select columns present in long data set from cn
cnLong2 = cn[cn %in% colnames(dfLong)]
# select remainder which will be in 3m dataset
cn3m2 = cn[!(cn %in% colnames(dfLong))] 

# use samples that did not match the 3m data, i.e. they are only present in long data
lSamples.long = lapply(samples[-(na.omit(i))], function(x){
  # copy columns that are present in long data, and any extra ones from 3m
  df = cbind(dfLong[x,cnLong2], df3m[x,cn3m2])
})

dfMerge.long = do.call(rbind, lSamples.long)

# order the columns in the same order
cn = unique(c('id', 'samples', 'age', cn))
rownames(dfMerge.3m) = dfMerge.3m$samples
rownames(dfMerge.long) = dfMerge.long$samples

dfMerge = rbind(dfMerge.3m[,cn], dfMerge.long[,cn])

write.csv(dfMerge, file='Data_external/merged_dataset.csv')

