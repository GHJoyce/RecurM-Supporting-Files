# For running in R-Studio or Jupyter. 

### Use PropR to generate CLR data and identify hosts to investigate

## Set up the workspace --------------------------------------------------------------------------------------

# Set Working Directory
setwd("/srv/home/s4479877")

dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  
# add to the path
.libPaths(Sys.getenv("R_LIBS_USER"))  
#install.packages("propr")
#install.packages('ggdendro')

library("propr")
library(dplyr)
library(data.table)
library(ggplot2)

## Set up the DataTable ---------------------------------------------------------------------------------------

# Read in CoverM output
r = fread('Jupyter/reads/reads.tsv')

# OPTIONAL: IF READING IN MULTIPLE FILES
# Combine coverM output
r = merge(r1, r2, by="Genome", all=TRUE )
r = merge(r, r3, by="Genome", all=TRUE)
r = merge(r, r4, by="Genome", all=TRUE)
r = merge(r, r5, by="Genome", all=TRUE)

m = melt.data.table(r, id.vars='Genome', variable.name='Sample',value.name='value')

# OPTIONAL: IF COVERM OUT INCLUDES 2 MEASURES, SEPARATE THEM OUT
# Add a new column, measure, with a label for sample type.
m[grep("Trimmed", Sample), measure:= 'trimmed_mean']
m[grep("Relative", Sample), measure:= 'rel_abundance']

# Remove the excess text from sample names
m[, Sample := gsub('_.*','',Sample)]

# Extract just the trimmed mean/rel_abundance data and adjust the dataframe so its suitable for Propr

trim_df = m[measure =="trimmed_mean"]
trim_df2 = trim_df[Genome != "unmapped"]

# Change to matrix function
df2matrix <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}
# Deal with NAs and 0s

df[is.na(df)] <- 0
df2 = df[apply(df[,-1], 1, function(x) !all(x<400)),]
df2[1:3]
M_df <- df2matrix(as.data.frame(df2))
df3 = M_df + 0.001

## Use Propr. 

# Use Propr
rho <- propr(df3, metric = "rho", ivar = "clr", symmetrize = TRUE)
lr = rho@logratio

# Export CLR data (Use this for Find-Hosts.py)
       
write.csv(lr,'Jupyter/spreadsheets/clr_27-08-20.csv')

### Plot the CLR data using qplot ----------------------------------------------------------------------------------------


library(data.table)
library(ggplot2)


# Read in PropR output
df = fread("Jupyter/spreadsheets/clr_V2.csv")

m = melt.data.table(df, id.vars='V1', variable.name='sample',value.name='value')

colnames(m)[1] = "Sample"
colnames(m)[2] = "Genome"

## Identify the data of interest

# Group the MAGS of inetrest (example below). Identities are from gtdb-tk

plasmidA='CLUSTER_size_3_avlen_33206_avcov_215'

## Subset dataframe for plasmids/MAGs of interest
# PlasmidA

df_plasA_1 = m[Genome == "f.PSM7J193.6" | Genome==plasmidA]
df_plasA_1[1:3]

# OPTIONAL: Add together the abundances of the microbes in same taxonomic group
#df_plasA_2 = df_plasA_1 [Genome %in% b_vulgatus, .(value=sum(value), Genome ='b_vulgatus'), by =Sample]
#df_plasA_3 = rbind(df_plasA_2, df_plasA_1 [Genome==plasmidA, .(value, Genome, Sample)])

# Change cluster name to PlasmidA
df_plasA_1[Genome=='CLUSTER_size_3_avlen_33206_avcov_215', Genome := 'plasmidA']

# Rearrange datatable so that its suitable for plotting
m1 = dcast.data.table(df_plasA_1, Sample~Genome, value.var='value', fill=0)

m2 = melt.data.table(m1, id.vars=c('Sample','plasmidA'), variable.name='b. thetaiotaomicron', value.name='value')

# Add participantIDs to datatable
id_df = fread('Jupyter/spreadsheets/participant-sample.csv')
m3 = merge(m2, id_df, by = "Sample", allow.cartesian = TRUE)
m3[1:3]

## Plot data with qplot
# IF OPTIONAL STEP WAS DONE: plots will be side by side, by organism in group
# IF OPTIONAL STEP NOT DONE: one plot presented, abudances of all organisms combined.

plot = qplot(data=m3, plasmidA, value, colour = Participant)#+facet_wrap(~b_vulgatus)
plot # + theme(legend.position='none')
