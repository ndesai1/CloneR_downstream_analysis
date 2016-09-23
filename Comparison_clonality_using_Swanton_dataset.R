source("~/Lavoro/CloneR/functions.R")

load("Literature/STM_data/McGrhanam.Rdata")
x = read.mutation.file(file="Example/ACC_mutations.no_filter_freq.72.tsv")
# library(data.table)
# cR.swanton = fread("Results/PANCANCER_SWANTON_1741.Job18.2016-06-10.10_24_32/clonality.tsv")
# swanton = read.delim("/Volumes/temp2/CloneR/Results/PANCANCER_SWANTON_1741.Job18.2016-06-10.10_24_32/clonality.tsv")
# Remove (as Swanton) mutations on autosomal chromosome
# modified version from Cloner to get cell between 01

get.tc.correction.somatic = function( obs, tc, CNt, CNn=2){
  obs = as.numeric(obs)
  tc = as.numeric(tc)
  CNt = as.numeric(CNt)
  CNn = as.numeric(CNn)
  return( min(
    
    obs * ( 1 + (  ( CNn*(1-tc) )/( CNt * tc) ) )  ,
    
    1))
}


swanton$cR.freq.tc = NA
swanton$cR.freq.tc = apply(swanton, 1, function(x)  get.tc.correction.somatic(  x[13], x[2], x[16], x[14] ))
swanton$cR.clonality =  swanton$cR.freq.tc*2
swanton$cR.clonality[which(swanton$cR.clonality>1)] = 1

# corr.test CR clonality v. CCF
with(swanton, cor.test( cR.clonality, absolute.ccf ))
# with(swanton, scatter.smooth( cR.clonality, absolute.ccf ))

f = list.files("/Volumes/temp2/CloneR/Pyclone/Results", full.names = T)
library(data.table)
options(datatable.fread.datatable=FALSE)
f = lapply(f,  fread, stringsAsFactors=F)
l = which(sapply(f, nrow)>0)
f = f[l]
f = lapply(f, function(x) { colnames(x)=c("mutation_id", "clonality","std", "cluster_id" ); return(x) })
f = do.call(rbind, f)

id = match(swanton$mutation_id, f$mutation_id)

swanton$pyclone.clonality  = f$clonality[id]
swanton$pyclone.std        = f$std[id]
swanton$pyclone.cluster.id = f$cluster_id[id]

with(swanton, cor.test( cR.clonality, pyclone.clonality ))
with(swanton, cor.test( absolute.ccf, pyclone.clonality ))

save(swanton,file="Rdata/Comparison_clonality_Swanton_dataset.Rdata")
