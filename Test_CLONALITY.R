load("Rdata/Comparison_clonality_Swanton_dataset.Rdata")
library(plyr)

# test case COAD
ct = "COAD"

# mutations from SWANTON =================

mut = subset(swanton, cancer_type==ct)

#CNVs from SWANTON =================

load("Literature/STM_data/Input/tcga.coad.seg.hg19.rdata")
ascat = seg.mat.copy.list$segments[,c('SampleID','Chr','Start','End','cn','nA','nB')]

# CVN from TCGA =================

load("~/FC/DB/TCGA/01_03_2015/COAD/Tumor/CNV/COAD_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_GainLoss_BLAT_overlaps_nofilters.Rdata")

sc = lapply(somatic_cnv, function(x) unique(x[,c('Patient','Chromosome','Start','End','CNV_type','Copy_number','Segment_Mean')]))
to_exclude=which(sapply(somatic_cnv, function(x) sum(is.na(x$Chromosome)))==1)
if(length(to_exclude)>0) sc = sc[-to_exclude]
sc = do.call("rbind", sc)
# *********** CN RAW **********************
sc$Copy_number_raw = 2^sc$Segment_Mean * 2 
# *****************************************
sc=sc[,c("Patient" , 'Chromosome','Start','End','CNV_type','Copy_number','Copy_number_raw')]

# sample in common between from tcga 255 =================
l.ascat = dlply(ascat, .(SampleID))
l.tcga  = dlply(sc, .(Patient))
l.mut = dlply(mut, .(Patient))

sampleID = intersect(intersect(names(l.tcga), names(l.mut)),names(l.ascat))
print(length(sampleID))

l.tcga = l.tcga[sampleID]
l.ascat = l.ascat[sampleID]
l.mut = l.mut[sampleID]

# Assign CN_RAW TCGA to SWANTON =================
set_CN_raw_from_TCGA = function(m, c=NULL, chr=chr_levels){
  require(GenomicRanges)
  if(!is.null(c) & !is.null(m)){
    m$CN_raw = NA
    im  = with(m, GRanges(Chromosome, IRanges(Position,Position) )); seqlevels(im, force=T)=chr
    ic  = with(c, GRanges(Chromosome, IRanges(Start,End) )); seqlevels(ic, force=T)=chr
    ov  = findOverlaps(im,ic)
    if(length(ov)>0){
      m$CN_raw[queryHits(ov)] = c$Copy_number_raw[subjectHits(ov)]
    }else{
    }
  }
  return(m)
}
chr_levels = paste0('chr',c(as.character(1:22),'X','Y'))

l.mut = mapply( set_CN_raw_from_TCGA, l.mut, l.tcga, SIMPLIFY = F)


# MERGE CNV DATA FROM ASCAT AND TCGA ================

# CNV data from TCGA from patient TCGA-A6-2671
head(l.tcga[["TCGA-A6-2671"]])
# CNV data from ASCAT for patient TCGA-A6-2671
head(l.ascat[["TCGA-A6-2671"]])
# SNV data from SWANTON for patient TCGA-A6-2671
head(l.mut[["TCGA-A6-2671"]][,c('Chromosome',"Position","ref_counts","var_counts","obs.VAF","normal_cn","minor_cn",'major_cn',"CN_raw")])

save.image(file="~/Lavoro/CloneR_downstream_analysis/Rdata/TEST_CLONALITY_COAD.Rdata")

