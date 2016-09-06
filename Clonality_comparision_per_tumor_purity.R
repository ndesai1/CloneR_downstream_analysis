# -------------------------------------------------------------------------------------------
# PANCANCER SWANTON 1741 ASCAT TC
cancer_types = c(
  'ACC','BLCA','BRCA','CESC','CHOL','COAD','ESCA','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','UCEC','UCS','UVM'
)

f  = list.files("Results/", "clonality.tsv", recursive=T, full.names = T)
f = f[-grep("PANCAN",f)]
f = f[-grep("GAR",f)]
c.tcga = lapply(f, read.delim)
c.tcga = do.call(rbind,c.tcga)

c.tcga = c.tcga[ which(c.tcga$sample%in%samples_in_common), ]
save(c.tcga,file="Results/clonality_TCGA.1582.Rdata")

# -------------------------------------------------------------------------------------------
# PANCANCER SWANTON 1741 ASCAT TC

c.abs = read.delim("Results/PANCANCER_ABSOLUTE_2271.Job17.2016-06-08.17_33_02/clonality.tsv")
c.abs = c.abs[ which(c.abs$sample%in%samples_in_common), ]
save(c.abs,file="Results/clonality_ABSOLUTE.1582.Rdata")

# -------------------------------------------------------------------------------------------
# PANCANCER SWANTON 1741 ASCAT TC

c.swanton = read.delim("Results/PANCANCER_SWANTON_1741.Job18.2016-06-10.10_24_32/clonality.tsv")
c.swanton = c.swanton[ which(c.swanton$sample%in%samples_in_common), ]
save(c.swanton,file="Results/clonality_ASCAT.1582.Rdata")

# -------------------------------------------------------------------------------------------

load("Results/clonality_TCGA.1582.Rdata")
load("Results/clonality_ABSOLUTE.1582.Rdata")
load("Results/clonality_ASCAT.1582.Rdata")

# samples_in_common = tmp$sample
# save(samples_in_common,file="Results/1582_samples.Rdata")

c.abs$key     = paste0(c.abs$sample,".",c.abs$id)
c.tcga$key    = paste0(c.tcga$sample,".",c.tcga$id)
c.swanton$key = paste0(c.swanton$sample,".",c.swanton$id)

c.tcga     = c.tcga[which(!duplicated(c.tcga$key)),]
c.abs      = c.abs[which(!duplicated(c.abs$key)),]
c.swanton  = c.swanton[which(!duplicated(c.swanton$key)),]

id = intersect(c.abs$key,c.tcga$key)
id = intersect(id, c.swanton$key)
x = c.tcga[ which( c.tcga$key%in%id), c('key','cell') ]
x$cell2 = c.abs[ match( c.tcga$key, c.abs$key), 'cell']
x$cell3 = c.swanton[ match( c.tcga$key, c.swanton$key), 'cell']

corr=cor.test(x$cell, x$cell2)
ggplot(x, aes(x=cell,y=cell2))+geom_point(size=0.1)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("TCGA")+ylab("ABSOLUTE")+
  ggtitle(paste0("1582 samples, CloneR TCGA vs ABSOLUTE\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))

corr=cor.test(x$cell, x$cell3)
ggplot(x, aes(x=cell,y=cell3))+geom_point(size=0.1)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("TCGA")+ylab("ASCAT")+
  ggtitle(paste0("1582 samples, CloneR TCGA vs ASCAT\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))

corr=cor.test(x$cell2, x$cell3)
ggplot(x, aes(x=cell2,y=cell3))+geom_point(size=0.05)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("ABSOLUTE")+ylab("ASCAT")+
  ggtitle(paste0("1582 samples, CloneR ABSOLUTE vs ASCAT\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))



