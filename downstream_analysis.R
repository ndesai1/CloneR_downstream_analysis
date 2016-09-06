library("TCGAbiolinks")
cancer_types = c(
  'ACC','BLCA','BRCA','CESC','CHOL','COAD','ESCA','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','UCEC','UCS','UVM'
)

# COVERAGE Comparison ==============
cov=list()
for(i in cancer_types){
  print(i)
  load(paste0("/Volumes/FC/DB/TCGA/01_03_2015/",i,"/Tumor/Somatic_Mutations/",i,"_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata"))
  x=lapply(somatic_mutations, function(x){
    x$cov1=as.numeric(x$TTotCov)
    x$cov2=as.numeric(x$t_alt_count)+as.numeric(x$t_ref_count)
    x$cov3=as.numeric(x$tumor_ref_reads)+as.numeric(x$tumor_var_reads)
    x$cov4=as.numeric(x$Tumor_Alt_Count)+as.numeric(x$Tumor_Mut_Freq)
    x$mean_cov=apply(x[,c("cov1","cov2","cov3","cov4")],1,mean, na.rm=T)
    x=unrowname(x)
    x[,c("Patient","Center",'key',"cov1","cov2","cov3","cov4","mean_cov")]
  })

  x = mapply(function(x,y){x$cancer=y; return(x)},x, rep(i,length(x)), SIMPLIFY = F  )

  cov[[i]]=x
}
# 
# x = lapply(cov, function(x) do.call(c,lapply(x, function(y) y$mean_cov)))
# x = lapply(x, function(y) summary(y)[c('Min.','1st Qu.','Median','Mean','3rd Qu.','Max.')])
# x = do.call(rbind,x)
# coverage_cancer_type=x
# save(coverage_cancer_type, file="Results/converage_cancer_types.Rdata")

load("Results/converage_cancer_types.Rdata")
load("Results/pancancer_clone_composition.27.Rdata")

# poly = subset(pcc, composition=="P")
# pdf(file="Results/scatter_plot_Median_cov_Perc_polyclonal.pdf",w=6,h=6)
# ggplot(poly, aes(x=Median,y=polyclonal))+geom_point()+geom_smooth(method="lm", se=F)+theme_bw()+ylab("Polyclonal(%)")+xlab("Median coverage")
# dev.off()

y = ddply(pcc, .(cancer_type), function(a) table(a$composition)[c("M","B","P")])
y[,2:4] = y[,2:4]/ apply(y[,2:4], 1, sum, na.rm=T)*100
y$cancer_type = sapply( strsplit( y$cancer_type," "), function(x) x[1])

cancer_type = cbind( y, coverage_cancer_type[match(y$cancer_type, rownames(coverage_cancer_type)),])

f = list.files("Example/", pattern="mutations", full.names = T)
f = f[which(!f%in%c("Example//PANCANCER_ABSOLUTE_mutations.no_filter_freq.2271.tsv","Example//PANCANCER_SWANTON_mutations.no_filter_freq.1741.tsv"))]

n = lapply(f, read.delim)
n = do.call(rbind,lapply(n, nrow))
f = gsub("Example//","", f)
f = sapply(strsplit(f, split="\\_"), function(x) x[1] )
rownames(n) = f

cancer_type$mutations = n[ match(y$cancer_type, rownames(n)), 1 ]

f = list.files("Example/", pattern="cnvs", full.names = T)
f = f[which(!f%in%c("Example//PANCANCER_ABSOLUTE_cnvs.2271.tsv","Example//PANCANCER_SWANTON_cnvs.1741.tsv"))]
n = lapply(f, read.delim)
n = do.call(rbind,lapply(n, nrow))
f = gsub("Example//","", f)
f = sapply(strsplit(f, split="\\_"), function(x) x[1] )
rownames(n) = f

cancer_type$cnvs = n[ match(y$cancer_type, rownames(n)), 1 ]

save(cancer_type, file="Results/pancancer_clone_composition.27.summary_by_cancer_type.Rdata")

pdf(file="Results/scatter_plot_Median_cov_Perc_polyclonal_per_cancer_type.pdf",w=6,h=6)
ggplot(cancer_type, aes(x=Median,y=P))+geom_point()+geom_smooth(method="lm", se=F)+theme_bw()+ylab("Polyclonal(%)")+xlab("Median coverage")
dev.off()

# pdf(file="Results/scatter_plot_Mutations_cov_Perc_polyclonal_per_cancer_type.pdf",w=6,h=6)
# ggplot(cancer_type, aes(x=mutations,y=P))+geom_point()+geom_smooth(method="lm", se=F)+theme_bw()+ylab("Polyclonal(%)")+xlab("Mutations")
# dev.off()
# 
# pdf(file="Results/scatter_plot_cnvs_cov_Perc_polyclonal_per_cancer_type.pdf",w=6,h=6)
# ggplot(cancer_type, aes(x=cnvs,y=P))+geom_point()+geom_smooth(method="lm", se=F)+theme_bw()+ylab("Polyclonal(%)")+xlab("CNVs")
# dev.off()

f = list.files("Example/", pattern="samples", full.names = T)
f = f[which(!f%in%c("Example//PANCANCER_ABSOLUTE_cnvs.2271.tsv","Example//PANCANCER_SWANTON_cnvs.1741.tsv"))]
n = do.call(rbind,lapply(f, read.delim, header=F))
colnames(n) = c('sample','male','tc')

pcc = cbind(pcc, n[match(pcc$sample, n$sample), 2:3] )

x = do.call(rbind,lapply(cov, function(x) do.call(rbind,lapply(x, function(y) c('sample' = y$Patient[1], summary(y$mean_cov)[c('Min.','1st Qu.','Median','Mean','3rd Qu.','Max.')])))))
x = as.data.frame(x)
for(i in 2:ncol(x)) x[,i] = as.numeric(x[,i])
x=unrowname(x)

pcc = cbind(pcc, x[match(pcc$sample, x$sample), 2:ncol(x)] )
pcc=pcc[,c(1,10:18,2:9)]

pcc = unrowname(pcc)

# add ABSOLUTE
# -------------
absolute = read.delim("Purity/pancan12.sample_info.txt")
absolute$patient = substr(absolute$tcga_id,1,12)
colnames(absolute) = paste0("ABSOLUTE_", colnames(absolute))

pcc = cbind(pcc, absolute[match(pcc$sample, absolute$ABSOLUTE_patient),5:10] )

a = read.delim("/Volumes/ceredam/CloneR/Results/PANCANCER_ABSOLUTE_2271.Job17.2016-06-08.17_33_02/composition.tsv")
colnames(a)=paste0("ABSOLUTE_2271_",colnames(a))

pcc = cbind(pcc, a[match(pcc$sample, a$ABSOLUTE_2271_sample),2:9] )


# add PANCANCER 1157
# -------------------
andor = read.xlsx("Literature/Andor_2015_Pan-cancer_analysis_results .xlsx", 1)
colnames(andor) = paste0("ANDOR_1157_",colnames(andor) )
pcc = cbind(pcc, andor[match(pcc$sample, andor$ANDOR_1157_TCGA.ID),2:9] )

# with(subset(pcc, !is.na(ANDOR_1157_CloneNumber.PurityNormalized.)), cor.test( as.numeric(composition), ANDOR_1157_CloneNumber.PurityNormalized., meth="spearman" ))
# with(subset(pcc, !is.na(ANDOR_1157_CloneNumber.PurityNormalized.)), scatter.smooth( as.numeric(composition), ANDOR_1157_CloneNumber.PurityNormalized., meth="spearman" ))


# add SWANTON 1741
# -------------------
load("Literature/STM_data/McGrhanam.Rdata")

tmp = unique(swanton[,c('Patient','TCGA.purity')])
colnames(tmp)[2] = "ASCAT.purity"
colnames(tmp) = paste0("SWANTON_1741_",colnames(tmp))

pcc$SWANTON_1741_ASCAT.purity = tmp[match(pcc$sample, tmp$SWANTON_1741_Patient),2] 

a = read.delim("/Volumes/ceredam/CloneR/Results/PANCANCER_SWANTON_1741.Job18.2016-06-10.10_24_32/composition.tsv")
colnames(a)=paste0("SWANTON_1741_",colnames(a))

pcc = cbind(pcc, a[match(pcc$sample, a$SWANTON_1741_sample),2:9] )

pcc$composition=factor(pcc$composition, levels=c("M","B","P"))
pcc$SWANTON_1741_composition=factor(pcc$SWANTON_1741_composition, levels=c("M","B","P"))
pcc$ABSOLUTE_2271_composition=factor(pcc$ABSOLUTE_2271_composition, levels=c("M","B","P"))

save(pcc, file="Results/pancancer_clone_composition.27.Rdata")


#FIGURES TUMOUR CONTENT ===========

pdf(file="Results/boxplot_composition_tumour_content_per_sample.pdf", h=6, w=6)
ggplot(pcc, aes(x=composition, y=tc))+geom_boxplot()+theme_bw() +ylab("")+ggtitle("6331 TCGA tumour Content")                      
dev.off()

pdf(file="Results/boxplot_composition_tumour_content_per_sample.SWANTON_1741.pdf", h=6, w=6)
ggplot(subset(pcc, !is.na(SWANTON_1741_n)), aes(x=composition, y=tc))+geom_boxplot()+theme_bw() +xlab("")+ggtitle("1741 TCGA tumour Content")                      
dev.off()

pdf(file="Results/boxplot_SWANTON_composition_ASCAT_purity_per_sample.SWANTON_1741.pdf", h=6, w=6)
ggplot(subset(pcc, !is.na(SWANTON_1741_n)), aes(x=SWANTON_1741_composition, y=SWANTON_1741_ASCAT.purity))+geom_boxplot()+theme_bw() +xlab("")+ylab("ASCAT purity")+ggtitle("1741 ASCAT purity")                      
dev.off()

pdf(file="Results/boxplot_composition_tumour_content_per_sample.ABSOLUTE_2271.pdf", h=6, w=6)
ggplot(subset(pcc, !is.na(ABSOLUTE_2271_composition)), aes(x=composition, y=tc))+geom_boxplot()+theme_bw() +xlab("")+ggtitle("2271 TCGA tumour Content")                      
dev.off()

pdf(file="Results/boxplot_ABSOLUTE_composition_ABSOLUTE_purity_per_sample.ABSOLUTE_2271.pdf", h=6, w=6)
ggplot(subset(pcc, !is.na(ABSOLUTE_2271_composition)), aes(x=ABSOLUTE_2271_composition, y=ABSOLUTE_abs_purity))+geom_boxplot()+theme_bw() +xlab("")+ylab("ABSOLUTE purity")+ggtitle("2271 ABSOLUTE purity")                      
dev.off()

tmp = subset(pcc, !is.na(ABSOLUTE_2271_composition) & (!is.na(SWANTON_1741_composition)) )

pdf(file="Results/boxplot_composition_TCGA_purity_per_sample.TCGA_1582.pdf", h=6, w=6)
ggplot(tmp , aes(x=composition, y=tc))+geom_boxplot()+theme_bw() +xlab("")+ggtitle("1582 TCGA purity")   +ylab("TCGA purity")                   
dev.off()

pdf(file="Results/boxplot_composition_ABSOLUTE_purity_per_sample.ASCAT_ABSOLUTE_1582.pdf", h=6, w=6)
ggplot(tmp , aes(x=ABSOLUTE_2271_composition, y=ABSOLUTE_abs_purity))+geom_boxplot()+theme_bw() +xlab("")+ggtitle("1582 ABSOLUTE purity")   +ylab("ABSOLUTE purity")                   
dev.off()

pdf(file="Results/boxplot_composition_ASCAT_purity_per_sample.ASCAT_ABSOLUTE_1582.pdf", h=6, w=6)
ggplot(tmp, aes(x=SWANTON_1741_composition, y=SWANTON_1741_ASCAT.purity))+geom_boxplot()+theme_bw() +xlab("")+ylab("ASCAT purity")+ggtitle("1582 ASCAT purity")                      
dev.off()

pdf(file="Results/correlation_TCGA_tumour_content_vs_ASCAT_purity.1741.pdf", h=6, w=6)
corr=with(subset(pcc, !is.na(SWANTON_1741_n)),cor.test(tc, SWANTON_1741_ASCAT.purity))
ggplot(subset(pcc, !is.na(SWANTON_1741_n)), aes(x=tc,y=SWANTON_1741_ASCAT.purity))+geom_point(size=0.5)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("TCGA tumour content")+ylab("ASCAT purity")+
  ggtitle(paste0("1741 samples\nTCGA tc vs. ASCAT purity\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))
dev.off()

pdf(file="Results/correlation_TCGA_tumour_content_vs_ABSOLUTE_purity.2271.pdf", h=6, w=6)
corr=with(subset(pcc, !is.na(ABSOLUTE_2271_composition)),cor.test(tc, ABSOLUTE_abs_purity))
ggplot(subset(pcc, !is.na(ABSOLUTE_2271_composition)), aes(x=tc,y=ABSOLUTE_abs_purity))+geom_point(size=0.5)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("TCGA tumour content")+ylab("ABSOLUTE purity")+
  ggtitle(paste0("2271 samples\nTCGA tc vs. ABSOLUTE purity\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))
dev.off()

pdf(file="Results/correlation_ABSOLUTE_purity_vs_ASCAT_purity.1741.pdf", h=6, w=6)
corr=with(subset(pcc, !is.na(SWANTON_1741_n)),cor.test(ABSOLUTE_abs_purity, SWANTON_1741_ASCAT.purity))
ggplot(subset(pcc, !is.na(SWANTON_1741_n)), aes(x=ABSOLUTE_abs_purity,y=SWANTON_1741_ASCAT.purity))+geom_point(size=0.5)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("ABSOLUTE purity")+ylab("ASCAT purity")+
  ggtitle(paste0("1741 samples\nABSOLUTE purity vs. ASCAT purity\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))
dev.off()

#FIGURES COMPOSITON =========

pdf(file="Results/correlation_composition_TCGA_tumour_content_vs_composition_ASCAT_purity.1741.pdf", h=6, w=6)
corr=with(subset(pcc, !is.na(SWANTON_1741_n)),cor.test(as.numeric(composition), as.numeric(SWANTON_1741_composition)))
ggplot(subset(pcc, !is.na(SWANTON_1741_n)), aes(x=(composition),y=(SWANTON_1741_composition)))+geom_point(size=0.5)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("compostion TCGA tc")+ylab("composition ASCAT purity")+
  ggtitle(paste0("1741 samples\n composition TCGA tc vs. ASCAT purity\n Pearson R=", round(corr$estimate,2), " pv=", corr$p.value))
dev.off()

with(subset(pcc, !is.na(SWANTON_1741_n)), table(composition, SWANTON_1741_composition))
with(subset(pcc, !is.na(ABSOLUTE_2271_composition)),cor.test(as.numeric(composition), as.numeric(ABSOLUTE_2271_composition)))
with(subset(pcc, !is.na(ABSOLUTE_2271_composition)), table(composition, ABSOLUTE_2271_composition))


tmp = subset(pcc, !is.na(SWANTON_1741_composition) & !is.na(ABSOLUTE_2271_composition))
pcc.1582 = tmp 
save(pcc.1582, file="Results/pancancer_clone_composition.27.1582_samples.Rdata")

with(tmp,cor.test(as.numeric(ABSOLUTE_2271_composition), as.numeric(SWANTON_1741_composition)))
with(tmp,cor.test(as.numeric(composition), as.numeric(SWANTON_1741_composition)))
with(tmp,cor.test( as.numeric(composition), as.numeric(ABSOLUTE_2271_composition)))

with(tmp, table(ABSOLUTE_2271_composition, SWANTON_1741_composition))
with(tmp, table(composition, SWANTON_1741_composition))
with(tmp, table(composition, ABSOLUTE_2271_composition))


tmp$jump = with(tmp, paste0(ABSOLUTE_2271_composition,"(a) - ", SWANTON_1741_composition,"(s)"))
tmp$jump = factor(tmp$jump, levels = c( "M(a) - M(s)", "B(a) - B(s)", "P(a) - P(s)", "M(a) - B(s)", "M(a) - P(s)", "B(a) - M(s)", "P(a) - M(s)", "P(a) - B(s)"  ))

pdf(file="Results/correlation_composition_ABSOLUTE_vs_ASCAT_purity.by_jump.1582.pdf", h=6, w=9)
ggplot(tmp, aes(x=ABSOLUTE_abs_purity,y=SWANTON_1741_ASCAT.purity, shape=jump, color=jump, size=jump, fill=jump))+geom_point(size=2.5)+
  scale_color_manual(values = c("grey70","grey70","grey70","red","blue","#998ec3","#dd1c77","#008837"))+
  scale_fill_manual(values = c("grey70","grey70","grey70","red","blue","#998ec3","#dd1c77","#008837"))+
  scale_shape_manual(values=c(0:2, 21:25))+geom_abline(slope=1,intercept = 0, lty='dashed')+
  scale_size_manual(values=c(rep(1,3),rep(3,5)))+
  theme_bw()+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,1))+
  ylab("ASCAT purity")+xlab("ABSOLUTE purity")+theme(panel.grid.minor=element_blank())
dev.off()


tmp$jump = with(tmp, paste0(composition,"(a) - ", SWANTON_1741_composition,"(s)"))
tmp$jump = factor(tmp$jump, levels = c( "M(a) - M(s)", "B(a) - B(s)", "P(a) - P(s)", "M(a) - B(s)", "M(a) - P(s)", "B(a) - M(s)", "P(a) - M(s)", "P(a) - B(s)"  ))

pdf(file="Results/correlation_composition_TCGA_vs_ASCAT_purity.by_jump.1582.pdf", h=6, w=9)
ggplot(tmp, aes(x=ABSOLUTE_abs_purity,y=SWANTON_1741_ASCAT.purity, shape=jump, color=jump, size=jump, fill=jump))+geom_point(size=2.5)+
  scale_color_manual(values = c("grey70","grey70","grey70","red","blue","#998ec3","#dd1c77","#008837"))+
  scale_fill_manual(values = c("grey70","grey70","grey70","red","blue","#998ec3","#dd1c77","#008837"))+
  scale_shape_manual(values=c(0:2, 21:25))+geom_abline(slope=1,intercept = 0, lty='dashed')+
  scale_size_manual(values=c(rep(1,3),rep(3,5)))+
  theme_bw()+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,1))+
  ylab("ASCAT purity")+xlab("TCGA purity")+theme(panel.grid.minor=element_blank())
dev.off()

tmp$jump = with(tmp, paste0(composition,"(a) - ", ABSOLUTE_2271_composition,"(s)"))
tmp$jump = factor(tmp$jump, levels = c( "M(a) - M(s)", "B(a) - B(s)", "P(a) - P(s)", "M(a) - B(s)", "M(a) - P(s)", "B(a) - M(s)", "P(a) - M(s)", "P(a) - B(s)"  ))

pdf(file="Results/correlation_composition_TCGA_vs_ABS_purity.by_jump.1582.pdf", h=6, w=9)
ggplot(tmp, aes(x=ABSOLUTE_abs_purity,y=SWANTON_1741_ASCAT.purity, shape=jump, color=jump, size=jump, fill=jump))+geom_point(size=2.5)+
  scale_color_manual(values = c("grey70","grey70","grey70","red","blue","#998ec3","#dd1c77","#008837"))+
  scale_fill_manual(values = c("grey70","grey70","grey70","red","blue","#998ec3","#dd1c77","#008837"))+
  scale_shape_manual(values=c(0:2, 21:25))+geom_abline(slope=1,intercept = 0, lty='dashed')+
  scale_size_manual(values=c(rep(1,3),rep(3,5)))+
  theme_bw()+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,1))+
  ylab("ABSOLUTE purity")+xlab("TCGA purity")+theme(panel.grid.minor=element_blank())
dev.off()


pdf(file="Results/correlation_composition_ABSOLUTE_vs_ASCAT_purity.by_jump.no_same.1582.pdf", h=6, w=9)
ggplot(subset(tmp, !jump%in%c("M(a) - M(s)", "B(a) - B(s)", "P(a) - P(s)")),
       aes(x=ABSOLUTE_abs_purity,y=SWANTON_1741_ASCAT.purity, shape=jump, color=jump, size=jump, fill=jump))+geom_point(size=2.5)+
  scale_color_manual(values = c("red","blue","#998ec3","#dd1c77","#008837"))+
  scale_fill_manual(values = c("red","blue","#998ec3","#dd1c77","#008837"))+
  scale_shape_manual(values=c( 21:25))+geom_abline(slope=1,intercept = 0, lty='dashed')+
  scale_size_manual(values=c(rep(3,5)))+
  theme_bw()+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,1))+
  ylab("ASCAT purity")+xlab("ABSOLUTE purity")+theme(panel.grid.minor=element_blank())
dev.off()

pdf(file="Results/correlation_TCGA_purity_vs_ASCAT_purity.ABSOLUTE_purity_1.pdf", h=6, w=6)
corr=with(subset(tmp, ABSOLUTE_abs_purity==1),cor.test(tc, SWANTON_1741_ASCAT.purity))
ggplot(subset(tmp, ABSOLUTE_abs_purity==1), aes(x=tc,y=SWANTON_1741_ASCAT.purity))+geom_point(size=0.5)+geom_smooth(method="lm", se=F)+theme_bw()+xlab("TCGA purity")+ylab("ASCAT purity")+
  ggtitle(paste0("64 samples ABSOLUTE purity=1\nTCGA purity vs. ASCAT purity\n Pearson R=", round(corr$estimate,2), " pv=", round(corr$p.value,2)))
dev.off()


#==== PANCANCER

x = subset(pcc, cancer_type=="CHOL (35)")
# z = z[,c('sample',)]
x = x[,c('sample',"composition",'monoclonal',  'biclonal', 'polyclonal')]
x$composition = factor(as.character(x$composition), levels=c("P","B","M") )
# x = x[order((x$composition),x$monoclonal,x$biclonal,x$polyclonal, decreasing = T),]
x = x[order(x$monoclonal,x$biclonal,x$polyclonal, decreasing = T),]
# ids = (x$sample)
y = melt(x, id.var=c('sample',"composition"))

y$variable = factor( as.character(y$variable), levels=c( 'polyclonal' , 'biclonal' , 'monoclonal'))

# y$sample = factor(y$sample, levels=ids)
cl = color_clone_composition
names(cl)=c("monoclonal",'biclonal','polyclonal')


ggplot(y, aes(x=sample, y=value,fill=variable))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values=cl)

plot_ly(y, x=sample, y=value, color=(variable), colors=rev(cl), type='bar') %>%
        layout(barmode = 'stack',
               yaxis=list(title='', tickfont=list(family = "Arial, sans-serif")),
               xaxis=list(title='Alterations (%)', titlefont=list(family = "Arial, sans-serif"), tickfont=list(family = "Arial, sans-serif"))
        )

y$sample = factor(y$sample, levels=x$sample)
ggplot(y, aes(x=sample, y=value,fill=variable))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values=cl)+theme(axis.text.x=element_text(angle=90))





# FIGURES VENN ================
library(gplots)

venn(list("CloneR"=pcc$sample, "SWANTON"=unique(swanton$Patient), "ABSOLUTE" = unique(absolute$ABSOLUTE_patient )))

pdf(file="Results/venn_composition_CLONER_SWANTON_ABSOLUTE.pdf", h=6, w=6)
venn(list("CloneR"=pcc$sample, "SWANTON"=subset(pcc, !is.na(SWANTON_1741_composition))$sample, "ABSOLUTE" = subset(pcc, !is.na(ABSOLUTE_2271_composition))$sample))
dev.off()

