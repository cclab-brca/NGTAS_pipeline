library(VariantAnnotation)
args <- commandArgs()

# import vcfs
file.hc<-paste0(args[2],args[3],".hc.vcf",sep="")
param <- ScanVcfParam()
vcf.hc <- readVcf(file.hc, "hg19", param)
file.sort<-paste0(args[2],args[3],".sort.vcf",sep="")
param <- ScanVcfParam()
vcf.sort <- readVcf(file.sort, "hg19", param)

# compute average VAF
vaf.all<-numeric()
for(i in 1:ncol(vcf.hc)) {
	ad<-data.frame(geno(vcf.hc)$AD[,i])
	ref.dp<-ad[1,]
	alt.dp<-ad[2,]
	vaf<-alt.dp/(alt.dp+ref.dp)
	vaf.all<-rbind(vaf.all,vaf)
}	
rownames(vaf.all)<-colnames(vcf.hc)
vaf.t<-as.matrix(vaf.all[args[4:6],])
vaf.t.avg<-colMeans(vaf.t,na.rm=T)
vaf.n<-as.matrix(vaf.all[args[7:9],])
vaf.n.avg<-colMeans(vaf.n,na.rm=T)

# ratio
tn_ratio<-vaf.t.avg/vaf.n.avg

# C>T or G>A
mutID<-as.matrix(rownames(vcf.hc))
ref_alt<-as.matrix(apply(mutID,1,function(x){strsplit(x, "_")[[1]][2]}))
ref<-apply(ref_alt,1,function(x){strsplit(x, "/")[[1]][1]})
alt<-apply(ref_alt,1,function(x){strsplit(x, "/")[[1]][2]})
ct_ga<-(ref=="C"	&alt=="T")|(ref=="G"	&alt=="A")

# FILTER
if(args[10]=="FFPE")	{
	filter<-vaf.n.avg<as.double(args[11])&vaf.t.avg>as.double(args[12])&tn_ratio>as.double(args[13])&(ct_ga==FALSE | (ct_ga==T & vaf.t.avg>as.double(args[14])))
} else {
	filter<-vaf.n.avg<as.double(args[11])&vaf.t.avg>as.double(args[12])&tn_ratio>as.double(args[13])
}
vcf.f<-vcf.sort[filter]

if(nrow(vcf.f)>0) {
	# mutation list txt
	mut<-rownames(vcf.f)
	sample<-apply(as.matrix(args[3]),1,function(x){strsplit(x,"\\.")[[1]][1]})
	ampli<-apply(as.matrix(args[3]),1,function(x){strsplit(x,"\\.")[[1]][2]})
	called<-data.frame(mut,sample,ampli)
	write.table(called,file=paste0(args[2],args[3],".called.txt"),quote=F,sep="\t",row.names=F,col.names=F)

	# save vcf
	writeVcf(obj = vcf.f, filename = paste0(args[2],args[3],".sort.filt.vcf"))
}


#args<-c("","/tmp/mut_10786/","0872D2_20.PIK3CA_D0069_056F", "0872D2_20_a","0872D2_20_b","0872D2_20_c","0872N_5_a","0872N_5_b","0872N_5_c","FFPE","0.01","0.01","5","0.2")
	