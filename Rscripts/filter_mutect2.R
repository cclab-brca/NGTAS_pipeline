library(VariantAnnotation)

args <- commandArgs()
inputVCF<-paste0(args[2],args[3],".mut2.vcf",sep="")
# import
param <- ScanVcfParam()
vcf <- readVcf(inputVCF, "hg19", param)
if(nrow(vcf)>0) {
	#filter
	af<-geno(vcf)$AF
	tn_af<-af[,"TUMOR"]/af[,"NORMAL"]
	n_cov<-as.vector(as.matrix(data.frame(geno(vcf)$AD[,"NORMAL"])[1,]))+as.vector(as.matrix(data.frame(geno(vcf)$AD[,"NORMAL"])[2,]))
	t_cov<-as.vector(as.matrix(data.frame(geno(vcf)$AD[,"TUMOR"])[1,]))+as.vector(as.matrix(data.frame(geno(vcf)$AD[,"TUMOR"])[2,]))
	keep<-c("PASS",
		"alt_allele_in_normal",
		"alt_allele_in_normal;clustered_events",
		"alt_allele_in_normal;clustered_events;homologous_mapping_event",
		"alt_allele_in_normal;homologous_mapping_event",
		"clustered_events",
		"clustered_events;homologous_mapping_event",
		"clustered_events;homologous_mapping_event;multi_event_alt_allele_in_normal",
		"clustered_events;multi_event_alt_allele_in_normal",
		"homologous_mapping_event",
		"homologous_mapping_event;multi_event_alt_allele_in_normal",
		"multi_event_alt_allele_in_normal")
	filter<-rowRanges(vcf)$FILTER%in%keep&af[,"NORMAL"]<as.double(args[4])&tn_af>as.double(args[5])&n_cov>as.double(args[6])&t_cov>as.double(args[7])
	vcf.f<-vcf[filter]
	if(nrow(vcf.f)>0) {
		# save
		writeVcf(obj = vcf.f, filename = paste0(args[2],args[3],".mut2.filt.vcf"))
	}
} 

quit()
