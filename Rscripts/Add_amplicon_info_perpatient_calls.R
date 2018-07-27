library(VariantAnnotation)
args <- commandArgs()

param <- ScanVcfParam()

# import
vcf <- readVcf(args[2], "hg19", param)
if(nrow(vcf)>0)	{
	# generate AMPLI INFO
	ampli<-strsplit(args[2],"\\.")[[1]][2]
	AMPLI<-rep(ampli,nrow(vcf))		
	# add to vcf (header and info)
	names<-c(rownames(info(header(vcf))),"AMPLI")
	info(header(vcf))<-rbind(info(header(vcf)), data.frame(Number="1",Type="String",Description="Amplicon_name"))
	rownames(info(header(vcf)))<-names
	info(vcf)$AMPLI<-AMPLI
	writeVcf(obj = vcf, filename = paste0(args[2],".ampli"),index=F)
}

quit()