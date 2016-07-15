#######
# WES #
#######

#Base URL of the current project
base_path="/mnt/Data4/working_directory/garan/validation_nextseq"

#Java executable path
#Usually this set by a symbolic link in /usr/bin to the most recent version available of Java
#However previously this has not always been up-to-date enough for Picard
#Explicitly state the location of the Java executable to use
java_path="/usr/share/java-1.8.0/java"

#GATK URL executable
gatk_version="3.6.0"
gatk_path="/usr/share/gatk/GenomeAnalysisTK-3.6.0/GenomeAnalysisTK.jar"


rtg_path="/mnt/Data4/working_directory/garan/rtg-tools-3.6.2/rtg"

#Exome filtering files
nkmi_path="/mnt/Data1/resources/common_no_known_medical_impact_20140929.vcf"
artefacts_path="/mnt/Data2/exome_sequencing/HaplotypeCallerVariants/WE107-WE298.repeat.HC_3.1-1.variants_to_exclude.vcf"
exac_path="/mnt/Data1/resources/ExAC.r0.1.sites.vep.AF5.vcf"

#Reference URLs
#Reference genome path, currently 1000genomes GRCh37
#http://www.1000genomes.org/category/assembly
reference_path="/mnt/Data1/resources/human_g1k_v37.fasta"

sdf_path="/mnt/Data1/resources/1000g_v37_phase2.sdf"

tmp_path="/mnt/Data3/tmp"

function SelectWesInterval(){
	#Select variants from GIAB reference by SureSelect_All_Exon_V5_S0430110
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/resources/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz \
	-o ${base_path}/vcfs/selected/NA12878_GIAB_wes_intervals_all_exon.vcf \
	-L ${base_path}/intervals/SureSelect_All_Exon_V5_S04380110_Regions.interval_list
	
	#Select variants from WE317 by SureSelect_All_Exon_V5_S0430110
	#WE317 should already have been processed with these intervals but this to standardise the inputs
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/WE0317.hiseq.vcf \
	-o ${base_path}/vcfs/selected/WE0317_hiseq_wes_intervals_all_exon.vcf \
	-L ${base_path}/intervals/SureSelect_All_Exon_V5_S04380110_Regions.interval_list
	
	#Select variants from WE417 by SureSelect_All_Exon_V5_S0430110
	#WE417 should already have been processed with these intervals but this to standardise the inputs
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/WE0417.nextseq.vcf \
	-o ${base_path}/vcfs/selected/WE0417_nextseq_wes_intervals_all_exon.vcf \
	-L ${base_path}/intervals/SureSelect_All_Exon_V5_S04380110_Regions.interval_list
	
	#Select variants from WE317 extra reads by SureSelect_All_Exon_V5_S0430110
	#WE317 should already have been processed with these intervals but this to standardise the inputs
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/WE0317.nextseq.extra.vcf \
	-o ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_intervals_all_exon.vcf \
	-L ${base_path}/intervals/SureSelect_All_Exon_V5_S04380110_Regions.interval_list
	
}

function SelectHighConfInterval(){
	#Select variants from GIAB reference by SureSelect_All_Exon_V5_S0430110
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/selected/NA12878_GIAB_wes_intervals_all_exon.vcf \
	-o ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf \
	-L ${base_path}/intervals/GIAB_highconf.bed
	
	#Select variants from WE317 by SureSelect_All_Exon_V5_S0430110
	#WE317 should already have been processed with these intervals but this to standardise the inputs
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/selected/WE0317_hiseq_wes_intervals_all_exon.vcf \
	-o ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf \
	-L ${base_path}/intervals/GIAB_highconf.bed
	
	#Select variants from WE417 by SureSelect_All_Exon_V5_S0430110
	#WE417 should already have been processed with these intervals but this to standardise the inputs
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/selected/WE0417_nextseq_wes_intervals_all_exon.vcf \
	-o ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf \
	-L ${base_path}/intervals/GIAB_highconf.bed
	
	#Select variants from WE317 extra by SureSelect_All_Exon_V5_S0430110
	#WE317 should already have been processed with these intervals but this to standardise the inputs
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T SelectVariants \
	-R ${reference_path} \
	--variant ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_intervals_all_exon.vcf \
	-o ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.vcf \
	-L ${base_path}/intervals/GIAB_highconf.bed
	
}

function GenerateStats(){
	
		${rtg_path} vcfstats ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf > ${base_path}/results/NA12878_GIAB.stats
		
		${rtg_path} vcfstats ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf > ${base_path}/results/WE0317.stats
		
		${rtg_path} vcfstats ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf > ${base_path}/results/WE0417.stats
		
		${rtg_path} vcfstats ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.vcf > ${base_path}/results/WE0317.stats
}

function BgzipVcfs() {
	bgzip ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf
	tabix -p vcf ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf.gz
	
	bgzip ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf
	tabix -p vcf ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf.gz
	
	bgzip ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf
	tabix -p vcf ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf.gz
	
	bgzip ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.vcf
	tabix -p vcf ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.vcf.gz
}

function CompareGiab(){
	
	${rtg_path} vcfeval -t ${sdf_path} -b ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf.gz -c ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf.gz -o ${base_path}/results/GIAB_vs_WE0317-hiseq
	
	${rtg_path} vcfeval -t ${sdf_path} -b ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf.gz -c ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf.gz -o ${base_path}/results/GIAB_vs_WE0417-nextseq

	${rtg_path} vcfeval -t ${sdf_path} -b ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf.gz -c ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.vcf.gz -o ${base_path}/results/GIAB_vs_WE0317-nextseq-extra
}

function CompareNextseqVsHiseq(){
	
	${rtg_path} vcfeval -t ${sdf_path} -b ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf.gz -c ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf.gz -o ${base_path}/results/WE0317-hiseq_vs_WE0417-nextseq

}

function AddFilters(){
	#Add Filter annotation
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T VariantFiltration \
	-R ${reference_path} \
	--filterExpression "QD < 2.0" --filterName "QD2" \
	--filterExpression "MQ < 40.0" --filterName "MQ40" \
	--filterExpression "ReadPosRankSum < -8.0" --filterName "RPRS-8" \
	--filterExpression "FS > 60.0" --filterName "FS60" \
	--filterExpression "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
	--variant ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.vcf.gz \
	-o ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.filtered.vcf.gz
	
	
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T VariantFiltration \
	-R ${reference_path} \
	--filterExpression "QD < 2.0" --filterName "QD2" \
	--filterExpression "MQ < 40.0" --filterName "MQ40" \
	--filterExpression "ReadPosRankSum < -8.0" --filterName "RPRS-8" \
	--filterExpression "FS > 60.0" --filterName "FS60" \
	--filterExpression "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
	--variant ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.vcf.gz \
	-o ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.filtered.vcf.gz
	
	
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T VariantFiltration \
	-R ${reference_path} \
	--filterExpression "QD < 2.0" --filterName "QD2" \
	--filterExpression "MQ < 40.0" --filterName "MQ40" \
	--filterExpression "ReadPosRankSum < -8.0" --filterName "RPRS-8" \
	--filterExpression "FS > 60.0" --filterName "FS60" \
	--filterExpression "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
	--variant ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.vcf.gz \
	-o ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.filtered.vcf.gz
	
	${java_path} -Xmx12g -jar ${gatk_path} \
	-T VariantFiltration \
	-R ${reference_path} \
	--filterExpression "QD < 2.0" --filterName "QD2" \
	--filterExpression "MQ < 40.0" --filterName "MQ40" \
	--filterExpression "ReadPosRankSum < -8.0" --filterName "RPRS-8" \
	--filterExpression "FS > 60.0" --filterName "FS60" \
	--filterExpression "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
	--variant ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.vcf.gz \
	-o ${base_path}/vcfs/selected/WE0317_nextseq_extra_wes_high-conf_intervals.filtered.vcf.gz
}

function SelectByExomeFilters(){
	
	#NKMI
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.filtered.vcf.gz --discordance ${nkmi_path} -o ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.nkmi.filtered.vcf.gz
	
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.filtered.vcf.gz --discordance ${nkmi_path} -o ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.nkmi.filtered.vcf.gz
	
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.filtered.vcf.gz --discordance ${nkmi_path} -o ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.nkmi.filtered.vcf.gz
	
	#ARTEFACTS
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.nkmi.filtered.vcf.gz --discordance ${artefacts_path} -o ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.artefacts.nkmi.filtered.vcf.gz
	
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.nkmi.filtered.vcf.gz --discordance ${artefacts_path} -o ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.artefacts.nkmi.filtered.vcf.gz
	
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.nkmi.filtered.vcf.gz --discordance ${artefacts_path} -o ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.artefacts.nkmi.filtered.vcf.gz
	
	#EXAC
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.artefacts.nkmi.filtered.vcf.gz --discordance ${exac_path} -o ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz
	
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.artefacts.nkmi.filtered.vcf.gz --discordance ${exac_path} -o ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz
	
	${java_path} -Djava.io.tmpdir=${tmp_path}-Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -jar ${gatk_path} -T SelectVariants -R ${reference_path} -V ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.artefacts.nkmi.filtered.vcf.gz --discordance ${exac_path} -o ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz	
}

function GenerateFilteredStats(){
	
		${rtg_path} vcfstats ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz > ${base_path}/results/NA12878_GIAB.filtered.stats
		
		${rtg_path} vcfstats ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz > ${base_path}/results/WE0317.filtered.stats
		
		${rtg_path} vcfstats ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz > ${base_path}/results/WE0417.filtered.stats
}

function CompareFilteredGiab(){
	
	${rtg_path} vcfeval -t ${sdf_path} -b ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz -c ${base_path}/vcfs/selected/WE0317_hiseq_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz -o ${base_path}/results/GIAB_vs_WE0317-hiseq_filtered
	
	${rtg_path} vcfeval -t ${sdf_path} -b ${base_path}/vcfs/selected/NA12878_GIAB_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz -c ${base_path}/vcfs/selected/WE0417_nextseq_wes_high-conf_intervals.exac.artefacts.nkmi.filtered.vcf.gz -o ${base_path}/results/GIAB_vs_WE0417-nextseq_filtered

}

function DepthOfCoverage(){	
		 ${java_path} -Xmx8g -jar ${gatk_path} -T DepthOfCoverage \
		-L ${base_path}/intervals/fn.bed \
		-R ${reference_path} \
		-I /mnt/Data5/exome_sequencing/nextseq/WE0417/r01_assembly/WE0417.bam \
		-o ${base_path}/results/coverage/WE0417.by_base_coverage \
		-mmq 30 -mbq 30 -dels -ct 1 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -omitLocusTable
}


#Initial analysis of unfiltered variants from GIAB, Hiseq and Nextseq
SelectWesInterval && SelectHighConfInterval && GenerateStats && BgzipVcfs && CompareGiab 
#&& CompareNextseqVsHiseq

#Second part of analysis. Variants filtered by hard-filters used in current exome pipeline (17/06/2016) *note not GATK best practice!
#This will effect the metric weighting during ROC analysis (no VQSLOD)
#AddFilters && SelectByExomeFilters && GenerateFilteredStats && CompareFilteredGiab
#DepthOfCoverage
