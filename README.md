Codes and scripts for key steps in HSCR WES analysis

Sample selection
-	case: selectVariant_bysample.sh
-	control: control selection based on PCA: FlashPCA.sh & UKB_ctrl_selection_PCdistance.ipynb

Bioinformatics procedures
-	WES alignment with Burrows-Wheeler Aligner (BWA-MEM): align_bwa_mem_broad.sh
-	Variant calling with GATK-HaplotypeCaller
o	gVCF.sh
o	different ploidy for males and females on sex chromosomes when generating gVCF: gVCF_sexchr_male_chrX_nonPAR.sh
-	Variant QC
o	VQSR: VQSR.sh
o	Hard filter: variantannot_hard_filter.sh
o	10% missing genotype filter: missingness_filter.sh
o	Contamination
•	verifyBamID: verifybam.sh
•	Haplocheck: haplocheck.sh
o	GATK’s CollectVariantCallingMetrics: vcfQC.sh
o	GAKT’s GenotypeConcordance: GenotypeConcordance.sh
-	Family-based de novo Mutation Calling: per_fam_DNM.sh
-	Sample QC for Sex. 
o	read depth on the sex differential SRY gene region: sex_check_SRY.sh
o	reads ratio on chrY versus chrX, normalized by the total reads for each sample: sex_check_ratio.sh
-	Sample QC for Relatedness: relatedness.R & relatedness.R.sh
-	Ancestry and Admixture Analysis: plink_admixture.sh
-	Variant functional annotation 
o	Ensembl Variant Effect Predictor (VEP): vep.sh
o	Pathogenicity prediction
•	phylop241way: phylop241.sh
•	VEST4: https://www.cravat.us/CRAVAT/help.jsp 
•	LOFTEE : vep_LOFTEE.sh & vep_annot_parse.ipynb
•	spliceAI: https://github.com/Illumina/SpliceAI
•	REVEL: https://sites.google.com/site/revelgenomics/downloads?authuser=0
•	metaRNN: 
•	gnomadAD AF: gnomadv4_AF_by_population.sh
-	Human Embryonic Gut Gene Expression: Seurat: scRNAseq_gut.R

Statistical analysis. 
-	Firth logistic regression: logit_firth.R & logit_firth.sh
-	Bootstrapping: bootstrapping.R & bootstrapping.R.sh
-	SKAT-O: SKAT.R & SKAT.R.sh
-	Inflation factor: qqperm.R & qqperm.R.sh
-	de novo Analysis: DNM.R
-	Joint-analysis with extTADA: https://github.com/hoangtn/extTADA
<img width="468" height="564" alt="image" src="https://github.com/user-attachments/assets/b7175e9e-18f1-441b-80f1-ed536ef755c1" />

