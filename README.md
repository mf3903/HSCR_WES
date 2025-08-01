**Codes and scripts for key steps in HSCR WES analysis**

**Sample selection**
-	case: selectVariant_bysample.sh <br>
-	control: FlashPCA.sh & UKB_ctrl_selection_PCdistance.ipynb (control selection based on PCA) <br>

**Bioinformatics procedures**
-	WES alignment with Burrows-Wheeler Aligner (BWA-MEM): align_bwa_mem_broad.sh <br>
-	Variant calling with GATK-HaplotypeCaller <br>
  o	gVCF.sh <br>
  o	different ploidy for males and females on sex chromosomes when generating gVCF: gVCF_sexchr_male_chrX_nonPAR.sh <br>
-	Variant QC <br>
o	VQSR: VQSR.sh <br>
o	Hard filter: variantannot_hard_filter.sh <br>
o	10% missing genotype filter: missingness_filter.sh <br>
o	Contamination <br>
  &nbsp;&nbsp;&nbsp;•	verifyBamID: verifybam.sh <br>
  &nbsp;&nbsp;&nbsp;•	Haplocheck: haplocheck.sh <br>
o	GATK’s CollectVariantCallingMetrics: vcfQC.sh <br>
o	GAKT’s GenotypeConcordance: GenotypeConcordance.sh <br>
-	Family-based de novo Mutation Calling: per_fam_DNM.sh <br>
-	Sample QC for Sex. <br>
o	read depth on the sex differential SRY gene region: sex_check_SRY.sh <br>
o	reads ratio on chrY versus chrX, normalized by the total reads for each sample: sex_check_ratio.sh <br>
-	Sample QC for Relatedness: relatedness.R & relatedness.R.sh <br>
-	Ancestry and Admixture Analysis: plink_admixture.sh <br>
-	Variant functional annotation <br>
o	Ensembl Variant Effect Predictor (VEP): vep.sh <br>
o	Pathogenicity prediction <br>
  &nbsp;&nbsp;&nbsp;•	phylop241way: phylop241.sh <br>
  &nbsp;&nbsp;&nbsp;•	VEST4: https://www.cravat.us/CRAVAT/help.jsp <br>
  &nbsp;&nbsp;&nbsp;•	LOFTEE : vep_LOFTEE.sh & vep_annot_parse.ipynb <br>
  &nbsp;&nbsp;&nbsp;•	spliceAI: https://github.com/Illumina/SpliceAI <br>
  &nbsp;&nbsp;&nbsp;•	REVEL: https://sites.google.com/site/revelgenomics/downloads?authuser=0 <br>
  &nbsp;&nbsp;&nbsp;•	metaRNN: metaRNN.sh <br>
&nbsp;&nbsp;&nbsp;•	gnomadAD AF: gnomadv4_AF_by_population.sh <br>
-	Human Embryonic Gut Gene Expression: Seurat: scRNAseq_gut.R <br>

**Statistical analysis**
-	Firth logistic regression: logit_firth.R & logit_firth.sh <br>
-	Bootstrapping: bootstrapping.R & bootstrapping.R.sh <br>
-	SKAT-O: SKAT.R & SKAT.R.sh <br>
-	Inflation factor: qqperm.R & qqperm.R.sh <br>
-	de novo Analysis: DNM.R <br>
-	Joint-analysis with extTADA: https://github.com/hoangtn/extTADA <br>
