# AGenDiv
### Source code for the work intended to study in depth the genomic diversity of the European amohioxus

This code is structured on a snakemake pipeline. It was generated and executed on the UNIL CURNAGL computing cluster. It is heavily dependent on the CURNAGL folder structure, data storage system and its Slurm workload manager.  

## Basic environment
- conda 25.11.0
- snakemake 9.13.7

---
## Before starting: 
`conda activate snakemake`

---
## Snakemake usage:

### for a first dry run
`snakemake --cores 1 --software-deployment-method conda -p -n`

### for a cluster run
`snakemake -p --software-deployment-method conda -j 30 --executor cluster-generic --cluster-generic-submit-cmd "sbatch -J {params.name} -N 1 -o ./logs/.slurm/%x.out -e ./logs/.slurm/%x.err --cpus-per-task={params.threads} --mem={params.mem} -t {params.time}"`

---
## Folder organization
- **config**  *--> Configuration files, to be changed if changing genome or dataset*
- **rules**  *--> Rules that call scripts to generate results*
- **scripts**  *--> Scripts to analyze raw data and generate results*
- **envs**  *--> Software environements*	
- **metadata**  *--> Metadata of the samples*
- **results**  *--> Intermediate results & plots*
- **data**  *--> Raw data (fastq files, genome files...)* 
- **logs**
- **benchmarks**  

---
# Workflow:
## Genomic data analysis
### 1. VariantCalling_DNA.smk
This part of the pipeline goes from the sample fastq files to the GATK recalibrated and hard filtered small variant calls in VCF (one per chr). It performs 2 rounds of base recalibration. 
- The GATK hard filtering applied (both for base recalibration and for final calls):
	- We filtered out SNPs that meet at least one of these conditions:
		- QD < 2.0
		- QUAL < 30.0
		- SOR > 3.0
		- FS > 60.0
		- MQ < 40.0
		- MQRankSum < -12.5
		- ReadPosRankSum < -8.0
	- We filtered out INDELs that meet at least one of these conditions:
		- QD < 2.0
		- QUAL < 30.0
		- FS > 200.0
		- ReadPosRankSum < -20.0

*Rules & specific steps:*
- rule copy_fastq_NAS_data:
	> Copying from NAS storage fastq folder to the working data folder
- rule r1_index_genome_for_BWA:
	> Indexing genome in fasta format for mapping with BWA (bwa, samtools, picard)
- rule r2_map_dna_to_reference_BWA_PerSample:
	> Map to reference, sort & mark duplicates per sample (BWA MEM & samtools sort & gatk MarkDuplicatesSpark)
- rule r2_split_bam_PerChr_PerSample:
	> Split sample bams per chr with samtools (from now on almost everything is done per chr)
- rule r2_merge_3_lanes_bams_markdup_PerChr
	> Sort by read name the three lane bams of each sample & merge them while marking duplicates per chr (samtools sort -n & gatk MarkDuplicatesSpark). Merges sequentially all samples. This was done like this instead of parellelizung by chr and sample to get rid of the sample&lane identifier and use only sample identifiers from now on (except for the recalifications that are done on the sample lane bams).
	>       WARNING! This rule cannot be executed on more than 15 parallel jobs at a time.
- rule r3_base_recalibration_GATK_PerChr_round1:
	> First round of GATK BaseRecalibrator using as known sites the hard filtered, uncalibrated calls
- rule r3_base_recalibration_GATK_PerChr_round2:
	> Second round of GATK BaseRecalibrator using as known sites the hard filtered, once-calibrated calls (calibrated in the first round)
- rule r3_base_recalibration_GATK_PerChr_Unification:
	> Base Recalibration Bootstrapping unification & plotting (check convergence of covariates)
- rule r4_call_short_variants_GATK_PerChr_PerSample:
	> Calling short variants for a given chr & sample from sorted & markduplicated bam (GATK HaplotypeCaller)
- rule r5_join_sample_calls_PerChr:
	> Joining sample short variants calls with GATK GenomicsDBImport & GenotypeGVCFs per each chr
- rule r6_hardfilter_VCF_PerChr
	> Hard filter VCF per each chr
- rule r7_SNPable_Regions
	> Find SNPable regions of the reference genome. 
	> Following SNPable pipeline http://lh3lh3.users.sourceforge.net/snpable.shtml. 
	> Source code downloaded on 23/04/11 from https://lh3lh3.users.sourceforge.net/download/seqbility-20091110.tar.bz2
- rule r7_coverage_per_site_PerChr
	> Calculate coverage per chr for all samples
- rule r7_callable_regions
	> Determine callable regions. This is, regions that are SNPable, have at least a minimum coverage of mincoveragethreshold and have a minimum length of mincallablewsize (see config file for specific values)
- rule r8_filter_VCF_for_callable_regions_PerChr
	> Apply callable regions filter to the already hard filtered VCF per each chr

### 2. VariantAnalysis_DNA.smk
	- rule Prepare_basic_analysis_files_PerChr:
		> Prepare files for popgen analysis
	- rule Call_basic_stats_PopGen_perWindow_PerChr:
		> Calculate basic popgen stats x chr & x window
	- rule Unifying_Basic_stats_PopGen:
		> Unify all chr stats, calculate general basic popgen stats
	- rule Prepare_GenomicFeatures_BED:
		> Prepare a bed file with the list of genomic features
	- rule Divide_variants_perGenomicFeature_PerChr:
		> Classify each variant according to its Genomic Feature

### 3. Plotting_DNA.smk
	- rule plot_Figure1_GenomicDiversity
	- rule plot_Figure2_PopulationStructure
	- rule plot_Figure3_vsSimulations
	







