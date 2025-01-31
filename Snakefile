

configfile: "config/AGenDiv.yaml"
localrules: copy_fastq_NAS_data, DownloadProteinSeq_from_Ensembl, BuildOGs_withOrthoLoger

include: 'rules/VariantCalling_DNA.smk'
include: 'rules/VariantAnalysis_DNA.smk'
include: 'rules/Plotting_DNA.smk'
include: 'rules/Simulations.smk'
include: 'rules/SpeciesTreeBuilding.smk'


rule all:
	'''
	List of final output files.
	'''
	input:
		"results/VariantCalling_DNA/9_run_all_chr/ShortVariants_HardExtraCallableFiltered.vcf.gz",
		#"results/Plotting_DNA/plot_SharedPrivate_Populations/plot_SharedPrivate_Populations.pdf",
		#"results/Plotting_DNA/plot_PerFeatureStats/plot_PerFeatureStats.pdf", 
		#"results/Plotting_DNA/plot_PerSampleStats/plot_PerSampleStats.pdf", 
		#"results/Plotting_DNA/plot_PerWindowStats/plot_PerWindowStats.pdf", 
		#"results/Plotting_DNA/plot_join_PSMC/joined_PSMC.pdf", 
		#"results/Plotting_DNA/CircosPlots/plot_Circos_Callable_Varsites_Pi.png", 
		#"results/Plotting_DNA/plot_CummulativePropOfDiversityExplained/plot_CummulativePropOfDiversityExplained.pdf", 
		#"results/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations/plot_PerSite_VarTypes_FuncRegions_Populations.pdf",
		#"results/Simulations/plot_Pi_RangeNe_RangeMu_singlePopulation/plot_Pi_RangeNe_RangeMu_singlePopulation.pdf"




