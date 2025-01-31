################################################
## 1. Estimate and describe genomic diversity
################################################
# 
rule plot_Figure1_GenomicDiversity:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		PiPerRegFiles = expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, ObsExp="Observed_Data", BED=["Exons", "Introns", "Promoters", "Intergenic"], GroupSamples=["AtlSamples", "MedSamples"]),
		PiTotalFiles = expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, ObsExp="Observed_Data", BED=["Exons", "Introns", "Promoters", "Intergenic", "Callable"], GroupSamples=["AtlSamples", "MedSamples", "AllSamples"])
	output:
		PDF = "results/Plotting_DNA/plot_Figure1_GenomicDiversity.pdf"
	log:
		err = "logs/Plotting_DNA/plot_Figure1_GenomicDiversity.err",
		out = "logs/Plotting_DNA/plot_Figure1_GenomicDiversity.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_Figure1_GenomicDiversity.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pFig1",
		threads = 1,
		mem = 2000
	shell:
		"./scripts/Plotting_DNA/plot_Figure1_GenomicDiversity.R \"{input.PiPerRegFiles}\" \"{input.PiPerRegFiles}\" {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"
#


