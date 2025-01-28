
################################################
## Filtering analysis
################################################

rule VariantFilterStats_PerChr:
	'''
	Per chr variant-to-table to extract statistics on variant filters
	'''
	input:
		VCF = rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered
	output:
		FilterStats = "results/VariantAnalysis_DNA/VariantFilterStats_PerChr/ShortVariant_FilterStats.{chr}.tab.gz"
	log:
		err = "logs/VariantAnalysis_DNA/VariantFilterStats_PerChr/VariantFilterStats_PerChr.{chr}.err",
		out = "logs/VariantAnalysis_DNA/VariantFilterStats_PerChr/VariantFilterStats_PerChr.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/VariantFilterStats_PerChr/VariantFilterStats_PerChr.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '03:00:00',
		name = "FiltSt{chr}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantAnalysis_DNA/VariantFilterStats_PerChr.sh {input.VCF} {output.FilterStats} > {log.out} 2> {log.err}"

rule ProcessingCoverageData_PerChr:
	'''
	From raw coverage values to statistics per window
	'''
	input:
		Coverage = rules.r7_join_coverage_per_site.output.covgenome,
		ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
	output:
		CoverageStats = "results/VariantAnalysis_DNA/ProcessingCoverageData/CoverageWindowStats_{chr}.bed.gz"
	log:
		err = "logs/VariantAnalysis_DNA/ProcessingCoverageData/ProcessingCoverageData_{chr}.err",
		out = "logs/VariantAnalysis_DNA/ProcessingCoverageData/ProcessingCoverageData_{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/ProcessingCoverageData/ProcessingCoverageData_{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '03:00:00',
		name = "pCov{chr}",
		threads = 1,
		mem = 20000,
		slidingwindowsize=config["slidingwindowsize"], # window size for sliding window calculations 
		slidingwindowstep=config["slidingwindowstep"] # step size for sliding window calculations
	shell:
		"./scripts/VariantAnalysis_DNA/ProcessingCoverageData.sh {input.Coverage} {input.ChrLengths} {output.CoverageStats} {params.slidingwindowsize} {params.slidingwindowstep} {wildcards.chr} > {log.out} 2> {log.err}"

# To be written
rule Coocurrence_of_BED_Regions_with_Several_features:
	'''
	Statistically study the coocurrence of bed formated regions (eg. callable regions) with other features:
	- Coverage
	- GC content
	- Repeats
	- Functional features (exons, introns, promoters...)
	'''
	input:
		#Regions = ,
		#SNPableRegions = rules.r7_callable_regions.output.SNPableRegions,
		#NonCallableRegions = rules.r7_callable_regions.output.NonCallableRegions,
		#NonExtraCallableRegions = rules.r7_extra_callable_regions.output.NonExtraCallableRegions
	output:
		CoocurrenceOutput = "results/VariantAnalysis_DNA/Coocurrence_of_BED_Regions_with_Several_features/Coocurrence_of_BED_Regions_with_Several_features.{type}.txt"
	log:
		err = "logs/VariantAnalysis_DNA/Coocurrence_of_BED_Regions_with_Several_features/Coocurrence_of_BED_Regions_with_Several_features.{type}.err",
		out = "logs/VariantAnalysis_DNA/Coocurrence_of_BED_Regions_with_Several_features/Coocurrence_of_BED_Regions_with_Several_features.{type}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Coocurrence_of_BED_Regions_with_Several_features/Coocurrence_of_BED_Regions_with_Several_features.{type}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '03:00:00',
		name = "Coo{type}",
		threads = 1,
		mem = 5000,
		lowermincovthres = config["lowermincovthres"], # for a given position, the minimum coverage among all samples has to be at least lowermincovthres
		uppermincovthres = config["uppermincovthres"]  # for a given position, the minimum coverage among all samples has to be less or equal to uppermincovthres
	shell:
		"./scripts/VariantAnalysis_DNA/Coocurrence_of_BED_Regions_with_Several_features.sh {input.ChrsFilterStats} {input.extraVCFs} {output.statsTXT} > {log.out} 2> {log.err}"

################################################
## Variant analysis
################################################

rule mask_callable_regions_from_genome:
	'''
	Mask all non-callable regions of the genome
	'''
	input:
		NonExtraCallableRegions = rules.r7_join_extra_callable_regions.output.NonExtraCallableRegions,
		genome = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		MaskedGenome = "results/VariantAnalysis_DNA/mask_callable_regions_from_genome/Branchiostoma_lanceolatum.BraLan3_genome_CallableRegionsMasked.fa"
	log:
		err = "logs/VariantAnalysis_DNA/mask_callable_regions_from_genome/mask_callable_regions_from_genome.err",
		out = "logs/VariantAnalysis_DNA/mask_callable_regions_from_genome/mask_callable_regions_from_genome.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/mask_callable_regions_from_genome/mask_callable_regions_from_genome.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '03:00:00',
		name = "MaskGenome",
		threads = 1,
		mem = 5000
	shell:
		'''
		out={output.MaskedGenome}
		rm -r $(dirname ${{out}})
		mkdir -p $(dirname ${{out}})
		bedtools maskfasta -fi {input.genome} -bed {input.NonExtraCallableRegions} -fo ${{out%fa}}unnorm.fa > {log.out} 2> {log.err}
		picard NormalizeFasta -I ${{out%fa}}unnorm.fa -O ${{out}} >> {log.out} 2>> {log.err}
		'''

rule Prepare_basic_analysis_files_PerChr:
	'''
	Prepare files for further analysis:
		- .genotype.bed file
		- .genotype.numericmatrix file (#alleles-1 lines for each site)
		- .0hom1het.bed file
	'''
	input:
		chrVCF = rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered,
		SamplesOrderInVCF = rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF
	output:
		chrGENOTYPE = "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.genotype.bed.gz",
		chrGENOTYPEMatrix = "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.genotype.numericmatrix.gz",
		HeterozygosityFile = "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.0hom1het.bed.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/Prepare_basic_analysis_files_PerChr.{chr}.err",
		out = "logs/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/Prepare_basic_analysis_files_PerChr.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/Prepare_basic_analysis_files_PerChr.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '5:00:00',
		name = "{chr}Genotype",
		threads = 1,
		mem = 25000
	shell:
		"./scripts/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr.sh {input.chrVCF} {output.chrGENOTYPE} > {log.out} 2> {log.err}"

rule Prepare_GenomicFeatures_BED:
	'''
	Prepare a bed file with the list of genomic features
	'''
	input:
		ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt",
		GTF = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_strong.gtf",
		ExtraCallableRegions = rules.r7_join_extra_callable_regions.output.ExtraCallableRegions
	output:
		GenomicFeatures = "results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures.bed",
		GenomicFeaturesExtraCallable = "results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures_ExtraCallableRegions.bed"
	log:
		err = "logs/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Prepare_GenomicFeatures_BED.err",
		out = "logs/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Prepare_GenomicFeatures_BED.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Prepare_GenomicFeatures_BED.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '5:00:00',
		name = "GFbed",
		threads = 1,
		mem = 10000
	shell:
		"./scripts/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED.sh {input.ChrLengths} {input.GTF} {input.ExtraCallableRegions} {output.GenomicFeatures} {output.GenomicFeaturesExtraCallable} > {log.out} 2> {log.err}"

rule Divide_variants_perGenomicFeature_PerChr:
	'''
	Divide variants per Genomic Feature, being biallelic or not, being polymorphisms or not and its major allele grequency
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		chrGENOTYPE = rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE,
		GenomicFeaturesExtraCallable = rules.Prepare_GenomicFeatures_BED.output.GenomicFeaturesExtraCallable
	output:
		chrDividedGENOTYPEs = "results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.{chr}_classified.genotype.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/Divide_variants_perGenomicFeature_PerChr.{chr}.err",
		out = "logs/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/Divide_variants_perGenomicFeature_PerChr.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/Divide_variants_perGenomicFeature_PerChr.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '10:00:00',
		name = "{chr}Div",
		threads = 1,
		mem = 10000,		
		maxfreqploymorphism=config["maxfreqploymorphism"], # maximum frequency of the majoritary allele to be considered a polymorphism (0.97 --> 3/72 is considered polym but 2/72 no // 0.98 --> 2/72 polym 1/72 no)
		samples=expand("{sample}", sample=config["samples"])
	shell:
		"./scripts/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr.sh {input.Metadata} {input.chrGENOTYPE} {input.GenomicFeaturesExtraCallable} {output.chrDividedGENOTYPEs} {params.maxfreqploymorphism} {params.samples} > {log.out} 2> {log.err}"

rule Sample_clustering:
	'''
	Whole genome sample clustering prepared to be called for all possible combinations of:
		- PCA type of clustering {wildcards.typeClust}
	    - all / exons / introns / promoters / intergenic / genic variants {wildcards.feat}
	    - all / polymorphic variants {wildcards.polym}
	    - all / biallelic variants {wildcards.numall}
	    - all / specific major frequency: {wildcards.freqint}
	    	- intervals e.g. 0.5to0.6 will include all variants with major freq >= 0.5 and < 0.6
	    	- 71  
	'''
	input:
		chrsGENOTYPEMatrices = expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPEMatrix, chr=config["bralan3chrs"]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt", 
		DividedSites = expand(rules.Divide_variants_perGenomicFeature_PerChr.output.chrDividedGENOTYPEs, chr=config["bralan3chrs"])
	output:
		PCsites = "results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCvalues_PerVariant_{typeClust}_{feat}_{polym}_{numall}_{freqint}.txt.gz",
		PCsamples = "results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCvalues_PerSample_{typeClust}_{feat}_{polym}_{numall}_{freqint}.txt.gz",
		PCprop = "results/VariantAnalysis_DNA/Sample_clustering/Numericmatrix_PCprop_{typeClust}_{feat}_{polym}_{numall}_{freqint}.txt.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Sample_clustering/Sample_clustering{typeClust}_{feat}_{polym}_{numall}_{freqint}.err",
		out = "logs/VariantAnalysis_DNA/Sample_clustering/Sample_clustering{typeClust}_{feat}_{polym}_{numall}_{freqint}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Sample_clustering/Sample_clustering{typeClust}_{feat}_{polym}_{numall}_{freqint}.txt"
	conda:
		'../envs/SampleClustering.yaml'
	params:
		time = '10:00:00',
		name = "SClust",
		threads = 1,
		mem = 500000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/VariantAnalysis_DNA/Sample_clustering.R \"{input.chrsGENOTYPEMatrices}\" {input.SamplesOrderInVCF} \"{input.DividedSites}\" {output.PCsites} {output.PCsamples} {output.PCprop} \"{params.chrs}\" {wildcards.typeClust} {wildcards.feat} {wildcards.polym} {wildcards.numall} {wildcards.freqint} > {log.out} 2> {log.err}"

# Check the translations were done propperly
# Check NSD matrix
rule Synonyms_and_nonsynonymous_fromGTF:
	'''
	Extract which sites are synonymous and which are nonsynonymous from GTF 
	'''
	input:
		GTF = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_strong.gtf",
		ExtraCallableRegions = rules.r7_join_extra_callable_regions.output.NonExtraCallableRegions,
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		SynNonSynBED = "results/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/SND_positions.bed",
		SynNonSynBED_callable = "results/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/SND_positions_callable.bed"
	log:
		err = "logs/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/Synonyms_and_nonsynonymous_fromGTF.err",
		out = "logs/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/Synonyms_and_nonsynonymous_fromGTF.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/Synonyms_and_nonsynonymous_fromGTF.txt"
	conda:
		'../envs/VariantAnnotation.yaml'
	params:
		time = '5:00:00',
		name = "SnSbed",
		threads = 1,
		mem = 10000,
		SnSfromSeq = "scripts/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromSeq.pl",
		ProtfromSeq = "scripts/VariantAnalysis_DNA/ProteinSeq_fromTranscriptSeq.pl"
	shell:
		"./scripts/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF.sh {input.GTF} {input.ExtraCallableRegions} {input.genomeFA} {output.SynNonSynBED} {output.SynNonSynBED_callable} {params.SnSfromSeq} {params.ProtfromSeq} > {log.out} 2> {log.err}"

rule PerSite_StatsMatrix_PerChr:
	'''
	Calculate popgen stats x site, put them in a big matrix file
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		chrVCF = rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered,
		chrGENOTYPE = rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE,
		HeterozygosityFile = rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile,
		MinCoverage = rules.r7_join_coverage_per_site.output.mincov,
		SamplesOrderInVCF = rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF,
		GenomicFeatureDivision = rules.Divide_variants_perGenomicFeature_PerChr.output.chrDividedGENOTYPEs,
		SynNonSynBED_callable = rules.Synonyms_and_nonsynonymous_fromGTF.output.SynNonSynBED_callable,
		PCA_All = expand(rules.Sample_clustering.output.PCsites, typeClust="PCA", feat="all", polym="all", numall="all", freqint="all"),
		PCA_Polym = expand(rules.Sample_clustering.output.PCsites, typeClust="PCA", feat="all", polym="polymorphic", numall="all", freqint="all"),
		PCA_Bial = expand(rules.Sample_clustering.output.PCsites, typeClust="PCA", feat="all", polym="all", numall="biallelic", freqint="all"),
		PCA_Exons = expand(rules.Sample_clustering.output.PCsites, typeClust="PCA", feat="exon", polym="all", numall="all", freqint="all"),
		PCA_freqint = expand(rules.Sample_clustering.output.PCsites, typeClust="PCA", feat="all", polym="all", numall="all", freqint=("0.5to0.6", "0.6to0.7", "0.7to0.8", "0.8to0.9", "0.9to1"))
	output:
		SiteMatrix="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.{chr}.tab.gz"
	log:
		err = "logs/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSite_StatsMatrix.{chr}.err",
		out = "logs/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSite_StatsMatrix.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSite_StatsMatrix.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "{chr}Site",
		threads = 1,
		mem = 150000,
		maxfreqploymorphism=config["maxfreqploymorphism"] # maximum frequency of the majoritary allele to be considered a polymorphism (0.97 --> 3/72 is considered polym but 2/72 no // 0.98 --> 2/72 polym 1/72 no)
	shell:
		"./scripts/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr.sh {input.Metadata} {input.chrVCF} {input.chrGENOTYPE} {input.MinCoverage} {input.SamplesOrderInVCF} {input.GenomicFeatureDivision} {input.SynNonSynBED_callable} {input.PCA_All} {output.SiteMatrix} {params.maxfreqploymorphism} {wildcards.chr} > {log.out} 2> {log.err}"

rule PerSample_StatsMatrix:
	'''
	Calculate popgen stats x sample, put them in a matrix file
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		VCF = rules.r9_run_all_chr.output.unifiedVCF,
		chrshomhetFiles = expand(rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile, chr=config["bralan3chrs"]),
		ExtraCallableRegions = rules.r7_join_extra_callable_regions.output.ExtraCallableRegions,
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt",
		chrsSiteMatrices=expand(rules.PerSite_StatsMatrix_PerChr.output.SiteMatrix, chr=config["bralan3chrs"]),
		PCA_All = expand(rules.Sample_clustering.output.PCsamples, typeClust="PCA", feat="all", polym="all", numall="all", freqint="all"),
		PCA_Polym = expand(rules.Sample_clustering.output.PCsamples, typeClust="PCA", feat="all", polym="polymorphic", numall="all", freqint="all"),
		PCA_Bial = expand(rules.Sample_clustering.output.PCsamples, typeClust="PCA", feat="all", polym="all", numall="biallelic", freqint="all"),
		PCA_Exons = expand(rules.Sample_clustering.output.PCsamples, typeClust="PCA", feat="exon", polym="all", numall="all", freqint="all"),
		PCA_freqint = expand(rules.Sample_clustering.output.PCsamples, typeClust="PCA", feat="all", polym="all", numall="all", freqint=("0.5to0.6", "0.6to0.7", "0.7to0.8", "0.8to0.9", "0.9to1")),
		UMAP_All = expand(rules.Sample_clustering.output.PCsamples, typeClust="UMAP", feat="all", polym="all", numall="all", freqint="all")
	output:
		SampleMatrix="results/VariantAnalysis_DNA/PerSample_StatsMatrix/PerSampleStats.tab.gz"
	log:
		err = "logs/VariantAnalysis_DNA/PerSample_StatsMatrix/PerSample_StatsMatrix.err",
		out = "logs/VariantAnalysis_DNA/PerSample_StatsMatrix/PerSample_StatsMatrix.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/PerSample_StatsMatrix/PerSample_StatsMatrix.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "WGSample",
		threads = 1,
		mem = 25000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/VariantAnalysis_DNA/PerSample_StatsMatrix.sh {input.Metadata} {input.VCF} \"{input.chrshomhetFiles}\" {input.ExtraCallableRegions} {input.SamplesOrderInVCF} \"{input.chrsSiteMatrices}\" {input.PCA_All} {input.UMAP_All} {output.SampleMatrix} \"{params.chrs}\" > {log.out} 2> {log.err}"

rule PerSample_StatsMatrix_PerChr:
	'''
	Calculate popgen stats x sample, put them in a matrix file
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		chrVCF = rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered,
		chrGENOTYPE = rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE,
		chrGENOTYPEMatrix = rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPEMatrix,
		HeterozygosityFile = rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile,
		ExtraCallableRegions = rules.r7_join_extra_callable_regions.output.ExtraCallableRegions,
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.{chr}.txt",
		SiteMatrix = rules.PerSite_StatsMatrix_PerChr.output.SiteMatrix
	output:
		SampleMatrix="results/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSampleStats.{chr}.tab"
	log:
		err = "logs/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSample_StatsMatrix.{chr}.err",
		out = "logs/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSample_StatsMatrix.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr/PerSample_StatsMatrix.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "{chr}Sample",
		threads = 1,
		mem = 25000
	shell:
		"./scripts/VariantAnalysis_DNA/PerSample_StatsMatrix_PerChr.sh {input.Metadata} {input.chrVCF} {input.chrGENOTYPE} {input.ExtraCallableRegions} {input.SamplesOrderInVCF} {input.SiteMatrix} {output.SampleMatrix} {wildcards.chr} > {log.out} 2> {log.err}"

rule PerWindow_StatsMatrix_PerChr:
	'''
	Calculate popgen stats x window, put them in a big matrix file
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt",
		chrVCF = rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered,
		chrGENOTYPE = rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE,
		HeterozygosityFile = rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile,
		ExtraCallableRegions = rules.r7_extra_callable_regions_PerChr.output.ExtraCallableRegions,
		ValidWindows = rules.r7_extra_callable_regions_PerChr.output.ValidWindows,
		GenomicFeaturesExtraCallable = rules.Prepare_GenomicFeatures_BED.output.GenomicFeaturesExtraCallable,
		SamplesOrderInVCF="metadata/SamplesOrderInVCF.{chr}.txt",
		SiteMatrix=rules.PerSite_StatsMatrix_PerChr.output.SiteMatrix,
		MaskedGenome = rules.mask_callable_regions_from_genome.output.MaskedGenome
	output:
		WindowMatrix = expand("results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_{size}_{step}.{{chr}}.tab.gz", size=config["slidingwindowsize"], step=config["slidingwindowstep"])
	log:
		err = expand("logs/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindow_StatsMatrix_{size}_{step}.{{chr}}.err", size=config["slidingwindowsize"], step=config["slidingwindowstep"]),
		out = expand("logs/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindow_StatsMatrix_{size}_{step}.{{chr}}.out", size=config["slidingwindowsize"], step=config["slidingwindowstep"])
	benchmark:
		"benchmarks/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindow_StatsMatrix.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "{chr}Window",
		threads = 1,
		mem = 25000,
		slidingwindowsize=config["slidingwindowsize"], # window size for sliding window calculations 
		slidingwindowstep=config["slidingwindowstep"], # step size for sliding window calculations
		maxfreqploymorphism=config["maxfreqploymorphism"] # maximum frequency of the majoritary allele to be considered a polymorphism (0.97 --> 3/72 is considered polym but 2/72 no // 0.98 --> 2/72 polym 1/72 no)
	shell:
		"./scripts/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr.sh {input.Metadata} {input.ChrLengths} {input.chrVCF} {input.chrGENOTYPE} {input.ExtraCallableRegions} {input.ValidWindows} {input.GenomicFeaturesExtraCallable} {input.SamplesOrderInVCF} {input.SiteMatrix} {input.MaskedGenome} {output.WindowMatrix} {params.slidingwindowsize} {params.slidingwindowstep} {params.maxfreqploymorphism} {wildcards.chr} > {log.out} 2> {log.err}"

rule PerFeature_StatsMatrix_PerChr:
	'''
	Calculate popgen stats x genomic feature, put them in a big matrix file
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt",
		chrVCF = rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered,
		chrGENOTYPE = rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE,
		HeterozygosityFile = rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile,
		ExtraCallableRegions = rules.r7_extra_callable_regions_PerChr.output.ExtraCallableRegions,
		SNPableRegions ="results/VariantCalling_DNA/7_callable_regions/SNPableRegions.bed.gz",
		GenomicFeatures = "results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures.bed",
		GenomicFeaturesExtraCallable = "results/VariantAnalysis_DNA/Prepare_GenomicFeatures_BED/Branchiostoma_lanceolatum.BraLan3_strong_GenomicFeatures_ExtraCallableRegions.bed",
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.{chr}.txt",
		SynNonSynBED = "results/VariantAnalysis_DNA/Synonyms_and_nonsynonymous_fromGTF/SND_positions_callable.bed",
		SiteMatrix="results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.{chr}.tab.gz"
	output:
		FeatureMatrix="results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.{chr}.tab.gz"
	log:
		err = "logs/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeature_StatsMatrix.{chr}.err",
		out = "logs/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeature_StatsMatrix.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeature_StatsMatrix.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "{chr}Feature",
		threads = 1,
		mem = 25000,
		maxfreqploymorphism=config["maxfreqploymorphism"] # maximum frequency of the majoritary allele to be considered a polymorphism (0.97 --> 3/72 is considered polym but 2/72 no // 0.98 --> 2/72 polym 1/72 no)
	shell:
		"./scripts/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr.sh {input.Metadata} {input.ChrLengths} {input.chrVCF} {input.chrGENOTYPE} {input.ExtraCallableRegions} {input.SNPableRegions} {input.GenomicFeatures} {input.GenomicFeaturesExtraCallable} {input.SamplesOrderInVCF} {input.SynNonSynBED} {input.SiteMatrix} {output.FeatureMatrix} {params.maxfreqploymorphism} {wildcards.chr} > {log.out} 2> {log.err}"

rule Prepare_FSTAT_file:
	'''
	Goudet, J. (2005) Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes. 5: 184-186 [pdf]. Latest version available at github .
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		chrsGENOTYPEMatrices = expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPEMatrix, chr=config["bralan3chrs"]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt" 
	output:
		FstatDataFile="results/VariantAnalysis_DNA/Prepare_FSTAT_file/ShortVariants.filtered.genotype.fstat.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Prepare_FSTAT_file/Prepare_FSTAT_file.err",
		out = "logs/VariantAnalysis_DNA/Prepare_FSTAT_file/Prepare_FSTAT_file.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Prepare_FSTAT_file/Prepare_FSTAT_file.txt"
	conda:
		'../envs/hierfstat.yaml'
	params:
		time = '15:00:00',
		name = "prepFSTAT",
		threads = 1,
		mem = 10000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/VariantAnalysis_DNA/Prepare_FSTAT_file.sh {input.Metadata} \"{input.chrsGENOTYPEMatrices}\" {input.SamplesOrderInVCF} {output.FstatDataFile} \"{params.chrs}\" > {log.out} 2> {log.err}"

rule Hierfstat_statistics:
	'''
	Goudet, J. (2005) Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes. 5: 184-186 [pdf]. Latest version available at github .
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		chrGENOTYPEs = expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE, chr=config["bralan3chrs"]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt" 
	output:
		HierFstatInput = "results/VariantAnalysis_DNA/Hierfstat_statistics/Hierfstat_input.tbl.gz",
		OverallStats = "results/VariantAnalysis_DNA/Hierfstat_statistics/Hierfstat_OverallStats.tbl"
	log:
		err = "logs/VariantAnalysis_DNA/Hierfstat_statistics/Hierfstat_statistics.err",
		out = "logs/VariantAnalysis_DNA/Hierfstat_statistics/Hierfstat_statistics.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Hierfstat_statistics/Hierfstat_statistics.txt"
	conda:
		'../envs/hierfstat.yaml'
	params:
		time = '70:00:00',
		name = "OverallStats",
		threads = 1,
		mem = 200000
	shell:
		"""
		mkdir -p $(dirname {output.HierFstatInput}) 2> {log.err}
		zcat {input.chrGENOTYPEs} | sed 's/:/\\t/g' | cut -f5- | awk '{{num=0;for(i=1;i<=NF;i++){{if(!a[$i]){{a[$i]=num;num++}}}}; str=""; for(i=1;i<=NF;i++){{if(a[$i]<9){{str=str"\\t0"a[$i]}}else{{str=str"\\t"a[$i]}}}}print str; delete a}}' | sed 's/^\\t//g' | awk '{{str=$1$2; for(i=3;i<NF;i+=2){{j=i+1; str=str"\\t"$i$j}}print str}}' | gzip > {output.HierFstatInput} 2>> {log.err}
		./scripts/VariantAnalysis_DNA/Hierfstat_statistics.R {input.Metadata} {input.SamplesOrderInVCF} {output.HierFstatInput} {output.OverallStats} > {log.out} 2>> {log.err}
		"""

rule Build_ConsensusSequences_PerSample:
	'''
	Only SNPs
	IUPAC codes
	'''
	input:
		genomeFA = rules.mask_callable_regions_from_genome.output.MaskedGenome,
		VCF = rules.r9_run_all_chr.output.passSNPsVCF
	output:
		ConsensusSeq="results/VariantAnalysis_DNA/Build_ConsensusSequences_PerSample/ConsensusSequence_{sample}.fa.gz",
		ConcatenatedConsensusSeq="results/VariantAnalysis_DNA/Build_ConsensusSequences_PerSample/ConsensusSequence_concatenated_{sample}.fa.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Build_ConsensusSequences_PerSample/ConsensusSequence_{sample}.err",
		out = "logs/VariantAnalysis_DNA/Build_ConsensusSequences_PerSample/ConsensusSequence_{sample}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Build_ConsensusSequences_PerSample/ConsensusSequence_{sample}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '10:00:00',
		name = "CS{sample}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantAnalysis_DNA/Build_ConsensusSequences_PerSample.sh {input.genomeFA} {input.VCF} {output.ConsensusSeq} {output.ConcatenatedConsensusSeq} {wildcards.sample} > {log.out} 2> {log.err}"

rule Join_ConsensusSequences:
	'''
	Make sure all the sequences have the same length
	Print them in a single file as aligned sequences
	'''
	input:
		ConcatenatedConsensusSeqs = expand(rules.Build_ConsensusSequences_PerSample.output.ConcatenatedConsensusSeq, sample=config["samples"])
	output:
		JoinedConsensusSeqs = "results/VariantAnalysis_DNA/Join_ConsensusSequences/WholeGenomeAlignmentConsensusSeqs.fa.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Join_ConsensusSequences/Join_ConsensusSequences.err",
		out = "logs/VariantAnalysis_DNA/Join_ConsensusSequences/Join_ConsensusSequences.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Join_ConsensusSequences/Join_ConsensusSequences.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '20:00:00',
		name = "jCS",
		threads = 10,
		mem = 100000
	shell:
		"./scripts/VariantAnalysis_DNA/Join_ConsensusSequences.sh \"{input.ConcatenatedConsensusSeqs}\" {output.JoinedConsensusSeqs} > {log.out} 2> {log.err}"

rule Build_DistanceTree_Samples:
	'''
	'''
	input:
		ConsensusSeq = rules.Join_ConsensusSequences.output.JoinedConsensusSeqs
	output:
		tree = "results/VariantAnalysis_DNA/Build_DistanceTree_Samples/Build_DistanceTree_Samples.treefile"
	log:
		err = "logs/VariantAnalysis_DNA/Build_DistanceTree_Samples/Build_DistanceTree_Samples.err",
		out = "logs/VariantAnalysis_DNA/Build_DistanceTree_Samples/Build_DistanceTree_Samples.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Build_DistanceTree_Samples/Build_DistanceTree_Samples.txt"
	conda:
		'../envs/Phylogenetics.yaml'
	params:
		time = '40:00:00',
		name = "treeSamp",
		threads = 10,
		mem = 100000
	shell:
		"""
		outfile={output.tree}
		rm -r $(dirname  ${{outfile}})
		mkdir -p $(dirname  ${{outfile}})
		iqtree -redo -nt {params.threads} -s {input.ConsensusSeq} -m GTR --prefix ${{outfile%.*}} > {log.out} 2> {log.err}
		"""

rule PSMC_PerSample:
	'''
	PSMC Pairwise Sequentially Markovian Coalescent
	Li H, Durbin R. Inference of human population history from individual whole-genome sequences. Nature. 2011 Jul 13;475(7357):493-6. doi: 10.1038/nature10231. PMID: 21753753; PMCID: PMC3154645.
	Code from Heng Li's GitHub: https://github.com/lh3/psmc downloaded on 05/05/2023

	Build only with SNPs (no INDELs)
	'''
	input:
		ConsensusSeq = rules.Build_ConsensusSequences_PerSample.output.ConcatenatedConsensusSeq
	output:
		psmc="results/VariantAnalysis_DNA/PSMC_PerSample/{sample}.psmc"
	log:
		err = "logs/VariantAnalysis_DNA/PSMC_PerSample/{sample}_PSMC_PerSample.err",
		out = "logs/VariantAnalysis_DNA/PSMC_PerSample/{sample}_PSMC_PerSample.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/PSMC_PerSample/{sample}_PSMC_PerSample.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "PSMC{sample}",
		threads = 1,
		mem = 10000,
		fq2psmcfa = "scripts/VariantAnalysis_DNA/psmc-master/utils/fq2psmcfa",
		psmc = "scripts/VariantAnalysis_DNA/psmc-master/psmc",
		psmc2history = "scripts/VariantAnalysis_DNA/psmc-master/utils/psmc2history.pl",
		history2ms = "scripts/VariantAnalysis_DNA/psmc-master/utils/history2ms.pl",
		psmc_plot = "scripts/VariantAnalysis_DNA/psmc-master/utils/psmc_plot.pl",
		Blangenerationtime = config["Blangenerationtime"]
	shell:
		"./scripts/VariantAnalysis_DNA/PSMC_PerSample.sh {input.ConsensusSeq} {output.psmc} {wildcards.sample} {params.fq2psmcfa} {params.psmc} {params.psmc2history} {params.history2ms} {params.psmc_plot} > {log.out} 2> {log.err}"










rule calling_stats_GATK_PerChr:
	'''
	Calcualte and print the calling stats in a table format with VariantsToTable (GATK)
	'''
	input:
		chrVCF = "results/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/ShortVariants_HardCallableFiltered.{chr}.vcf.gz"
	output:
		statsTXT = "results/VariantAnalysis_DNA/calling_stats_GATK_PerChr/ShortVariant_FilterStats.{chr}.txt"
	log:
		err = "logs/VariantAnalysis_DNA/calling_stats_GATK_PerChr/calling_stats_GATK_PerChr.{chr}.err",
		out = "logs/VariantAnalysis_DNA/calling_stats_GATK_PerChr/calling_stats_GATK_PerChr.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/calling_stats_GATK_PerChr/calling_stats_GATK_PerChr.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '03:00:00',
		name = "V2T{chr}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantAnalysis_DNA/calling_stats_GATK_PerChr.sh {input.chrVCF} {output.statsTXT} > {log.out} 2> {log.err}"

rule GONE_Recent_Demographic_History:
	'''
	GONE analysis of Recent Demographic History
	Santiago, E., Novo, I., Pardiñas, A. F. Saura, M., Wang, J., Caballero, A. (2020). Recent demographic history inferred by high-resolution analysis of linkage disequilibrium. Molecular Biology and Evolution 37: 3642–3653.
	Code from https://github.com/esrud/GONE downloaded on 07/09/2023
	'''
	input:
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa",
		VCF = "results/VariantCalling_DNA/9_run_all_chr/ShortVariants_HardCallableFiltered.vcf.gz"
	output:
		goneout="results/VariantAnalysis_DNA/GONE_Recent_Demographic_History/OUTPUT_GONE_Recent_Demographic_History"
	log:
		err = "logs/VariantAnalysis_DNA/GONE_Recent_Demographic_History/GONE.err",
		out = "logs/VariantAnalysis_DNA/GONE_Recent_Demographic_History/GONE.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/GONE_Recent_Demographic_History/GONE.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '20:00:00',
		name = "GONE",
		threads = 1,
		mem = 5000000,
		gone="scripts/VariantAnalysis_DNA/GONE/script_GONE_debug.sh",
		gone_INPUT_params="scripts/VariantAnalysis_DNA/GONE/INPUT_PARAMETERS_FILE",
		gone_PROGRAM_folder="scripts/VariantAnalysis_DNA/GONE/PROGRAMMES"
	shell:
		"./scripts/VariantAnalysis_DNA/GONE_Recent_Demographic_History.sh {input.genomeFA} {input.VCF} {output.goneout} {params.gone} {params.gone_INPUT_params} {params.gone_PROGRAM_folder} > {log.out} 2> {log.err}"




#################################################
## Deprecated rules


