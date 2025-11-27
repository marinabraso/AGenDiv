

rule copy_fastq_NAS_data:
	'''
	Copying from NAS storage dna fastq folder to the working data folder
	'''
	input:
		NASfastqfile = "/nas/FAC/FBM/DEE/mrobinso/default/D2c/mbrasovi/Banyuls_Roscoff/DNAseqFASTQ/{idfile}.fastq.gz"
	output:
		WORKfastqfile = "data/DNAseqFASTQ/{idfile}.fastq.gz"
	log:
		err = "logs/VariantCalling_DNA/copy_fastq_NAS_data/{idfile}.err",
		out = "logs/VariantCalling_DNA/copy_fastq_NAS_data/{idfile}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/copy_fastq_NAS_data/{idfile}.txt"
	shell:
		"cp {input.NASfastqfile} {output.WORKfastqfile}"
rule r1_index_genome_for_BWA:
	'''
	Indexing genome in fasta format for mapping with BWA (bwa, samtools, picard)
	'''
	input:
		genomeFA = "data/genomes/{genome}/{genomefile}.fa"
	output:
		genomeFAI = "data/genomes/{genome}/{genomefile}.fa.fai",
		genomeDICT = "data/genomes/{genome}/{genomefile}.dict"
	log:
		err = "logs/VariantCalling_DNA/1_index_genome_for_BWA/{genome}_{genomefile}.err",
		out = "logs/VariantCalling_DNA/1_index_genome_for_BWA/{genome}_{genomefile}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/1_index_genome_for_BWA/{genome}_{genomefile}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '01:00:00',
		name = "Ix{genome}",
		threads = 1,
		mem = 2000
	shell:
		"./scripts/VariantCalling_DNA/1_index_genome_for_BWA.sh {input.genomeFA} {output.genomeDICT} > {log.out} 2> {log.err}"
rule r2_map_dna_to_reference_BWA_PerSample:
	'''
	Map to reference, sort & mark duplicates per sample (BWA MEM & samtools sort & gatk MarkDuplicatesSpark)
	'''
	input:
		sampleFASTQ1 = "data/DNAseqFASTQ/{sample}_{lane}_{platform}_R1.fastq.gz",
		sampleFASTQ2 = "data/DNAseqFASTQ/{sample}_{lane}_{platform}_R2.fastq.gz",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa",
		genomeFAI = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa.fai",
		readgroupsfile = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000_formated.txt"
	output:
		sampleBAM = "results/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{sample}_{lane}_{platform}_sorted.bam",
		sampleBAMmarkdup = "results/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{sample}_{lane}_{platform}_sorted_markdup.bam"
	log:
		err = "logs/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{sample}_{lane}_{platform}.err",
		out = "logs/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{sample}_{lane}_{platform}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{sample}_{lane}_{platform}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '10:00:00',
		name = "BWA{sample}",
		threads = 1,
		mem = 15000
	shell:
		"./scripts/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample.sh {wildcards.sample} {wildcards.lane} {wildcards.platform} {input.genomeFA} {input.sampleFASTQ1} {input.sampleFASTQ2} {input.readgroupsfile} {output.sampleBAM} {output.sampleBAMmarkdup} > {log.out} 2> {log.err}"
rule r2_split_bam_PerChr_PerSample:
	'''
	Split sample bams per chr with samtools
	'''
	input:
		sampleBAMmarkdup = "results/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{samplelaneplatform}_sorted_markdup.bam"
	output:
		samplechrBAMmarkdup = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/{samplelaneplatform}_round0.{chr}.bam"
	log:
		err = "logs/VariantCalling_DNA/2_split_bam_PerChr_PerSample/2_split_bam_PerChr_PerSample_{samplelaneplatform}.{chr}.err",
		out = "logs/VariantCalling_DNA/2_split_bam_PerChr_PerSample/2_split_bam_PerChr_PerSample_{samplelaneplatform}.{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/2_split_bam_PerChr_PerSample/2_split_bam_PerChr_PerSample_{samplelaneplatform}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '1:00:00',
		name = "S{chr}{samplelaneplatform}",
		threads = 1,
		mem = 4000
	shell:
		"mkdir -p $(dirname {output.samplechrBAMmarkdup}); samtools view -b {input.sampleBAMmarkdup} {wildcards.chr} > {output.samplechrBAMmarkdup} 2> {log.err}"
rule r2_merge_3_lanes_bams_markdup_PerChr:
	'''
	Sort by read name the three lane bams & merge them while marking duplicates per chr (samtools sort -n & gatk MarkDuplicatesSpark)
	WARNING! This rule cannot run more than 15 parallel tasks at a time.
	'''
	input:
		sampleschrBAMHiSeq = expand("results/VariantCalling_DNA/{{folder}}/{samplelaneHiSeq}_HiSeq4000_{{round}}.{{chr}}.bam", samplelaneHiSeq=config["samplelanesHiSeq4000"]),
		sampleschrBAMNovaSeq1 = expand("results/VariantCalling_DNA/{{folder}}/{sample}_L1_NovaSeq6000_{{round}}.{{chr}}.bam", sample=config["samples"]),
		sampleschrBAMNovaSeq2 = expand("results/VariantCalling_DNA/{{folder}}/{sample}_L2_NovaSeq6000_{{round}}.{{chr}}.bam", sample=config["samples"])
	output:
		sampleschrBAMmarkdup = expand("results/VariantCalling_DNA/{{folder}}/{sample}_{{round}}_merged.{{chr}}.bam", sample=config["samples"])
	log:
		err = "logs/VariantCalling_DNA/{folder}/2_merge_3_lanes_bams_markdup_PerChr_{chr}_{round}.err",
		out = "logs/VariantCalling_DNA/{folder}/2_merge_3_lanes_bams_markdup_PerChr_{chr}_{round}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/{folder}/2_merge_3_lanes_bams_markdup_PerChr_{chr}_{round}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '15:00:00',
		name = "M3{chr}{round}",
		threads = 1,
		mem = 5000,
		samples = expand("{sample}", sample=config["samples"])
	resources:
		potatoes=10
	shell:
		"./scripts/VariantCalling_DNA/2_merge_3_lanes_bams_markdup_PerChr.sh \"{input.sampleschrBAMHiSeq}\" \"{input.sampleschrBAMNovaSeq1}\" \"{input.sampleschrBAMNovaSeq2}\" \"{output.sampleschrBAMmarkdup}\" \"{params.samples}\" > {log.out} 2> {log.err}"
rule r3_base_recalibration_GATK_PerChr_PerSample_round1:
	'''
 	Round 1 of base recalibration per chr, sample & platform
 	'''
	input:
		chrVCFHardFiltered = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/ShortVariants_HardFiltered_round0.{chr}.vcf.gz",
		samplechrplatformBAM = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/{samplelaneplatform}_round0.{chr}.bam",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		BRecal = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/BaseRecalibration_GATK_{samplelaneplatform}_round1.{chr}.txt",
		samplechrplatformBAMRecal = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/{samplelaneplatform}_round1.{chr}.bam"
	log:
		err = "logs/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round1.{chr}.err",
		out = "logs/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round1.{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round1.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BR1{chr}{samplelaneplatform}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample.sh {input.chrVCFHardFiltered} {input.samplechrplatformBAM} {input.genomeFA} {output.BRecal} {output.samplechrplatformBAMRecal} > {log.out} 2> {log.err}"
rule r3_base_recalibration_GATK_PerChr_PerSample_round2:
	'''
 	Round 2 of base recalibration per chr, sample & platform
 	'''
	input:
		chrVCFHardFiltered = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/ShortVariants_HardFiltered_round1.{chr}.vcf.gz",
		samplechrplatformBAM = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/{samplelaneplatform}_round1.{chr}.bam",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		BRecal = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/BaseRecalibration_GATK_{samplelaneplatform}_round2.{chr}.txt",
		samplechrplatformBAMRecal = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round3/{samplelaneplatform}_round2.{chr}.bam"
	log:
		err = "logs/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round2.{chr}.err",
		out = "logs/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round2.{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round2.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BR2{chr}{samplelaneplatform}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample.sh {input.chrVCFHardFiltered} {input.samplechrplatformBAM} {input.genomeFA} {output.BRecal} {output.samplechrplatformBAMRecal} > {log.out} 2> {log.err}"
rule r3_base_recalibration_GATK_PerChr_Unification:
	'''
 	Base Recalibration Bootstrapping Unification per all samples
 	'''
	input:
		BRecalR1HiSeq = expand("results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/BaseRecalibration_GATK_{samplelaneplatform}_round1.{{chr}}.txt", samplelaneplatform=config["samplelanesHiSeq4000"]),
		BRecalR2HiSeq = expand("results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/BaseRecalibration_GATK_{samplelaneplatform}_round2.{{chr}}.txt", samplelaneplatform=config["samplelanesHiSeq4000"]),
		BRecalR1NovaSeq1 = expand("results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/BaseRecalibration_GATK_{samplelaneplatform}_round1.{{chr}}.txt", samplelaneplatform=config["samplelanesNovaSeq6000_1"]),
		BRecalR2NovaSeq1 = expand("results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/BaseRecalibration_GATK_{samplelaneplatform}_round2.{{chr}}.txt", samplelaneplatform=config["samplelanesNovaSeq6000_1"]),
		BRecalR1NovaSeq2 = expand("results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1/BaseRecalibration_GATK_{samplelaneplatform}_round1.{{chr}}.txt", samplelaneplatform=config["samplelanesNovaSeq6000_2"]),
		BRecalR2NovaSeq2 = expand("results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2/BaseRecalibration_GATK_{samplelaneplatform}_round2.{{chr}}.txt", samplelaneplatform=config["samplelanesNovaSeq6000_2"])
	output:
		BaseRecalTXT = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_Unification/BaseRecalibration_GATK_covariates.{chr}.txt"
	log:
		err = "logs/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_Unification/3_base_recalibration_GATK_PerChr_Unification_{chr}.err",
		out = "logs/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_Unification/3_base_recalibration_GATK_PerChr_Unification_{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_Unification/3_base_recalibration_GATK_PerChr_Unification_{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BRU{chr}",
		threads = 1,
		mem = 2000,
		samples = expand("{sample}", sample=config["samples"])
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_Unification.sh \"{input.BRecalR1HiSeq}\" \"{input.BRecalR2HiSeq}\" \"{input.BRecalR1NovaSeq1}\" \"{input.BRecalR2NovaSeq1}\" \"{input.BRecalR1NovaSeq2}\" \"{input.BRecalR2NovaSeq2}\" {output.BaseRecalTXT} \"{params.samples}\" {wildcards.chr} > {log.out} 2> {log.err}"
rule r4_call_short_variants_GATK_PerChr_PerSample:
	'''
	Calling short variants for a given chr from a merged bam file for the 3 lanes, sorted & markduplicated with GATK HaplotypeCaller
	'''
	input:
		samplechrBAMmarkdup = "results/VariantCalling_DNA/{folder}/{sample}_{round}_merged.{chr}.bam",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		samplechrGVCF = "results/VariantCalling_DNA/{folder}/{sample}_{round}.{chr}.chrg.vcf.gz"
	log:
		err = "logs/VariantCalling_DNA/{folder}/4_call_short_variants_GATK_PerChr_PerSample_{sample}_{round}.{chr}.err",
		out = "logs/VariantCalling_DNA/{folder}/4_call_short_variants_GATK_PerChr_PerSample_{sample}_{round}.{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/{folder}/4_call_short_variants_GATK_PerChr_PerSample_{sample}_{round}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '20:00:00',
		name = "VC{chr}{sample}",
		threads = 32,
		mem = 10000
	shell:
		"./scripts/VariantCalling_DNA/4_call_short_variants_GATK_chromosome.sh {wildcards.chr} {input.samplechrBAMmarkdup} {input.genomeFA} {output.samplechrGVCF} > {log.out} 2> {log.err}"
rule r5_join_sample_calls_PerChr:
	'''
	Joining sample short variants calls with GATK GenomicsDBImport & GenotypeGVCFs per each chr
	'''
	input:
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa",
		samplesGVCF = expand("results/VariantCalling_DNA/{{folder}}/{sample}_{{round}}.{{chr}}.chrg.vcf.gz", sample=config["samples"])
	output:
		SamplesMap = "results/VariantCalling_DNA/{folder}/SamplesMap_{round}.{chr}.txt",
		VCF = "results/VariantCalling_DNA/{folder}/ShortVariants_{round}.{chr}.vcf.gz"
	wildcard_constraints:
		round = 'round[0-9]+'
	log:
		err = "logs/VariantCalling_DNA/{folder}/5_join_sample_calls_PerChr_{round}.{chr}.err",
		out = "logs/VariantCalling_DNA/{folder}/5_join_sample_calls_PerChr_{round}.{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/{folder}/5_join_sample_calls_PerChr_{round}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '20:00:00',
		name = "JS{chr}",
		threads = 5,
		mem = 4000,
		samples = expand("{sample}", sample=config["samples"])
	shell:
		"./scripts/VariantCalling_DNA/5_join_sample_calls_PerChr.sh {wildcards.chr} {input.genomeFA} \"{input.samplesGVCF}\" {output.SamplesMap} {output.VCF} \"{params.samples}\" > {log.out} 2> {log.err}"
rule r6_hardfilter_VCF_PerChr:
	'''
	Hard filter VCF per each chr
	'''
	input:
		chrVCF = "results/VariantCalling_DNA/{folder}/ShortVariants_{round}.{chr}.vcf.gz"
	output:
		chrVCFHardFiltered = "results/VariantCalling_DNA/{folder}/ShortVariants_HardFiltered_{round}.{chr}.vcf.gz"
	wildcard_constraints:
		round = 'round[0-9]+'
	log:
		err = "logs/VariantCalling_DNA/{folder}/6_hardfilter_VCF_PerChr_{round}.{chr}.err",
		out = "logs/VariantCalling_DNA/{folder}/6_hardfilter_VCF_PerChr_{round}.{chr}.out",
	benchmark:
		"benchmarks/VariantCalling_DNA/{folder}/6_hardfilter_VCF_PerChr_{round}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '3:00:00',
		name = "HF{chr}",
		threads = 1,
		mem = 4000
	shell:
		"./scripts/VariantCalling_DNA/6_hardfilter_VCF_PerChr.sh {input.chrVCF} {output.chrVCFHardFiltered} > {log.out} 2> {log.err}"
rule r7_SNPable_Regions:
	'''
	Following SNPable pipeline http://lh3lh3.users.sourceforge.net/snpable.shtml
	Source code downloaded on 23/04/11 from https://lh3lh3.users.sourceforge.net/download/seqbility-20091110.tar.bz2
	'''
	input:
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		maskedgenomeFA = "results/VariantCalling_DNA/7_SNPable_Regions/Branchiostoma_lanceolatum.BraLan3_genome_masked.fa",
		SNPableQuality = "results/VariantCalling_DNA/7_SNPable_Regions/Branchiostoma_lanceolatum.BraLan3_mask_35_50.fa",
		unSNPableRegions = "results/VariantCalling_DNA/7_SNPable_Regions/Branchiostoma_lanceolatum.BraLan3_unSNPableReg_35_50.bed"
	log:
		err = "logs/VariantCalling_DNA/7_SNPable_Regions/7_SNPable_Regions.err",
		out = "logs/VariantCalling_DNA/7_SNPable_Regions/7_SNPable_Regions.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/7_SNPable_Regions/7_SNPable_Regions.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '15:00:00',
		name = "SNPable",
		threads = 4,
		mem = 25000,
		splitfa = "scripts/VariantCalling_DNA/SNPable/splitfa",
		gen_raw_mask = "scripts/VariantCalling_DNA/SNPable/gen_raw_mask.pl",
		gen_mask = "scripts/VariantCalling_DNA/SNPable/gen_mask",
		apply_mask_s = "scripts/VariantCalling_DNA/SNPable/apply_mask_s"
	shell:
		"./scripts/VariantCalling_DNA/7_SNPable_Regions.sh {input.genomeFA} {output.maskedgenomeFA} {params.splitfa} {params.gen_raw_mask} {params.gen_mask} {params.apply_mask_s} > {log.out} 2> {log.err}"
rule r7_coverage_per_site_PerChr:
	'''
	Coverage per site per chr
	'''
	input:
		chrBAMsmarkdup = expand("results/VariantCalling_DNA/7_coverage_per_site_PerChr/{sample}_round0_merged.{{chr}}.bam", sample=config["samples"])
	output:
		chrcovgenome = "results/VariantCalling_DNA/7_coverage_per_site_PerChr/Depth.allsamples.{chr}.tar.gz",
		chrmincov = "results/VariantCalling_DNA/7_coverage_per_site_PerChr/Depth.MinCov.{chr}.tar.gz"
	log:
		err = "logs/VariantCalling_DNA/7_coverage_per_site_PerChr/7_coverage_per_site_PerChr.{chr}.err",
		out = "logs/VariantCalling_DNA/7_coverage_per_site_PerChr/7_coverage_per_site_PerChr.{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/7_coverage_per_site_PerChr/7_coverage_per_site_PerChr.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '20:00:00',
		name = "Cov{chr}",
		threads = 1,
		mem = 10000,
		samples = expand("{sample}", sample=config["samples"])
	shell:
		"./scripts/VariantCalling_DNA/7_coverage_per_site_PerChr.sh \"{input.chrBAMsmarkdup}\" {output.chrcovgenome} {output.chrmincov} \"{params.samples}\" {wildcards.chr} > {log.out} 2> {log.err}"
rule r7_join_coverage_per_site:
	'''
	Join coverage per all chr
	'''
	input:
		mincovs = expand("results/VariantCalling_DNA/7_coverage_per_site_PerChr/Depth.MinCov.{chr}.tar.gz", chr=config["bralan3chrs"]),
		covgenomes = expand("results/VariantCalling_DNA/7_coverage_per_site_PerChr/Depth.allsamples.{chr}.tar.gz", chr=config["bralan3chrs"])
	output:
		mincov = "results/VariantCalling_DNA/7_join_coverage_per_site/Depth.MinCov.tbl.gz",
		covgenome = "results/VariantCalling_DNA/7_join_coverage_per_site/Depth.allsamples.tbl.gz"
	log:
		err = "logs/VariantCalling_DNA/7_join_coverage_per_site/7_join_coverage_per_site.err",
		out = "logs/VariantCalling_DNA/7_join_coverage_per_site/7_join_coverage_per_site.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/7_join_coverage_per_site/7_join_coverage_per_site.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '10:00:00',
		name = "jCov",
		threads = 1,
		mem = 50000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/VariantCalling_DNA/7_join_coverage_per_site.sh \"{input.mincovs}\" \"{input.covgenomes}\" {output.mincov} {output.covgenome} \"{params.chrs}\" > {log.out} 2> {log.err}"
rule r7_callable_regions:
	'''
	Callable = SNPable & coverage threshold
	'''
	input:
		SNPableQuality = "results/VariantCalling_DNA/7_SNPable_Regions/Branchiostoma_lanceolatum.BraLan3_mask_35_50.fa",
		mincov = "results/VariantCalling_DNA/7_join_coverage_per_site/Depth.MinCov.tbl.gz"
	output:
		CallableRegions = "results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz",
		NonCallableRegions = "results/VariantCalling_DNA/7_callable_regions/NonCallableRegions.bed.gz",
		SNPableRegions ="results/VariantCalling_DNA/7_callable_regions/SNPableRegions.bed.gz"
	log:
		err = "logs/VariantCalling_DNA/7_callable_regions/7_callable_regions.err",
		out = "logs/VariantCalling_DNA/7_callable_regions/7_callable_regions.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/7_callable_regions/7_callable_regions.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '15:00:00',
		name = "CallableR",
		threads = 1,
		mem = 50000,
		lowermincovthres = config["lowermincovthres"], # for a given position, the minimum coverage among all samples has to be at least lowermincovthres
		uppermincovthres = config["uppermincovthres"], # for a given position, the minimum coverage among all samples has to be less or equal to uppermincovthres
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/VariantCalling_DNA/7_callable_regions.sh {input.SNPableQuality} {input.mincov} {output.CallableRegions} {output.NonCallableRegions} {params.lowermincovthres} {params.uppermincovthres} \"{params.chrs}\" > {log.out} 2> {log.err}"
rule r7_extra_callable_regions_PerChr:
	'''
	Preparing extra filter for callable regions for a single chr: 
		Whole genome sliding windows have to have a minimum of its length being callable to be considered really callable (extra callable)
		Minimum prportion of the window being callable: slidingwindowmincallable
	'''
	input:
		CallableRegions = "results/VariantCalling_DNA/7_callable_regions/CallableRegions.bed.gz",
		ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
	output:
		ExtraCallableRegions = "results/VariantCalling_DNA/7_extra_callable_regions_PerChr/ExtraCallableRegions_{chr}.bed.gz",
		NonExtraCallableRegions = "results/VariantCalling_DNA/7_extra_callable_regions_PerChr/NonExtraCallableRegions_{chr}.bed.gz",
		ValidWindows = "results/VariantCalling_DNA/7_extra_callable_regions_PerChr/ValidWindows_{chr}.bed.gz"
	log:
		err = "logs/VariantCalling_DNA/7_extra_callable_regions_PerChr/7_extra_callable_regions_{chr}.err",
		out = "logs/VariantCalling_DNA/7_extra_callable_regions_PerChr/7_extra_callable_regions_{chr}.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/7_extra_callable_regions_PerChr/7_extra_callable_regions_{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '15:00:00',
		name = "eCal{chr}",
		threads = 1,
		mem = 50000,
		slidingwindowsize=config["slidingwindowsize"], # window size for sliding window calculations 
		slidingwindowstep=config["slidingwindowstep"], # step size for sliding window calculations
		slidingwindowmincallable=config["slidingwindowmincallable"], # minimum proportion of the sliding window being callable to be accepted as really callable
	shell:
		"./scripts/VariantCalling_DNA/7_extra_callable_regions_PerChr.sh {input.CallableRegions} {input.ChrLengths} {output.ExtraCallableRegions} {output.NonExtraCallableRegions} {output.ValidWindows} {params.slidingwindowsize} {params.slidingwindowstep} {params.slidingwindowmincallable} {wildcards.chr} > {log.out} 2> {log.err}"
rule r7_join_extra_callable_regions:
	'''
	Join all chrs extra callable regions
	'''
	input:
		chrsExtraCallableRegions = expand("results/VariantCalling_DNA/7_extra_callable_regions_PerChr/ExtraCallableRegions_{chr}.bed.gz", chr=config["bralan3chrs"]),
		chrsNonExtraCallableRegions = expand("results/VariantCalling_DNA/7_extra_callable_regions_PerChr/NonExtraCallableRegions_{chr}.bed.gz", chr=config["bralan3chrs"])
	output:
		ExtraCallableRegions = "results/VariantCalling_DNA/7_join_extra_callable_regions/ExtraCallableRegions.bed.gz",
		NonExtraCallableRegions = "results/VariantCalling_DNA/7_join_extra_callable_regions/NonExtraCallableRegions.bed.gz"
	log:
		err = "logs/VariantCalling_DNA/7_join_extra_callable_regions/7_extra_callable_regions.err",
		out = "logs/VariantCalling_DNA/7_join_extra_callable_regions/7_extra_callable_regions.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/7_join_extra_callable_regions/7_extra_callable_regions.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '15:00:00',
		name = "jeCal",
		threads = 1,
		mem = 50000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"scripts/VariantCalling_DNA/7_join_extra_callable_regions.sh \"{input.chrsExtraCallableRegions}\" \"{input.chrsNonExtraCallableRegions}\" {output.ExtraCallableRegions} {output.NonExtraCallableRegions} \"{params.chrs}\"  > {log.out} 2> {log.err}"
rule r8_filter_VCF_for_callable_regions_PerChr:
	'''
	Apply callable regions filter to the already hard filtered VCF per each chr
	'''
	input:
		NonCallableRegions = "results/VariantCalling_DNA/7_callable_regions/NonCallableRegions.bed.gz",
		chrVCFHardFiltered = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round3/ShortVariants_HardFiltered_round2.{chr}.vcf.gz"
	output:
		chrVCFHardCallableFiltered = "results/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/ShortVariants_HardCallableFiltered.{chr}.vcf.gz",
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.{chr}.txt" 
	log:
		err = "logs/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/8_filter_VCF_for_callable_regions_PerChr.{chr}.err",
		out = "logs/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/8_filter_VCF_for_callable_regions_PerChr.{chr}.out",
	benchmark:
		"benchmarks/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/8_filter_VCF_for_callable_regions_PerChr.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '3:00:00',
		name = "FVCF{chr}",
		threads = 1,
		mem = 4000
	shell:
		"./scripts/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr.sh {input.NonCallableRegions} {input.chrVCFHardFiltered} {output.chrVCFHardCallableFiltered} {output.SamplesOrderInVCF} {wildcards.chr} > {log.out} 2> {log.err}"
rule r8_filter_VCF_for_extra_callable_regions_PerChr:
	'''
	Apply extra callable regions filter to the already callable-filtered VCF per each chr
	See r7_extra_callable_regions for explanation on extra-callable regions
	'''
	input:
		NonExtraCallableRegions = "results/VariantCalling_DNA/7_join_extra_callable_regions/NonExtraCallableRegions.bed.gz",
		NonCallableRegions = "results/VariantCalling_DNA/7_callable_regions/NonCallableRegions.bed.gz",
		SNPableRegions ="results/VariantCalling_DNA/7_callable_regions/SNPableRegions.bed.gz",
		chrVCFHardCallableFiltered = "results/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr/ShortVariants_HardCallableFiltered.{chr}.vcf.gz"
	output:
		chrVCFHardExtraCallableFiltered = "results/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/ShortVariants_HardCallableFiltered.{chr}.vcf.gz"
	log:
		err = "logs/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/8_filter_VCF_for_extra_callable_regions_PerChr.{chr}.err",
		out = "logs/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/8_filter_VCF_for_extra_callable_regions_PerChr.{chr}.out",
	benchmark:
		"benchmarks/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/8_filter_VCF_for_extra_callable_regions_PerChr.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '3:00:00',
		name = "eFVCF{chr}",
		threads = 1,
		mem = 4000
	shell:
		"./scripts/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr.sh {input.NonExtraCallableRegions} {input.NonCallableRegions} {input.SNPableRegions} {input.chrVCFHardCallableFiltered} {output.chrVCFHardExtraCallableFiltered} {wildcards.chr} > {log.out} 2> {log.err}"
rule r9_run_all_chr:
	'''
	Create VCF file for all chr
	'''
	input:
		VCFs = expand("results/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr/ShortVariants_HardCallableFiltered.{chr}.vcf.gz", chr=config["bralan3chrs"])
	output:
		unifiedVCF = "results/VariantCalling_DNA/9_run_all_chr/ShortVariants_HardExtraCallableFiltered.vcf.gz",
		passVCF = "results/VariantCalling_DNA/9_run_all_chr/ShortVariants_HardExtraCallableFiltered_PASS.vcf.gz",
		passSNPsVCF = "results/VariantCalling_DNA/9_run_all_chr/ShortVariants_HardExtraCallableFiltered_PASS_SNPs.vcf.gz"
	log:
		err = "logs/VariantCalling_DNA/9_run_all_chr/9_run_all_chr.err",
		out = "logs/VariantCalling_DNA/9_run_all_chr/9_run_all_chr.out"
	benchmark:
		"benchmarks/VariantCalling_DNA/9_run_all_chr/9_run_all_chr.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '10:00:00',
		name = "jVCF",
		threads = 1,
		mem = 50000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/VariantCalling_DNA/9_run_all_chr.sh \"{input.VCFs}\" {output.unifiedVCF} {output.passVCF} {output.passSNPsVCF} > {log.out} 2> {log.err}"
















































##






