

rule r2_map_dna_to_reference_minimap2_PerSample:
	'''
	Map to reference, sort per sample (minimap2 & samtools sort)
	'''
	input:
		sampleFASTQ1 = "data/DNAseqFASTQ/{sample}_{lane}_{platform}_R1.fastq.gz",
		sampleFASTQ2 = "data/DNAseqFASTQ/{sample}_{lane}_{platform}_R2.fastq.gz",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa",
		genomeFAI = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa.fai",
	output:
		sampleBAM = "results/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample/{sample}_{lane}_{platform}.bam",
	log:
		err = "logs/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample/{sample}_{lane}_{platform}.err",
		out = "logs/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample/{sample}_{lane}_{platform}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample/{sample}_{lane}_{platform}.txt"
	conda:
		'../envs/AlternativeVariantCalling_DNA.yaml'
	params:
		time = '10:00:00',
		name = "minimap2{sample}",
		threads = 4,
		mem = 15000
	shell:
		"./scripts/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample.sh {input.genomeFA} {input.sampleFASTQ1} {input.sampleFASTQ2} {output.sampleBAM} {params.threads} > {log.out} 2> {log.err}"
#
def get_bam_to_add_read_groups(wildcards):
	sample={wildcards.sample}
	lane={wildcards.lane}
	platform={wildcards.platform}
	if wildcards.mapper == "minimap2":
		return rules.r2_map_dna_to_reference_minimap2_PerSample.output.sampleBAM
	elif wildcards.mapper == "bwa_mem2":
		return "results/VariantCalling_DNA/2_map_dna_to_reference_BWA_PerSample/{sample}_{lane}_{platform}_sorted_markdup.bam"
	else:
		raise ValueError("Unknown value for mapper file in get_bam_to_add_read_groups: %s" % wildcards.mapper)

rule r2_add_read_groups_after_mapper_PerSample:
	'''
	run gatk AddOrReplaceReadGroups after minimap2 mapping, to add read groups information to the bam files, which is required for GATK BaseRecalibrator and HaplotypeCaller
	'''
	input:
		sampleBAM = get_bam_to_add_read_groups,
	output:
		sampleBAM_wrg = "results/AlternativeVariantCalling_DNA/2_add_read_groups_after_{mapper}_PerSample/{sample}_{lane}_{platform}_wrg.bam",
	log:
		err = "logs/AlternativeVariantCalling_DNA/2_add_read_groups_after_{mapper}_PerSample/{sample}_{lane}_{platform}.err",
		out = "logs/AlternativeVariantCalling_DNA/2_add_read_groups_after_{mapper}_PerSample/{sample}_{lane}_{platform}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/2_add_read_groups_after_{mapper}_PerSample/{sample}_{lane}_{platform}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '2:00:00',
		name = "addRG{sample}{mapper}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/AlternativeVariantCalling_DNA/2_add_read_groups_after_mapper_PerSample.sh {input.sampleBAM} {output.sampleBAM_wrg} {wildcards.sample} {wildcards.lane} {wildcards.platform} > {log.out} 2> {log.err}"
#
rule r2_split_bam_PerChr_PerSample_from_minimap2:
	'''
	Split sample bams per chr with samtools
	'''
	input:
		sampleBAM = expand(rules.r2_add_read_groups_after_mapper_PerSample.output.sampleBAM_wrg, sample="{sample}", lane="{lane}", platform="{platform}", mapper="minimap2")
	output:
		samplechrBAM = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/{sample}_{lane}_{platform}_round0.{chr}.bam"
	log:
		err = "logs/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample_from_minimap2/2_split_bam_PerChr_PerSample_{sample}_{lane}_{platform}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample_from_minimap2/2_split_bam_PerChr_PerSample_{sample}_{lane}_{platform}.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample_from_minimap2/2_split_bam_PerChr_PerSample_{sample}_{lane}_{platform}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '1:00:00',
		name = "S{chr}{sample}_{lane}_{platform}",
		threads = 1,
		mem = 4000
	shell:
		"mkdir -p $(dirname {output.samplechrBAM}); samtools view -b {input.sampleBAM} {wildcards.chr} > {output.samplechrBAM} 2> {log.err}"
#
def get_samplelaneplatform_minimap2_BAMs(round, sample, lane, platform, folder, chr):
	if round == "round0":
		return expand(rules.r2_split_bam_PerChr_PerSample_from_minimap2.output.samplechrBAM, sample=sample, lane=lane, platform=platform, chr=chr)
	elif round == "round1":
		return expand(rules.r3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2.output.samplechrplatformBAMRecal, sample=sample, lane=lane, platform=platform, chr=chr)
	elif round == "round2":
		return expand(rules.r3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2.output.samplechrplatformBAMRecal, sample=sample, lane=lane, platform=platform, chr=chr)
	else:
		raise ValueError("Unknown value for round file in get_samplelaneplatform_BAMs: %s" % wildcards.round)
def get_hiseq_minimap2_bams(wildcards):
    return [
        bam
        for s in config["samplelanesHiSeq4000"]
        for parts in [s.split("_")]
        for bam in get_samplelaneplatform_minimap2_BAMs(
            round=wildcards.round,
            sample=parts[0],
            lane=parts[1],
            platform="HiSeq4000",
            folder=wildcards.folder,
            chr=wildcards.chr
        )
    ]
def get_novaseq1_minimap2_bams(wildcards):
    return [
        bam
        for s in config["samples"]
        for parts in [s.split("_")]
        for bam in get_samplelaneplatform_minimap2_BAMs(
            round=wildcards.round,
            sample=s,
			lane="L1",
            platform="NovaSeq6000",
            folder=wildcards.folder,
            chr=wildcards.chr
        )
    ]
def get_novaseq2_minimap2_bams(wildcards):
    return [
        bam
        for s in config["samples"]
        for parts in [s.split("_")]
        for bam in get_samplelaneplatform_minimap2_BAMs(
            round=wildcards.round,
            sample=s,
			lane="L2",
            platform="NovaSeq6000",
            folder=wildcards.folder,
            chr=wildcards.chr
        )
    ]
rule r2_merge_3_lanes_bams_markdup_PerChr_from_minimap2:
	'''
	Sort by read name the three lane bams & merge them while marking duplicates per chr (samtools sort -n & gatk MarkDuplicatesSpark)
	WARNING! This rule cannot run more than 15 parallel tasks at a time.
	'''
	input:
		sampleschrBAMHiSeq = get_hiseq_minimap2_bams,
		sampleschrBAMNovaSeq1 = get_novaseq1_minimap2_bams,
		sampleschrBAMNovaSeq2 = get_novaseq2_minimap2_bams
	output:
		sampleschrBAMmarkdup = expand("results/AlternativeVariantCalling_DNA/{{folder}}/{sample}_{{round}}_merged.{{chr}}.bam", sample=config["samples"])
	wildcard_constraints:
		folder="3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2|3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2|3_base_recalibration_GATK_PerChr_PerSample_round3_from_minimap2|Merged_BAMs_for_freebayes_from_minimap2"
	log:
		err = "logs/AlternativeVariantCalling_DNA/{folder}/2_merge_3_lanes_bams_markdup_PerChr_{round}_{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/{folder}/2_merge_3_lanes_bams_markdup_PerChr_{round}_{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/{folder}/2_merge_3_lanes_bams_markdup_PerChr_{round}_{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '15:00:00',
		name = "M3{chr}",
		threads = 1,
		mem = 5000,
		samples = expand("{sample}", sample=config["samples"])
	resources:
		potatoes=10
	shell:
		"./scripts/VariantCalling_DNA/2_merge_3_lanes_bams_markdup_PerChr.sh \"{input.sampleschrBAMHiSeq}\" \"{input.sampleschrBAMNovaSeq1}\" \"{input.sampleschrBAMNovaSeq2}\" \"{output.sampleschrBAMmarkdup}\" \"{params.samples}\" > {log.out} 2> {log.err}"
#
rule r4_call_short_variants_GATK_PerChr_PerSample_from_minimap2:
	'''
	Calling short variants for a given chr from a merged bam file for the 3 lanes, sorted & markduplicated with GATK HaplotypeCaller
	'''
	input:
		samplechrBAMmarkdup = "results/AlternativeVariantCalling_DNA/{folder}/{sample}_{round}_merged.{chr}.bam",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		samplechrGVCF = "results/AlternativeVariantCalling_DNA/{folder}/{sample}_{round}.{chr}.chrg.vcf.gz"
	log:
		err = "logs/AlternativeVariantCalling_DNA/{folder}/4_call_short_variants_GATK_PerChr_PerSample_{sample}_{round}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/{folder}/4_call_short_variants_GATK_PerChr_PerSample_{sample}_{round}.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/{folder}/4_call_short_variants_GATK_PerChr_PerSample_{sample}_{round}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '20:00:00',
		name = "VC{chr}{sample}",
		threads = 4,
		mem = 10000
	shell:
		"./scripts/VariantCalling_DNA/4_call_short_variants_GATK_chromosome.sh {wildcards.chr} {input.samplechrBAMmarkdup} {input.genomeFA} {output.samplechrGVCF} {params.threads} > {log.out} 2> {log.err}"
rule r5_join_sample_calls_PerChr_from_minimap2:
	'''
	Joining sample short variants calls with GATK GenomicsDBImport & GenotypeGVCFs per each chr
	'''
	input:
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa",
		samplesGVCF = expand("results/AlternativeVariantCalling_DNA/{{folder}}/{sample}_{{round}}.{{chr}}.chrg.vcf.gz", sample=config["samples"])
	output:
		SamplesMap = "results/AlternativeVariantCalling_DNA/{folder}/SamplesMap_{round}.{chr}.txt",
		VCF = "results/AlternativeVariantCalling_DNA/{folder}/ShortVariants_{round}.{chr}.vcf.gz"
	wildcard_constraints:
		round = 'round[0-9]+'
	log:
		err = "logs/AlternativeVariantCalling_DNA/{folder}/5_join_sample_calls_PerChr_{round}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/{folder}/5_join_sample_calls_PerChr_{round}.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/{folder}/5_join_sample_calls_PerChr_{round}.{chr}.txt"
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
rule r6_hardfilter_VCF_PerChr_from_minimap2:
	'''
	Hard filter VCF per each chr
	'''
	input:
		chrVCF = rules.r5_join_sample_calls_PerChr_from_minimap2.output.VCF
	output:
		chrVCFHardFiltered = "results/AlternativeVariantCalling_DNA/{folder}/ShortVariants_HardFiltered_{round}.{chr}.vcf.gz"
	wildcard_constraints:
		round = 'round[0-9]+'
	log:
		err = "logs/AlternativeVariantCalling_DNA/{folder}/6_hardfilter_VCF_PerChr_{round}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/{folder}/6_hardfilter_VCF_PerChr_{round}.{chr}.out",
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/{folder}/6_hardfilter_VCF_PerChr_{round}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '3:00:00',
		name = "HF{chr}",
		threads = 1,
		mem = 4000
	shell:
		"./scripts/VariantCalling_DNA/6_hardfilter_VCF_PerChr.sh {input.chrVCF} {output.chrVCFHardFiltered} > {log.out} 2> {log.err}"

rule r3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2:
	'''
 	Round 1 of base recalibration per chr, sample & platform
 	'''
	input:
		chrVCFHardFiltered = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/ShortVariants_HardFiltered_round0.{chr}.vcf.gz",
		samplechrplatformBAM = rules.r2_split_bam_PerChr_PerSample_from_minimap2.output.samplechrBAM,
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		BRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/BaseRecalibration_GATK_{sample}_{lane}_{platform}_round1.{chr}.txt",
		samplechrplatformBAMRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/{sample}_{lane}_{platform}_round1.{chr}.bam"
	log:
		err = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{sample}_{lane}_{platform}_round1.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{sample}_{lane}_{platform}_round1.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{sample}_{lane}_{platform}_round1.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BR1{chr}_{sample}_{lane}_{platform}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample.sh {input.chrVCFHardFiltered} {input.samplechrplatformBAM} {input.genomeFA} {output.BRecal} {output.samplechrplatformBAMRecal} > {log.out} 2> {log.err}"
rule r3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2:
	'''
 	Round 2 of base recalibration per chr, sample & platform
 	'''
	input:
		chrVCFHardFiltered = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/ShortVariants_HardFiltered_round1.{chr}.vcf.gz",
		samplechrplatformBAM = rules.r3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2.output.samplechrplatformBAMRecal,
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		BRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/BaseRecalibration_GATK_{sample}_{lane}_{platform}_round2.{chr}.txt",
		samplechrplatformBAMRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round3_from_minimap2/{sample}_{lane}_{platform}_round2.{chr}.bam"
	log:
		err = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{sample}_{lane}_{platform}_round2.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{sample}_{lane}_{platform}_round2.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{sample}_{lane}_{platform}_round2.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BR2{chr}_{sample}_{lane}_{platform}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample.sh {input.chrVCFHardFiltered} {input.samplechrplatformBAM} {input.genomeFA} {output.BRecal} {output.samplechrplatformBAMRecal} > {log.out} 2> {log.err}"
def get_AlternativeVCFFiles(wildcards):
	if wildcards.mapper == "minimap2" and wildcards.caller == "freebayes":
		return expand(rules.r4_call_short_variants_freebayes_PerChr_PerSample.output.chrVCF, mapper="minimap2", chr=wildcards.chr)
	elif wildcards.mapper == "bwa_mem2" and wildcards.caller == "freebayes":
		return expand(rules.r4_call_short_variants_freebayes_PerChr_PerSample.output.chrVCF, mapper="bwa_mem2", chr=wildcards.chr)
	elif wildcards.mapper == "minimap2" and wildcards.caller == "GATK":
		return expand("results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round3_from_minimap2/ShortVariants_HardFiltered_round2.{chr}.vcf.gz", chr=wildcards.chr)
	else:
		raise ValueError("Unknown value for caller or mapper file in get_AlternativeVCFFiles: %s, %s" % (wildcards.caller, wildcards.mapper))
rule r8_filter_VCF_for_callable_regions_PerChr_from_AlternativeVariantCallings:
	'''
	Apply callable regions filter to the already hard filtered VCF per each chr
	'''
	input:
		NonCallableRegions = "results/VariantCalling_DNA/7_callable_regions/NonCallableRegions.bed.gz",
		chrVCFHardFiltered = get_AlternativeVCFFiles
	output:
		chrVCFHardCallableFiltered = "results/AlternativeVariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr_from_AlternativeVariantCallings/ShortVariants_HardCallableFiltered.{mapper}_{caller}.{chr}.vcf.gz",
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF_from_AlternativeVariantCallings.{mapper}_{caller}.{chr}.txt" 
	log:
		err = "logs/AlternativeVariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr_from_AlternativeVariantCallings/8_filter_VCF_for_callable_regions_PerChr.{mapper}_{caller}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr_from_AlternativeVariantCallings/8_filter_VCF_for_callable_regions_PerChr.{mapper}_{caller}.{chr}.out",
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr_from_AlternativeVariantCallings/8_filter_VCF_for_callable_regions_PerChr.{mapper}_{caller}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '3:00:00',
		name = "FVCF{chr}_{mapper}_{caller}",
		threads = 1,
		mem = 4000
	shell:
		"./scripts/VariantCalling_DNA/8_filter_VCF_for_callable_regions_PerChr.sh {input.NonCallableRegions} {input.chrVCFHardFiltered} {output.chrVCFHardCallableFiltered} {output.SamplesOrderInVCF} {wildcards.chr} > {log.out} 2> {log.err}"
rule r8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings:
	'''
	Apply extra callable regions filter to the already callable-filtered VCF per each chr
	See r7_extra_callable_regions for explanation on extra-callable regions
	'''
	input:
		NonExtraCallableRegions = "results/VariantCalling_DNA/7_join_extra_callable_regions/NonExtraCallableRegions.bed.gz",
		NonCallableRegions = "results/VariantCalling_DNA/7_callable_regions/NonCallableRegions.bed.gz",
		SNPableRegions ="results/VariantCalling_DNA/7_callable_regions/SNPableRegions.bed.gz",
		chrVCFHardCallableFiltered = rules.r8_filter_VCF_for_callable_regions_PerChr_from_AlternativeVariantCallings.output.chrVCFHardCallableFiltered
	output:
		chrVCFHardExtraCallableFiltered = "results/AlternativeVariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings/ShortVariants_HardCallableFiltered.{mapper}_{caller}.{chr}.vcf.gz"
	log:
		err = "logs/AlternativeVariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings/8_filter_VCF_for_extra_callable_regions_PerChr.{mapper}_{caller}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings/8_filter_VCF_for_extra_callable_regions_PerChr.{mapper}_{caller}.{chr}.out",
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings/8_filter_VCF_for_extra_callable_regions_PerChr.{mapper}_{caller}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '3:00:00',
		name = "eFVCF{chr}_{mapper}_{caller}",
		threads = 1,
		mem = 4000
	shell:
		"./scripts/VariantCalling_DNA/8_filter_VCF_for_extra_callable_regions_PerChr.sh {input.NonExtraCallableRegions} {input.NonCallableRegions} {input.SNPableRegions} {input.chrVCFHardCallableFiltered} {output.chrVCFHardExtraCallableFiltered} {wildcards.chr} > {log.out} 2> {log.err}"
#
rule r2_split_bam_PerChr_PerSample_from_bwa_mem2:
	'''
	Split sample bams per chr with samtools
	'''
	input:
		sampleBAMmarkdup = "results/AlternativeVariantCalling_DNA/2_add_read_groups_after_bwa_mem2_PerSample/{sample}_{lane}_{platform}_wrg.bam"
	output:
		samplechrBAMmarkdup = "results/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample_from_bwa_mem2/{sample}_{lane}_{platform}_round0.{chr}.bam"
	log:
		err = "logs/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample/2_split_bam_PerChr_PerSample_{sample}_{lane}_{platform}.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample/2_split_bam_PerChr_PerSample_{sample}_{lane}_{platform}.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample/2_split_bam_PerChr_PerSample_{sample}_{lane}_{platform}.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '1:00:00',
		name = "S{chr}_{sample}_{lane}_{platform}",
		threads = 1,
		mem = 4000
	shell:
		"samtools index {input.sampleBAMmarkdup}; mkdir -p $(dirname {output.samplechrBAMmarkdup}); samtools view -b {input.sampleBAMmarkdup} {wildcards.chr} > {output.samplechrBAMmarkdup} 2> {log.err}"
#
def get_samplelaneplatform_bwa_mem2_BAMs(round, sample, lane, platform, folder, chr):
		return expand(rules.r2_split_bam_PerChr_PerSample_from_bwa_mem2.output.samplechrBAMmarkdup, sample=sample, lane=lane, platform=platform, chr=chr)
def get_hiseq_bwa_mem2_bams(wildcards):
    return [
        bam
        for s in config["samplelanesHiSeq4000"]
        for parts in [s.split("_")]
        for bam in get_samplelaneplatform_bwa_mem2_BAMs(
            round=wildcards.round,
            sample=parts[0],
            lane=parts[1],
            platform="HiSeq4000",
            folder="Merged_BAMs_for_freebayes_from_bwa_mem2",
            chr=wildcards.chr
        )
    ]
def get_novaseq1_bwa_mem2_bams(wildcards):
    return [
        bam
        for s in config["samples"]
        for parts in [s.split("_")]
        for bam in get_samplelaneplatform_bwa_mem2_BAMs(
            round=wildcards.round,
            sample=s,
			lane="L1",
            platform="NovaSeq6000",
            folder="Merged_BAMs_for_freebayes_from_bwa_mem2",
            chr=wildcards.chr
        )
    ]
def get_novaseq2_bwa_mem2_bams(wildcards):
    return [
        bam
        for s in config["samples"]
        for parts in [s.split("_")]
        for bam in get_samplelaneplatform_bwa_mem2_BAMs(
            round=wildcards.round,
            sample=s,
			lane="L2",
            platform="NovaSeq6000",
            folder="Merged_BAMs_for_freebayes_from_bwa_mem2",
            chr=wildcards.chr
        )
    ]
rule r2_merge_3_lanes_bams_markdup_PerChr_from_bwa_mem2:
	'''
	Sort by read name the three lane bams & merge them while marking duplicates per chr (samtools sort -n & gatk MarkDuplicatesSpark)
	WARNING! This rule cannot run more than 15 parallel tasks at a time.
	'''
	input:
		sampleschrBAMHiSeq = get_hiseq_bwa_mem2_bams,
		sampleschrBAMNovaSeq1 = get_novaseq1_bwa_mem2_bams,
		sampleschrBAMNovaSeq2 = get_novaseq2_bwa_mem2_bams
	output:
		sampleschrBAMmarkdup = expand("results/AlternativeVariantCalling_DNA/Merged_BAMs_for_freebayes_from_bwa_mem2/{sample}_{{round}}_merged.{{chr}}.bam", sample=config["samples"])
	log:
		err = "logs/AlternativeVariantCalling_DNA/Merged_BAMs_for_freebayes_from_bwa_mem2/2_merge_3_lanes_bams_markdup_PerChr_{round}_{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/Merged_BAMs_for_freebayes_from_bwa_mem2/2_merge_3_lanes_bams_markdup_PerChr_{round}_{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/Merged_BAMs_for_freebayes_from_bwa_mem2/2_merge_3_lanes_bams_markdup_PerChr_{round}_{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '15:00:00',
		name = "M3{chr}",
		threads = 1,
		mem = 5000,
		samples = expand("{sample}", sample=config["samples"])
	resources:
		potatoes=10
	shell:
		"./scripts/VariantCalling_DNA/2_merge_3_lanes_bams_markdup_PerChr.sh \"{input.sampleschrBAMHiSeq}\" \"{input.sampleschrBAMNovaSeq1}\" \"{input.sampleschrBAMNovaSeq2}\" \"{output.sampleschrBAMmarkdup}\" \"{params.samples}\" > {log.out} 2> {log.err}"
#
def get_bam_for_freebayes(mapper, sample, round, folder, chr):
	if mapper == "minimap2":
		return expand(rules.r2_merge_3_lanes_bams_markdup_PerChr_from_minimap2.output.sampleschrBAMmarkdup, sample=sample, chr=chr, folder=folder, round=round)
	elif mapper == "bwa_mem2":
		return expand(rules.r2_merge_3_lanes_bams_markdup_PerChr_from_bwa_mem2.output.sampleschrBAMmarkdup, sample=sample, chr=chr, folder=folder, round=round)
	else:
		raise ValueError("Unknown value for BED file in get_bed_file_name: %s" % mapper)
def get_bams_for_freebayes_for_all_samples_bams(wildcards):
    return [
        bam
        for bam in get_bam_for_freebayes(
            mapper=wildcards.mapper,
            sample=config["samples"][1],
			round="round0",
            folder=f"Merged_BAMs_for_freebayes_from_{wildcards.mapper}",
            chr=wildcards.chr
        )
    ]

rule r4_call_short_variants_freebayes_PerChr_PerSample:
	'''
	Calling short variants for a given chr from a merged bam file for the 3 lanes with freebayes
	'''
	input:
		sampleschrBAMs = get_bams_for_freebayes_for_all_samples_bams,		
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		chrVCF = "results/AlternativeVariantCalling_DNA/4_call_short_variants_freebayes_PerChr_PerSample/{chr}.freebayes_from_{mapper}.vcf.gz"
	log:
		err = "logs/AlternativeVariantCalling_DNA/4_call_short_variants_freebayes_PerChr_PerSample/{chr}_from_{mapper}.err",
		out = "logs/AlternativeVariantCalling_DNA/4_call_short_variants_freebayes_PerChr_PerSample/{chr}_from_{mapper}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/4_call_short_variants_freebayes_PerChr_PerSample/{chr}_from_{mapper}.txt"
	conda:
		'../envs/AlternativeVariantCalling_DNA.yaml'
	params:
		time = '60:00:00',
		name = "freeb{chr}_{mapper}",
		threads = 1,
		mem = 10000
	shell:
		"./scripts/AlternativeVariantCalling_DNA/4_call_short_variants_freebayes_PerChr_PerSample.sh {wildcards.chr} \"{input.sampleschrBAMs}\" {input.genomeFA} {output.chrVCF} > {log.out} 2> {log.err}"
# Provide vcf files according to the wildcards caller and mapper
def get_AlternativeFilteredVCFFiles(wildcards):
	if wildcards.mapper == "minimap2" and wildcards.caller == "freebayes":
		return expand(rules.r8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings.output.chrVCFHardExtraCallableFiltered, mapper="minimap2", caller="freebayes", chr=wildcards.chr)
	elif wildcards.mapper == "bwa_mem2" and wildcards.caller == "freebayes":
		return expand(rules.r8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings.output.chrVCFHardExtraCallableFiltered, mapper="bwa_mem2", caller="freebayes", chr=wildcards.chr)
	elif wildcards.mapper == "minimap2" and wildcards.caller == "GATK":
		return expand(rules.r8_filter_VCF_for_extra_callable_regions_PerChr_from_AlternativeVariantCallings.output.chrVCFHardExtraCallableFiltered, mapper="minimap2", caller="GATK", chr=wildcards.chr)
	elif wildcards.mapper == "bwa_mem2" and wildcards.caller == "GATK":
		return expand(rules.r8_filter_VCF_for_extra_callable_regions_PerChr.output.chrVCFHardExtraCallableFiltered, chr=wildcards.chr)
	else:
		raise ValueError("Unknown value for caller or mapper file in get_AlternativeVCFFiles: %s, %s" % (wildcards.caller, wildcards.mapper))
# Read VCF and build genotype, numericmatrix and 0hom1het files
rule Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings:
	'''
	Prepare files for further analysis:
		- .genotype.bed file
		- .genotype.numericmatrix file (#alleles-1 lines for each site)
		- .0hom1het.bed file
	'''
	input:
		chrVCF = get_AlternativeFilteredVCFFiles,
		SamplesOrderInVCF = rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF
	output:
		chrGENOTYPE = "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings/ShortVariants.filtered.{mapper}_{caller}.{chr}.genotype.bed.gz",
		chrGENOTYPEMatrix = "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings/ShortVariants.filtered.{mapper}_{caller}.{chr}.genotype.numericmatrix.gz",
		HeterozygosityFile = "results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings/ShortVariants.filtered.{mapper}_{caller}.{chr}.0hom1het.bed.gz"
	log:
		err = "logs/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings.{mapper}_{caller}.{chr}.err",
		out = "logs/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings.{mapper}_{caller}.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings/Prepare_basic_analysis_files_PerChr_from_AlternativeVariantCallings.{mapper}_{caller}.{chr}.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '5:00:00',
		name = "{mapper}_{caller}Genotype{chr}",
		threads = 1,
		mem = 25000
	shell:
		"./scripts/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr.sh {input.chrVCF} {output.chrGENOTYPE} > {log.out} 2> {log.err}"
#
