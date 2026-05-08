

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
		sampleBAM = "results/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample/{sample}_{lane}_{platform}_sorted.bam",
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
rule r2_split_bam_PerChr_PerSample_from_minimap2:
	'''
	Split sample bams per chr with samtools
	'''
	input:
		sampleBAM = "results/AlternativeVariantCalling_DNA/2_map_dna_to_reference_minimap2_PerSample/{sample}_{lane}_{platform}_sorted.bam"
	output:
		samplechrBAM = "results/AlternativeVariantCalling_DNA/2_split_bam_PerChr_PerSample_from_minimap2/{sample}_{lane}_{platform}_round0.{chr}.bam"
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
rule r2_merge_3_lanes_bams_markdup_PerChr_from_minimap2:
	'''
	Sort by read name the three lane bams & merge them while marking duplicates per chr (samtools sort -n & gatk MarkDuplicatesSpark)
	WARNING! This rule cannot run more than 15 parallel tasks at a time.
	'''
	input:
		sampleschrBAMHiSeq = expand("results/AlternativeVariantCalling_DNA/{{folder}}/{samplelaneHiSeq}_HiSeq4000.{{chr}}.bam", samplelaneHiSeq=config["samplelanesHiSeq4000"]),
		sampleschrBAMNovaSeq1 = expand("results/AlternativeVariantCalling_DNA/{{folder}}/{sample}_L1_NovaSeq6000.{{chr}}.bam", sample=config["samples"]),
		sampleschrBAMNovaSeq2 = expand("results/AlternativeVariantCalling_DNA/{{folder}}/{sample}_L2_NovaSeq6000.{{chr}}.bam", sample=config["samples"])
	output:
		sampleschrBAMmarkdup = expand("results/AlternativeVariantCalling_DNA/{{folder}}/{sample}_{{round}}_merged.{{chr}}.bam", sample=config["samples"])
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
rule r3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2:
	'''
 	Round 1 of base recalibration per chr, sample & platform
 	'''
	input:
		chrVCFHardFiltered = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round1_from_minimap2/ShortVariants_HardFiltered_round0.{chr}.vcf.gz",
		samplechrplatformBAM = "results/AlternativeVariantCalling_DNA/2_merge_3_lanes_bams_markdup_PerChr_from_minimap2/{samplelaneplatform}_round0.{chr}.bam",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		BRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/BaseRecalibration_GATK_{samplelaneplatform}_round1.{chr}.txt",
		samplechrplatformBAMRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/{samplelaneplatform}_round1.{chr}.bam"
	log:
		err = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round1.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round1.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round1_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round1.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BR2{chr}{samplelaneplatform}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample.sh {input.chrVCFHardFiltered} {input.samplechrplatformBAM} {input.genomeFA} {output.BRecal} {output.samplechrplatformBAMRecal} > {log.out} 2> {log.err}"
rule r3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2:
	'''
 	Round 2 of base recalibration per chr, sample & platform
 	'''
	input:
		chrVCFHardFiltered = "results/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_round2_from_minimap2/ShortVariants_HardFiltered_round1.{chr}.vcf.gz",
		samplechrplatformBAM = "results/AlternativeVariantCalling_DNA/2_merge_3_lanes_bams_markdup_PerChr_from_minimap2/{samplelaneplatform}_round1.{chr}.bam",
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		BRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/BaseRecalibration_GATK_{samplelaneplatform}_round2.{chr}.txt",
		samplechrplatformBAMRecal = "results/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/{samplelaneplatform}_round2.{chr}.bam"
	log:
		err = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round2.{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round2.{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample_round2_from_minimap2/3_base_recalibration_GATK_PerChr_PerSample_{samplelaneplatform}_round2.{chr}.txt"
	conda:
		'../envs/VariantCalling_DNA.yaml'
	params:
		time = '02:00:00',
		name = "BR2{chr}{samplelaneplatform}",
		threads = 1,
		mem = 5000
	shell:
		"./scripts/VariantCalling_DNA/3_base_recalibration_GATK_PerChr_PerSample.sh {input.chrVCFHardFiltered} {input.samplechrplatformBAM} {input.genomeFA} {output.BRecal} {output.samplechrplatformBAMRecal} > {log.out} 2> {log.err}"
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

#
rule r4_call_short_variants_freebayes_PerChr_PerSample:
	'''
	Calling short variants for a given chr from a merged bam file for the 3 lanes with freebayes
	'''
	input:
		sampleschrBAMs = expand("results/AlternativeVariantCalling_DNA/r2_merge_3_lanes_bams_markdup_PerChr_from_minimap2/{sample}_merged.{{chr}}.bam", sample=config["samples"]),		
		genomeFA = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		chrVCF = "results/AlternativeVariantCalling_DNA/r4_call_short_variants_freebayes_PerChr_PerSample/{chr}.freebayes.vcf.gz"
	log:
		err = "logs/AlternativeVariantCalling_DNA/r4_call_short_variants_freebayes_PerChr_PerSample/{chr}.err",
		out = "logs/AlternativeVariantCalling_DNA/r4_call_short_variants_freebayes_PerChr_PerSample/{chr}.out"
	benchmark:
		"benchmarks/AlternativeVariantCalling_DNA/r4_call_short_variants_freebayes_PerChr_PerSample/{chr}.txt"
	conda:
		'../envs/AlternativeVariantCalling_DNA.yaml'
	params:
		time = '20:00:00',
		name = "freeb{chr}",
		threads = 1,
		mem = 10000
	shell:
		"./scripts/AlternativeVariantCalling_DNA/4_call_short_variants_freebayes_PerChr_PerSample.sh {wildcards.chr} \"{input.sampleschrBAMs}\" {input.genomeFA} {output.chrVCF} > {log.out} 2> {log.err}"
