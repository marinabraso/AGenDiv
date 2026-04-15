
import re
import random
import gzip


########################################
# Functions
########################################
# Provide bed file name according to the wildcard BED
def get_bed_file_name(wildcards):
	if wildcards.BED == "Callable":
		return rules.r7_join_extra_callable_regions.output.ExtraCallableRegions
	elif wildcards.BED == "All":
		return "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt"
	elif wildcards.BED == "SNPs":
		return rules.VariantType_BEDs_PerSubsetOfSamples.output.SNPs
	elif wildcards.BED == "INDELs":
		return rules.VariantType_BEDs_PerSubsetOfSamples.output.INDELs
	elif wildcards.BED == "Exons":
		return expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Exons")
	elif wildcards.BED == "Introns":
		return expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Introns")
	elif wildcards.BED == "Promoters":
		return expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Promoters")
	elif wildcards.BED == "Intergenic":
		return expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Intergenic")
	elif wildcards.BED == "Singletons":
		return rules.VariantFrequency_BEDs_PerSubsetOfSamples.output.Singletons
	elif wildcards.BED == "Doubletons":
		return rules.VariantFrequency_BEDs_PerSubsetOfSamples.output.Doubletons
	elif wildcards.BED == "Tripletons":
		return rules.VariantFrequency_BEDs_PerSubsetOfSamples.output.Tripletons
	elif wildcards.BED == "BiAllelic":
		return rules.VariantAlleleNumber_BEDs_PerSubsetOfSamples.output.BiAllelic
	elif wildcards.BED == "TriAllelic":
		return rules.VariantAlleleNumber_BEDs_PerSubsetOfSamples.output.TriAllelic
	elif wildcards.BED == "TetraAllelic":
		return rules.VariantAlleleNumber_BEDs_PerSubsetOfSamples.output.TetraAllelic
	elif wildcards.BED == "PentaPlusAllelic":
		return rules.VariantAlleleNumber_BEDs_PerSubsetOfSamples.output.PentaPlusAllelic
	elif wildcards.BED == "WindowsSimChrSize":
		return rules.SimualtedChrSizeWindows_BEDs.output.WindowsBED
	elif wildcards.BED == "SelectedWindowsSimChrSize":
		return rules.SelectRegions_SimualtedChrSize.output.SelectedBED
	elif re.match("AbsFreq", wildcards.BED):
		Afreq = wildcards.BED.split("_")[1]
		BiAll = wildcards.BED.split("_")[2]
		return 	expand(rules.VariantFrequency_BEDs_PerAllSamples.output.VariantList, Afreq=Afreq, BiAll=BiAll)
	elif re.match("FixedLengthRegions", wildcards.BED):
		return 	expand(rules.Windows_AccumulatedCallableSize_BEDs.output.WindowsBED, windowsize=config["windowsize"])
	elif re.match("Simulated", wildcards.BED):
		SimParam="_".join(wildcards.ObsExp.split("/")[2].split("_")[1:])
		Simulationtype=wildcards.ObsExp.split("/")[1]
		if re.search("InfiniteAlleles", wildcards.BED):
			FiniteInfiniteAlleles = "InfiniteAlleles"
		elif re.search("FiniteAlleles", wildcards.BED):
			FiniteInfiniteAlleles = "FiniteAlleles"
		Replica = wildcards.BED.split("_")[1]
		if SimParam != "" and  FiniteInfiniteAlleles != "" and Simulationtype != "" and Replica != "":
			return expand(rules.Formating_SimulationFiles.output.BED, SimFolder="%s/SimParam_%s" % (Simulationtype, SimParam), SimParam="%s_%s" % (SimParam, Replica), FiniteInfiniteAlleles=FiniteInfiniteAlleles)
		else:
			raise ValueError("Unknown value for SimParam, Simulationtype or FiniteInfiniteAlleles get_bed_file_name: %s, %s" % SimParam, FiniteInfiniteAlleles)
	else:
		raise ValueError("Unknown value for BED file in get_bed_file_name: %s" % wildcards.BED)

# Provide list of samples according to wildcard GroupSamples
def get_group_of_samples(wildcards):
	if re.match("Rand", wildcards.GroupSamples): 
		if re.search("AtlSamples", wildcards.GroupSamples):
			numsamples = len(expand("{sample}", sample=config["AtlSamples"]))
			samples = random.sample(expand("{sample}", sample=config["samples"]),numsamples)
			return samples
		elif re.search("MedSamples", wildcards.GroupSamples):
			numsamples = len(expand("{sample}", sample=config["MedSamples"]))
			samples = random.sample(expand("{sample}", sample=config["samples"]),numsamples)
			return samples
		else:
			raise ValueError("Unknown value for GroupSamples in get_group_of_samples: %s" % wildcards.GroupSamples)
	elif wildcards.GroupSamples == "AllSamples":
		return expand("{sample}", sample=config["samples"])
	elif wildcards.GroupSamples == "AtlSamples":
		return expand("{sample}", sample=config["AtlSamples"])
	elif wildcards.GroupSamples == "MedSamples":
		return expand("{sample}", sample=config["MedSamples"])
	elif wildcards.GroupSamples == "SimSamples":
		return expand("{sample}", sample=range(1,config["SimulationsSampleSize"]+1))
	else:
		raise ValueError("Unknown value for GroupSamples in get_group_of_samples: %s" % wildcards.GroupSamples)


# Provide heterozygosity files according to the wildcard ObsExp
def get_HetFiles(wildcards):
	if wildcards.ObsExp == "Observed_Data":
		return expand(rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile, chr=config["bralan3chrs"])
	elif wildcards.ObsExp == "Observed_SNPs":
		return expand(rules.VariantType_BEDs_PerSubsetOfSamples.output.SNPsHetFile, GroupSamples="AllSamples")
	elif wildcards.ObsExp == "Observed_INDELs":
		return expand(rules.VariantType_BEDs_PerSubsetOfSamples.output.INDELsHetFile, GroupSamples="AllSamples")
	elif re.match("Expected_AsInSimulations", wildcards.ObsExp):
		SimParam="_".join(wildcards.ObsExp.split("/")[2].split("_")[1:])
		Simulationtype=wildcards.ObsExp.split("/")[1]
		if re.search("InfiniteAlleles", wildcards.BED):
			FiniteInfiniteAlleles = "InfiniteAlleles"
		elif re.search("FiniteAlleles", wildcards.BED):
			FiniteInfiniteAlleles = "FiniteAlleles"
		Replica = wildcards.BED.split("_")[1]
		if SimParam != "" and  FiniteInfiniteAlleles != "" and Simulationtype != "" and Replica != "":
			return expand(rules.Formating_SimulationFiles.output.HetFile, SimFolder="%s/SimParam_%s" % (Simulationtype, SimParam), SimParam="%s_%s" % (SimParam, Replica), FiniteInfiniteAlleles=FiniteInfiniteAlleles)
		else:
			raise ValueError("Unknown value for SimParam, Simulationtype or FiniteInfiniteAlleles get_HetFiles: %s, %s" % SimParam, FiniteInfiniteAlleles)
	else:
		raise ValueError("Unknown value for ObsExp file in get_HetFiles: %s" % wildcards.ObsExp)

# Provide genotype files according to the wildcard ObsExp
def get_GenoFiles(wildcards):
	if wildcards.ObsExp == "Observed_Data":
		return expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE, chr=config["bralan3chrs"])
	elif wildcards.ObsExp == "Observed_SNPs":
		return expand(rules.VariantType_BEDs_PerSubsetOfSamples.output.SNPsGenotype, GroupSamples=wildcards.GroupSamples, BED=wildcards.BED)
	elif wildcards.ObsExp == "Observed_INDELs":
		return expand(rules.VariantType_BEDs_PerSubsetOfSamples.output.INDELsGenotype, GroupSamples=wildcards.GroupSamples, BED=wildcards.BED)
	elif re.match("Expected_AsInSimulations", wildcards.ObsExp):
		SimParam="_".join(wildcards.ObsExp.split("/")[2].split("_")[1:])
		Simulationtype=wildcards.ObsExp.split("/")[1]
		if re.search("InfiniteAlleles", wildcards.BED):
			FiniteInfiniteAlleles = "InfiniteAlleles"
		elif re.search("FiniteAlleles", wildcards.BED):
			FiniteInfiniteAlleles = "FiniteAlleles"
		Replica = wildcards.BED.split("_")[1]
		if SimParam != "" and  FiniteInfiniteAlleles != "" and Simulationtype != "" and Replica != "":
			return expand(rules.Formating_SimulationFiles.output.GenoFile, SimFolder="%s/SimParam_%s" % (Simulationtype, SimParam), SimParam="%s_%s" % (SimParam, Replica), FiniteInfiniteAlleles=FiniteInfiniteAlleles)
		else:
			raise ValueError("Unknown value for SimParam, Simulationtype or FiniteInfiniteAlleles get_GenoFiles: %s, %s" % SimParam, FiniteInfiniteAlleles)
	else:
		raise ValueError("Unknown value for ObsExp file in get_GenoFiles: %s" % wildcards.ObsExp)

# Provide file with list of samples according to wildcard GroupSamples
def get_group_of_samples_file(wildcards):
	if re.match("Rand", wildcards.GroupSamples): 
		if re.search("AtlSamples", wildcards.GroupSamples):
			b=wildcards.GroupSamples.split("_")[1]
			return expand(rules.Prepare_PopulationSamples_Files.output.outAtl, ObsExp=wildcards.ObsExp, RandObsSim="Rand", boots=b)
		elif re.search("MedSamples", wildcards.GroupSamples):
			b=wildcards.GroupSamples.split("_")[1]
			return expand(rules.Prepare_PopulationSamples_Files.output.outMed, ObsExp=wildcards.ObsExp, RandObsSim="Rand", boots=b)
		else:
			raise ValueError("Unknown value for GroupSamples in get_group_of_samples_file: %s" % wildcards.GroupSamples)
	elif wildcards.GroupSamples == "AllSamples":
		return expand(rules.Prepare_PopulationSamples_Files.output.outAll, ObsExp=wildcards.ObsExp, RandObsSim="Obs", boots=0)
	elif wildcards.GroupSamples == "AtlSamples":
		return expand(rules.Prepare_PopulationSamples_Files.output.outAtl, ObsExp=wildcards.ObsExp, RandObsSim="Obs", boots=0)
	elif wildcards.GroupSamples == "MedSamples":
		return expand(rules.Prepare_PopulationSamples_Files.output.outMed, ObsExp=wildcards.ObsExp, RandObsSim="Obs", boots=0)
	elif wildcards.GroupSamples == "SimSamples":
		return expand(rules.Prepare_PopulationSamples_Files.output.outAll, ObsExp=wildcards.ObsExp, RandObsSim="Sim", boots=0)
	else:
		raise ValueError("Unknown value for GroupSamples in get_group_of_samples_file: %s" % wildcards.GroupSamples)
# Provide file with Frequency per site depending on BED, GroupSamples and ObsExp
def get_freq_file_Pops_randPops(wildcards):
	if wildcards.ObsOrBoots == "Observed":
		return expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples=["AllSamples","AtlSamples","MedSamples"], ObsExp="Observed_Data")
	elif re.match("Rand", wildcards.ObsOrBoots):
		b=wildcards.ObsOrBoots.split("_")[1]
		return expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="AllSamples", ObsExp="Observed_Data") + expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples=["RandAtlSamples_%s" % b,"RandMedSamples_%s" % b], ObsExp="Observed_Data")
	else:
		raise ValueError("Unknown value for ObsOrBoots in get_freq_file_Pops_randPops: %s" % wildcards.ObsOrBoots)

# Provide 2 BED file names for intersection 
def get_bed_files_names_for_intersection(wildcards):
	if wildcards.Name1 == "Callable" and wildcards.Name2 == "Exons":
		return rules.r7_join_extra_callable_regions.output.ExtraCallableRegions, rules.FunctionalFeatures_BEDs.output.Exons
	if wildcards.Name1 == "Callable" and  wildcards.Name2 == "Introns":
		return rules.r7_join_extra_callable_regions.output.ExtraCallableRegions, rules.FunctionalFeatures_BEDs.output.Introns
	if wildcards.Name1 == "Callable" and  wildcards.Name2 == "Promoters":
		return rules.r7_join_extra_callable_regions.output.ExtraCallableRegions, rules.FunctionalFeatures_BEDs.output.Promoters
	if wildcards.Name1 == "Callable" and  wildcards.Name2 == "Intergenic":
		return rules.r7_join_extra_callable_regions.output.ExtraCallableRegions, rules.FunctionalFeatures_BEDs.output.Intergenic
	else:
		raise ValueError("Unknown values for BED files in get_bed_files_names_for_intersection: %s & %s" % (wildcards.Name1, wildcards.Name2))
#




########################################
#### Rules preparing BED files
########################################
# Compute the frequencey of each allele in each posistion for a subset of BED regions and group of samples 
rule FrequencyPerSite_inBEDregions_PerSubsetOfSamples:
	'''
	Calculate the frequencey of each allele in each posistion for a subset of regions (bed file) and a subset of samples (AllSamples, AtlSamples or MedSamples)
	'''
	input:
		GenoFiles = get_GenoFiles,
		BED = get_bed_file_name,
		SamplesFile = get_group_of_samples_file,
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		out = "results/Plotting_DNA/{ObsExp}/FrequencyPerSite_inBEDregions_PerSubsetOfSamples/FrequencyPerSite_in{BED}Regions_Per{GroupSamples}.txt.gz"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/FrequencyPerSite_inBEDregions_PerSubsetOfSamples/FrequencyPerSite_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/FrequencyPerSite_inBEDregions_PerSubsetOfSamples/FrequencyPerSite_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/FrequencyPerSite_inBEDregions_PerSubsetOfSamples/FrequencyPerSite_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "freq_{BED}_{GroupSamples}",
		threads = 1,
		mem = 50000
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output%.gz}}) 2> {log.err}
		if [[ {wildcards.ObsExp} =~ "Expected_AsInSimulations" ]]; then
			samplelines=$(cat {input.SamplesFile} | awk '{{str=str","$1+5}}END{{print str}}' | sed 's/^,//g' ) 2> {log.err}
		else
			samplelines=$(awk '{{if(NR==FNR){{a[$1]=1;next}}if(a[$1]){{str=str","FNR+4}}}}END{{print str}}' <(cat {input.SamplesFile} | sed 's/\\s/\\n/g') <(cat {input.SamplesOrderInVCF} | sed 's/\\s/\\n/g') | sed 's/^,//g') 2> {log.err}
		fi
		bedtools intersect -a <(for i in {input.GenoFiles}; do zcat < $i; done) -b <(zcat < {input.BED}) | cut -f${{samplelines}} | sed 's/:/\t/g' | perl -ne '@a = split(" ", $_); %h=(); $h{{$_}}++ for @a; $str=""; for(sort {{ $h{{$a}} <=> $h{{$b}} }} keys(%h)){{$str=$str."$_:$h{{$_}},"}} print $str."\n"' | sed 's/,$//g' 2> ~/.null > ${{output%.gz}}.tmp
		paste <(bedtools intersect -a <(for i in {input.GenoFiles}; do zcat < $i; done) -b <(zcat < {input.BED}) | cut -f1,2,3) ${{output%.gz}}.tmp | awk '{{n=split($4,a,","); if(n>1){{print $0}}}}' | gzip > ${{output}} 2> {log.err}
		rm ${{output%.gz}}.tmp
		"""

# to test 
rule SynonymousNonSynonymous_BEDs:
	'''
	Extract which sites are synonymous and which are nonsynonymous from GTF 
	'''
	input:
		DSN_BED=rules.Synonyms_and_nonsynonymous_fromGTF.output.SynNonSynBED_callable,
		Freqfile = expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="AllSamples", ObsExp="Observed_Data")
	output:
		Syn = "results/Plotting_DNA/SynonymousNonSynonymous_BEDs/Syn.bed.gz",
		NonSyn = "results/Plotting_DNA/SynonymousNonSynonymous_BEDs/NonSyn.bed.gz",
		AllVar_SynNonSyn = "results/Plotting_DNA/SynonymousNonSynonymous_BEDs/AllVar_SynNonSyn.bed.gz"
	log:
		err = "logs/Plotting_DNA/SynonymousNonSynonymous_BEDs/SynonymousNonSynonymous_BEDs.err",
		out = "logs/Plotting_DNA/SynonymousNonSynonymous_BEDs/SynonymousNonSynonymous_BEDs.out"
	benchmark:
		"benchmarks/Plotting_DNA/SynonymousNonSynonymous_BEDs/SynonymousNonSynonymous_BEDs.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '5:00:00',
		name = "SnSbed",
		threads = 1,
		mem = 50000,
		SnSfromSeq = "scripts/Plotting_DNA/Synonyms_and_nonsynonymous_fromSeq.pl",
		ProtfromSeq = "scripts/Plotting_DNA/ProteinSeq_fromTranscriptSeq.pl"
	shell:
		"""
		cat {input.DSN_BED} | grep 'S' | gzip > {output.Syn}
		cat {input.DSN_BED} | grep 'N' | gzip > {output.NonSyn}
		awk '{{if(NR==FNR){{a[$1"\\t"$2"\\t"$3]=$4; next}} if(a[$1"\t"$2"\t"$3]){{print $1"\\t"$2"\\t"$3"\\t"a[$1"\\t"$2"\\t"$3]}}else{{print $0"\tNA"}}}}' {input.DSN_BED} <(zcat {input.Freqfile} | cut -f1,2,3) | gzip > {output.AllVar_SynNonSyn}
		"""



# Create bed files for Exons, Introns, Promoters & Intergenic regions
rule FunctionalFeatures_BEDs:
	'''
	'''
	input:
		ChrLengths="data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt",
		GTF = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_strong.gtf.gz",
		Freqfile = expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="AllSamples", ObsExp="Observed_Data")
	output:
		Exons = "results/Plotting_DNA/FunctionalFeatures_BEDs/Exons.bed.gz",
		Introns = "results/Plotting_DNA/FunctionalFeatures_BEDs/Introns.bed.gz",
		Promoters = "results/Plotting_DNA/FunctionalFeatures_BEDs/Promoters.bed.gz",
		Intergenic = "results/Plotting_DNA/FunctionalFeatures_BEDs/Intergenic.bed.gz",
		PerSite_FeatureType = "results/Plotting_DNA/FunctionalFeatures_BEDs/PerSite_VariantType.txt.gz"
	log:
		err = "logs/Plotting_DNA/FunctionalFeatures_BEDs/FunctionalFeatures_BEDs.err",
		out = "logs/Plotting_DNA/FunctionalFeatures_BEDs/FunctionalFeatures_BEDs.out"
	benchmark:
		"benchmarks/Plotting_DNA/FunctionalFeatures_BEDs/FunctionalFeatures_BEDs.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "FuncBED",
		threads = 1,
		mem = 10000
	shell:
		"scripts/Plotting_DNA/FunctionalFeatures_BEDs.sh {input.ChrLengths} {input.GTF} {input.Freqfile} {output.Exons} {output.Introns} {output.Promoters} {output.Intergenic} {output.PerSite_FeatureType} > {log.out} 2> {log.err}"

rule SimualtedChrSizeWindows_BEDs:
	'''
	'''
	input:
		CallableBED = rules.r7_join_extra_callable_regions.output.ExtraCallableRegions
	output:
		WindowsBED = "results/Plotting_DNA/SimualtedChrSizeWindows_BEDs/SimualtedChrSizeWindows.bed.gz"
	log:
		err = "logs/Plotting_DNA/SimualtedChrSizeWindows_BEDs/SimualtedChrSizeWindows_BEDs.err",
		out = "logs/Plotting_DNA/SimualtedChrSizeWindows_BEDs/SimualtedChrSizeWindows_BEDs.out"
	benchmark:
		"benchmarks/Plotting_DNA/SimualtedChrSizeWindows_BEDs/SimualtedChrSizeWindows_BEDs.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "WindBED",
		threads = 1,
		mem = 10000,
		simChsize=config["SimulationsChrSize"]
	shell:
		"""
		zcat {input.CallableBED} | awk -v chrsize={params.simChsize} '{{if($4 >= chrsize){{st=$2; while(st+chrsize<=$3){{print $1"\t"st"\t"st+chrsize; st+=chrsize}}}}}}' | gzip > {output.WindowsBED} 2> {log.err}
		"""


# to check
rule VariantType_BEDs_PerSubsetOfSamples:
	'''
	'''
	input:
		GenoFiles = expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE, chr=config["bralan3chrs"]),
		HetFiles = expand(rules.Prepare_basic_analysis_files_PerChr.output.HeterozygosityFile, chr=config["bralan3chrs"]),
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		SNPs = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/SNPs_Per{GroupSamples}.bed.gz",
		INDELs = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/INDELs_Per{GroupSamples}.bed.gz",
		PerSite_VariantType = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/PerSite_VariantType_Per{GroupSamples}.txt.gz",
		SNPsGenotype = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/SNPs_Per{GroupSamples}.genotype.bed.gz",
		INDELsGenotype = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/INDELs_Per{GroupSamples}.genotype.bed.gz",
		SNPsHetFile = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/SNPs_Per{GroupSamples}.0hom1het.bed.gz",
		INDELsHetFile = "results/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/INDELs_Per{GroupSamples}.0hom1het.bed.gz",
	log:
		err = "logs/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/VariantType_BEDs_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/VariantType_BEDs_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/VariantType_BEDs_PerSubsetOfSamples/VariantType_BEDs_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "Vtype_{GroupSamples}",
		threads = 1,
		mem = 10000,
		samples = get_group_of_samples
	shell:
		"""
		mkdir -p $(dirname {output.SNPs}) 2> {log.err}
		samplelines=$(awk '{{if(NR==FNR){{a[$1]=1;next}}if(a[$1]){{str=str","FNR+4}}}}END{{print str}}' <(echo {params.samples} | sed 's/\\s/\\n/g') <(cat {input.SamplesOrderInVCF} | sed 's/\\s/\\n/g') | sed 's/^,//g') 2> {log.err}
		for i in {input.GenoFiles}; do zcat < $i; done | cut -f1,2,3,4,${{samplelines}}  | awk '{{for(i=5;i<=NF;i++){{split($i,a,":");for(j in a){{al[a[j]]=1}}}} min=999999999; max=0; for(i in al){{n=split(i,a,"");if(n>max){{max=n}}if(n<min){{min=n}}if(i=="*" || i=="."){{min=0}}}}; if(min==1 && max==1){{type="SNP"}}else{{type="INDEL"}} print $1"\t"$2"\t"$3"\t"type; delete al}}' | gzip > {output.PerSite_VariantType} 2> {log.err}
		zcat < {output.PerSite_VariantType} | awk '{{if($4 ~/SNP/){{print $1"\\t"$2"\\t"$3}}}}' | gzip > {output.SNPs} 2> {log.err}
		zcat < {output.PerSite_VariantType} | awk '{{if($4 ~/INDEL/){{print $1"\\t"$2"\\t"$3}}}}' | gzip > {output.INDELs} 2> {log.err}
		awk '{{if(NR==FNR){{a[$1"_"$2"_"$3]=1; next}}if(a[$1"_"$2"_"$3]==1){{print $0}}}}' <(zcat {output.SNPs}) <(for i in {input.GenoFiles}; do zcat < $i; done) | gzip > {output.SNPsGenotype} 2> {log.err}
		awk '{{if(NR==FNR){{a[$1"_"$2"_"$3]=1; next}}if(a[$1"_"$2"_"$3]==1){{print $0}}}}' <(zcat {output.INDELs}) <(for i in {input.GenoFiles}; do zcat < $i; done) | gzip > {output.INDELsGenotype} 2> {log.err}
		awk '{{if(NR==FNR){{a[$1"_"$2"_"$3]=1; next}}if(a[$1"_"$2"_"$3]==1){{print $0}}}}' <(zcat {output.SNPs}) <(for i in {input.HetFiles}; do zcat < $i; done) | gzip > {output.SNPsHetFile} 2> {log.err}
		awk '{{if(NR==FNR){{a[$1"_"$2"_"$3]=1; next}}if(a[$1"_"$2"_"$3]==1){{print $0}}}}' <(zcat {output.INDELs}) <(for i in {input.HetFiles}; do zcat < $i; done) | gzip > {output.INDELsHetFile} 2> {log.err}
		"""

rule SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples:
	'''
	'''
	input:
		Freqfile = expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="{GroupSamples}", ObsExp="Observed_Data")
	output:
		Singletons = "results/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/Singletons_Per{GroupSamples}.bed.gz",
		DoubletonsBial = "results/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/Doubletons_Biallelic_Per{GroupSamples}.bed.gz",
		DoubletonsMultial = "results/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/Doubletons_Multiallelic_Per{GroupSamples}.bed.gz",
		TripletonsBial = "results/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/Tripletons_Biallelic_Per{GroupSamples}.bed.gz",
		TripletonsMultial = "results/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/Tripletons_Multiallelic_Per{GroupSamples}.bed.gz",
	log:
		err = "logs/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/VariantFrequency_BEDs_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/VariantFrequency_BEDs_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/SingletonsDoubletonsTripletons_BEDs_PerSubsetOfSamples/VariantFrequency_BEDs_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "Singl_{GroupSamples}",
		threads = 1,
		mem = 10000,
		samples = get_group_of_samples
	shell:
		"""
		mkdir -p $(dirname {output.Singletons}) 2> {log.err}

		## Singletons
		nsamp=$(echo "{params.samples}" | awk '{{print NF*2-1}}')
		echo ${{nsamp}}
		zcat < {input.Freqfile} | awk '{{if(split($4,a,",")==2){{print $0}}}}' | sed 's/\\t\\S*\\:\\([0-9]*\\)$/\\t\\1/g' | grep -P "\\t${{nsamp}}$" | cut -f1,2,3 | gzip > {output.Singletons} 2> {log.err}


		## Doubletons
		nsamp=$(echo "{params.samples}" | awk '{{print NF*2-2}}')
		echo ${{nsamp}}
		zcat < {input.Freqfile} | awk '{{if(split($4,a,",")==2){{print $0}}}}' | sed 's/\\t\\S*\\:\\([0-9]*\\)$/\\t\\1/g' | grep -P "\\t${{nsamp}}$" | cut -f1,2,3 | gzip > {output.DoubletonsBial} 2> {log.err}
		zcat < {input.Freqfile} | awk '{{if(split($4,a,",")>2){{print $0}}}}' | sed 's/\\t\\S*\\:\\([0-9]*\\)$/\\t\\1/g' | grep -P "\\t${{nsamp}}$" | cut -f1,2,3 | gzip > {output.DoubletonsMultial} 2> {log.err}


		## Tripletons
		nsamp=$(echo "{params.samples}" | awk '{{print NF*2-3}}')
		echo ${{nsamp}}
		zcat < {input.Freqfile} | awk '{{if(split($4,a,",")==2){{print $0}}}}' | sed 's/\\t\\S*\\:\\([0-9]*\\)$/\\t\\1/g' | grep -P "\\t${{nsamp}}$" | cut -f1,2,3 | gzip > {output.TripletonsBial} 2> {log.err}
		zcat < {input.Freqfile} | awk '{{if(split($4,a,",")>2){{print $0}}}}' | sed 's/\\t\\S*\\:\\([0-9]*\\)$/\\t\\1/g' | grep -P "\\t${{nsamp}}$" | cut -f1,2,3 | gzip > {output.TripletonsMultial} 2> {log.err}
		"""
# Divide variants per frequency
rule VariantFrequency_BEDs_PerAllSamples:
	'''
	'''
	input:
		Freqfile = expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="AllSamples", ObsExp="Observed_Data")
	output:
		VariantList = "results/Plotting_DNA/VariantFrequency_BEDs_PerAllSamples/Variants_AbsFreq_{Afreq}_{BiAll}.bed.gz",
	log:
		err = "logs/Plotting_DNA/VariantFrequency_BEDs_PerAllSamples/Variants_AbsFreq_{Afreq}_{BiAll}.err",
		out = "logs/Plotting_DNA/VariantFrequency_BEDs_PerAllSamples/Variants_AbsFreq_{Afreq}_{BiAll}.out"
	benchmark:
		"benchmarks/Plotting_DNA/VariantFrequency_BEDs_PerAllSamples/Variants_AbsFreq_{Afreq}_{BiAll}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "Vfreq_{Afreq}_{BiAll}",
		threads = 1,
		mem = 10000,
	shell:
		"""
		mkdir -p $(dirname {output.VariantList}) 2> {log.err}
		stFreq=$(echo {wildcards.Afreq} | awk '{{split($1,a,"-"); print a[1]}}')
		endFreq=$(echo {wildcards.Afreq} | awk  -v stFreq=${{stFreq}} '{{split($1,a,"-"); if(a[2]){{print a[2]}}else{{print stFreq+1}}}}')
		echo ${{stFreq}} ${{endFreq}}
		if [ {wildcards.BiAll} == "Biallelic" ]
		then
			zcat < {input.Freqfile} | awk -v stFreq=${{stFreq}} -v endFreq=${{endFreq}} '{{n=split($4,a,","); split(a[length(a)],b,":"); if(n==2 && b[2]>=stFreq && b[2]<endFreq){{print $0}}}}' | cut -f1,2,3 | gzip > {output.VariantList} 2> {log.err}
		else
			zcat < {input.Freqfile} | awk -v stFreq=${{stFreq}} -v endFreq=${{endFreq}} '{{n=split($4,a,","); split(a[length(a)],b,":"); if(b[2]>=stFreq && b[2]<endFreq){{print $0}}}}' | cut -f1,2,3 | gzip > {output.VariantList} 2> {log.err}
		fi
		"""

rule VariantAlleleNumber_BEDs_PerSubsetOfSamples:
	'''
	'''
	input:
		Freqfile = expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="{GroupSamples}", ObsExp="Observed_Data")
	output:
		BiAllelic = "results/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/BiAllelic_Per{GroupSamples}.bed.gz",
		TriAllelic = "results/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/TriAllelic_Per{GroupSamples}.bed.gz",
		TetraAllelic = "results/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/TetraAllelic_Per{GroupSamples}.bed.gz",
		PentaPlusAllelic = "results/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/PentaPlusAllelic_Per{GroupSamples}.bed.gz",
	log:
		err = "logs/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/VariantAlleleNumber_BEDs_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/VariantAlleleNumber_BEDs_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/VariantAlleleNumber_BEDs_PerSubsetOfSamples/VariantAlleleNumber_BEDs_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "VAllNum_{GroupSamples}",
		threads = 1,
		mem = 10000
	shell:
		"""
		mkdir -p $(dirname {output.BiAllelic}) 2> {log.err}
		zcat < {input.Freqfile} | awk '{{n=split($4,a,","); if(n==2){{print $0}}}}' | gzip > {output.BiAllelic}
		zcat < {input.Freqfile} | awk '{{n=split($4,a,","); if(n==3){{print $0}}}}' | gzip > {output.TriAllelic}
		zcat < {input.Freqfile} | awk '{{n=split($4,a,","); if(n==4){{print $0}}}}' | gzip > {output.TetraAllelic}
		zcat < {input.Freqfile} | awk '{{n=split($4,a,","); if(n>4){{print $0}}}}' | gzip > {output.PentaPlusAllelic}
		"""

# Determine the proportion of Private variants, Shared...
rule SharedPrivatePopulation_BEDs:
	'''
	'''
	input:
		Freqfiles = get_freq_file_Pops_randPops
	output:
		Numbers = "results/Plotting_DNA/SharedPrivatePopulation_BEDs/Numbers_{ObsOrBoots}.txt"		
	log:
		err = "logs/Plotting_DNA/SharedPrivatePopulation_BEDs/SharedPrivatePopulation_BEDs_{ObsOrBoots}.err",
		out = "logs/Plotting_DNA/SharedPrivatePopulation_BEDs/SharedPrivatePopulation_BEDs_{ObsOrBoots}.out"
	benchmark:
		"benchmarks/Plotting_DNA/SharedPrivatePopulation_BEDs/SharedPrivatePopulation_BEDs_{ObsOrBoots}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "VShPriv_{ObsOrBoots}",
		threads = 1,
		mem = 100000
	shell:
		"""
		mkdir -p $(dirname {output.Numbers}) 2> {log.err}

		inputs=( $(echo {input.Freqfiles} | sed 's/ /\\n/g') )
		inputAll=${{inputs[0]}}
		inputAtl=${{inputs[1]}}
		inputMed=${{inputs[2]}}
		awk '{{if(FNR==1){{file++}}if(file==1){{a[$1"\\t"$2"\\t"$3]=$4; next}}if(file==2){{m[$1"\\t"$2"\\t"$3]=$4;next}} print $0"\\t"a[$1"\\t"$2"\\t"$3]"\\t"m[$1"\\t"$2"\\t"$3]}}' <(zcat ${{inputAtl}}) <(zcat ${{inputMed}}) <(zcat ${{inputAll}}) | gzip > $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp.gz
		zcat $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp.gz | sed 's/:[0-9]\\+//g' | sed 's/\\t\\t/\\t.\\t/g' | sed 's/\\t$/\\t./g' | awk '{{a=$5; m=$6; if(a=="." && m!="."){{print $0"\\tPrivateM"}}else if(a!="." && m=="."){{print $0"\\tPrivateA"}}else if(a=="." && m=="."){{print $0"\\tFixedDifferent"}}else{{split(a,aAl,",");split(m,mAl,",");commAl=0; for(i in aAl){{for(j in mAl){{if(mAl[j]==aAl[i]){{commAl++}}}}}}if(commAl<=1){{print $0"\\tVariantDifferent"}}else{{print $0"\\tVariantShared"}}}}}}' | gzip > $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp2.gz
		paste <(zcat $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp.gz | cut -f-6) <(zcat $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp2.gz | cut -f7) | sed 's/\\t\\t/\\t.\\t/g' | sed 's/\\t\\t/\\t.\\t/g' | gzip > $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp3.gz
		rm $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp.gz $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp2.gz
		zcat $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp3.gz | cut -f7 | sort | uniq -c | sed 's/\\s\\+/\\t/g' | sed 's/^\\t//g' > {output.Numbers}
		rm $(dirname {output.Numbers})/AllVariants_BothPopulations_{wildcards.ObsOrBoots}.tmp3.gz
		"""

rule Intersection2BEDs:
	'''
	'''
	input:
		BEDs = get_bed_files_names_for_intersection
	output:
		Intersection = "results/Plotting_DNA/Intersection2BEDs/{Name1}_{Name2}.bed.gz",
	log:
		err = "logs/Plotting_DNA/Intersection2BEDs/Intersection2BEDs_{Name1}_{Name2}.err",
		out = "logs/Plotting_DNA/Intersection2BEDs/Intersection2BEDs_{Name1}_{Name2}.out"
	benchmark:
		"benchmarks/Plotting_DNA/Intersection2BEDs/Intersection2BEDs_{Name1}_{Name2}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "i_{Name1}_{Name2}",
		threads = 1,
		mem = 10000
	shell:
		"""
		mkdir -p $(dirname {output.Intersection}) 2> {log.err}
		BED1=$(echo {input.BEDs} | cut -f1 -d' ')
		BED2=$(echo {input.BEDs} | cut -f2 -d' ')
		bedtools intersect -a <(zcat ${{BED1}} | cut -f1,2,3) -b <(zcat ${{BED2}} | cut -f1,2,3) | gzip > {output.Intersection}
		"""

# add group of samples!!!! we need co compute the two pop separately
rule Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions:
	'''
	'''
	input:
		ChrGenotypesFile=get_GenoFiles,
		BED = get_bed_file_name,
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt"
	output:
		OutMajorFreq = "results/Plotting_DNA/{ObsExp}/in{BED}Regions_Per{GroupSamples}.MajorFreq.gz",
		OutNumAlleles = "results/Plotting_DNA/{ObsExp}/in{BED}Regions_Per{GroupSamples}.NumAlleles.gz"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Calc_NumAlleles_MajorFreq_InEachRegion_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/Calc_NumAlleles_MajorFreq_InEachRegion_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Calc_NumAlleles_MajorFreq_InEachRegion_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "MajFNumAll_{ObsExp}_{BED}_{GroupSamples}",
		threads = 1,
		mem = 100000,
		samples = get_group_of_samples,
	shell:
		"""
		outMF={output.OutMajorFreq}
		if [[ {wildcards.GroupSamples} == "SimSamples" ]]
		then
			samplelines=$(awk '{{str=str","$1+4}}END{{print str}}' <(echo {params.samples} | sed 's/\\s/\\n/g') | sed 's/^,//g')
		else
			samplelines=$(awk '{{if(NR==FNR){{a[$1]=1;next}}if(a[$1]){{str=str","FNR+4}}}}END{{print str}}' <(echo {params.samples} | sed 's/\\s/\\n/g') <(cat {input.SamplesOrderInVCF} | sed 's/\\s/\\n/g') | sed 's/^,//g')
		fi
		echo ${{samplelines}} > {log.out}
 		bedtools intersect -a <(zcat {input.BED}) -b <(zcat {input.ChrGenotypesFile}) -wao | cut -f1,2,3,7- | rev | cut -f2- | rev | cut -f1,2,3,4,${{samplelines}} | gzip > ${{outMF%.MajorFreq.gz}}.bed.gz 2> {log.err};
		scripts/Plotting_DNA/Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.py ${{outMF%.MajorFreq.gz}}.bed.gz {output.OutMajorFreq} {output.OutNumAlleles} >> {log.out} 2>> {log.err};
		"""
# Build a BED file with windows with exactely windowsize callable length
rule Windows_AccumulatedCallableSize_BEDs:
	'''
	'''
	input:
		CallableBED = rules.r7_join_extra_callable_regions.output.ExtraCallableRegions,
	output:
		WindowsBED = "results/Plotting_DNA/Windows_AccumulatedCallableSize_BEDs/Windows_AccumulatedCallableSize_{windowsize}.bed.gz"
	log:
		err = "logs/Plotting_DNA/Windows_AccumulatedCallableSize_BEDs/Windows_AccumulatedCallableSize_{windowsize}.err",
		out = "logs/Plotting_DNA/Windows_AccumulatedCallableSize_BEDs/Windows_AccumulatedCallableSize_{windowsize}.out"
	benchmark:
		"benchmarks/Plotting_DNA/Windows_AccumulatedCallableSize_BEDs/Windows_AccumulatedCallableSize_{windowsize}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "WindBED",
		threads = 1,
		mem = 10000
	shell:
		"""
		output=$(echo {output.WindowsBED}) 2> {log.err}
		if [[ -s ${{output}} || -s ${{output%.gz}} ]]
		then
			rm ${{output%.gz}}*
		fi
		chrs=$(zcat {input.CallableBED} | cut -f1 | sort | uniq)
		for chr in ${{chrs}}
		do
			echo ${{chr}}
			zcat {input.CallableBED} | grep -w ${{chr}} | awk -v wlen={wildcards.windowsize} '{{for(i=$2;i<=$3;i++){{a[i]=1}}}}END{{st=1; for(i=1;i<=$3;i++){{if(a[i]){{accl++}}if(accl==wlen){{print $1"\t"st"\t"i; st=i+1; accl=0}}}}}}' >> ${{output%.gz}}
		done
		gzip ${{output%.gz}}
		"""
#
rule SelectRegions_SimualtedChrSize:
	'''
	'''
	input:
		BigBED = rules.SimualtedChrSizeWindows_BEDs.output.WindowsBED
	output:
		SelectedBED = "results/Plotting_DNA/CompareTreeShape/SelectRegions_SimualtedChrSize/Selected_SimChromSizeRegions.bed.gz"
	log:
		err = "logs/Plotting_DNA/CompareTreeShape/SelectRegions_SimualtedChrSize.err",
		out = "logs/Plotting_DNA/CompareTreeShape/SelectRegions_SimualtedChrSize.out"
	benchmark:
		"benchmarks/Plotting_DNA/CompareTreeShape/SelectRegions_SimualtedChrSize.txt"
	conda:
		'../envs/Phylogenetics.yaml'
	params:
		time = '02:00:00',
		name = "SelectBED",
		threads = 1,
		mem = 2000,
		TimesRandomTest = config["TimesRandomTest"]
	shell:
		"""
		outfile={output.SelectedBED}
		mkdir -p $(dirname ${{outfile}})
		if [[ -s ${{outfile}} ]]
		then
			rm ${{outfile}} 2> ~/null
		fi
		for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; do zcat {input.BigBED} | grep -w "chr"$c | shuf | head -{params.TimesRandomTest}; done | gzip -c > ${{outfile}}
		"""

########################################
# Rules computing statistics
########################################
# Compute total heterozygosity in a set of BED regions for each sample
rule Heterozygosity_Total_inBEDregions_PerSample:
	'''
	Calculate overall heterosygostity for each sample in a set of regions (bed file)
	Output examples:
		- results/Plotting_DNA/Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_7.5e-7_10000/Heterozygosity_Total_inSimulated_1_InfiniteAllelesRegions_PerSample.txt
	'''
	input:
		HetFiles = get_HetFiles,
		BED = get_bed_file_name
	output:
		out = "results/Plotting_DNA/{ObsExp}/Heterozygosity_Total_in{BED}Regions_PerSample.txt"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Heterozygosity_Total_in{BED}Regions_PerSample.err",
		out = "logs/Plotting_DNA/{ObsExp}/Heterozygosity_Total_in{BED}Regions_PerSample.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Heterozygosity_Total_in{BED}Regions_PerSample.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "het{BED}",
		threads = 1,
		mem = 10000
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output}}) 2> {log.err}
		awk '{{if(NR==FNR){{Tlen+=$3-$2;next}} if(FNR==1){{nsamples=NF}}for(i=1; i<=NF; i++){{a[i]+=$i}}}}END{{str=""; for(i=1; i<=nsamples;i++){{str=str"\t"a[i]/Tlen}} print str}}' <(zcat {input.BED}) <(bedtools intersect -a <(zcat {input.HetFiles}) -b <(zcat {input.BED}) | cut -f5-) | sed 's/^\t//g' > ${{output}} 2>> {log.err}
		"""

# Compute heterozygosity for each one of a set of BED regions for each sample
# to be tested
# Fix: take regions into account as a whole, not divided in its callable regions
rule Heterozygosity_InEachRegion_inBEDregions_PerSample:
	'''
	Calculate heterosygostity for each sample and region (bed file)
	Output examples:
		- results/Plotting_DNA/Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_7.5e-7_10000/HeterozygosityPerSample_inSimulated_1_InfiniteAllelesRegions.txt
	'''
	input:
		HetFiles = get_HetFiles,
		BED = get_bed_file_name
	output:
		out = "results/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_in{BED}Regions_PerSample.bed.gz"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_in{BED}Regions_PerSample.err",
		out = "logs/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_in{BED}Regions_PerSample.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_in{BED}Regions_PerSample.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "het{BED}",
		threads = 1,
		mem = 100000
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output}}) 2> {log.err}
		bedtools intersect -a <(zcat {input.BED}) -b <(zcat {input.HetFiles}) -wao | cut -f1,2,3,8- | rev | cut -f2- | rev | awk '{{if(NR!=1 && reg!=$1"\\t"$2"\\t"$3){{str=reg; for(i=4;i<=NF;i++){{str=str"\\t"a[i]/len}} print str; delete a}} reg=$1"\\t"$2"\\t"$3; len=$3-$2+1; for(i=4;i<=NF;i++){{a[i]+=$i}}}}END{{str=reg; for(i=4;i<=NF;i++){{str=str"\\t"a[i]/len}} print str;}}' | gzip > ${{output}} 2>> {log.err}
		"""

# Compute heterozygosity for a set of windows with a fixed amount of callable sequence on them, for a group of samples
# to check
rule Heterozygosity_InEachRegion_FixedLengthRegions_PerSample:
	input:
		HetFiles = get_HetFiles,
		BED = rules.Windows_AccumulatedCallableSize_BEDs.output.WindowsBED,
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		out = "results/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_FixedLengthRegions_PerSample/Heterozygosity_InEachRegion_FixedLengthRegions_{windowsize}.txt.gz"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_FixedLengthRegions_PerSample/Heterozygosity_InEachRegion_FixedLengthRegions_{windowsize}.err",
		out = "logs/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_FixedLengthRegions_PerSample/Heterozygosity_InEachRegion_FixedLengthRegions_{windowsize}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Heterozygosity_InEachRegion_FixedLengthRegions_PerSample/Heterozygosity_InEachRegion_FixedLengthRegions_{windowsize}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '20:00:00',
		name = "Het_FixLen",
		threads = 1,
		mem = 100000,
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output%.gz}}) 2> {log.err}
		bedtools intersect -a <(zcat {input.BED}) -b <(zcat {input.HetFiles}) -wao | \
			cut -f1,2,3,8- | rev | cut -f2- | rev | \
			awk  -v len={wildcards.windowsize} '{{if(NR!=1 && reg!=$1"\\t"$2"\\t"$3){{str=reg; for(i=4;i<=NF;i++){{str=str"\\t"a[i]/len}} print str; delete a}} reg=$1"\\t"$2"\\t"$3; for(i=4;i<=NF;i++){{a[i]+=$i}}}}END{{str=reg; for(i=4;i<=NF;i++){{str=str"\\t"a[i]/len}} print str;}}' | \
			gzip > ${{output}} 2>> {log.err}
		"""
# Compute total average pairwise differences in a set of BED regions and a group of samples
rule Pi_Total_inBEDregions_PerSubsetOfSamples:
	'''
	General average pairwise differences in given regions and group of samples
	'''
	input:
		GenoFiles = get_GenoFiles,
		BED = get_bed_file_name,
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		out = "results/Plotting_DNA/{ObsExp}/Pi_Total_inBEDregions_PerSubsetOfSamples/Pi_Total_in{BED}Regions_Per{GroupSamples}.txt"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Pi_Total_inBEDregions_PerSubsetOfSamples/Pi_Total_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/Pi_Total_inBEDregions_PerSubsetOfSamples/Pi_Total_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Pi_Total_inBEDregions_PerSubsetOfSamples/Pi_Total_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '20:00:00',
		name = "PIt_{BED}_{GroupSamples}",
		threads = 1,
		mem = 10000,
		samples = get_group_of_samples
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output}}) 2> {log.err}
		samplelines=$(awk '{{if(NR==FNR){{a[$1]=1;next}}if(a[$1]){{str=str","FNR+4}}}}END{{print str}}' <(echo {params.samples} | sed 's/\\s/\\n/g') <(cat metadata/SamplesOrderInVCF.chr19.txt | sed 's/\\s/\\n/g') | sed 's/^,//g')
		totalsites=$(zcat {input.BED} | awk '{{len+=$3-$2+1}}END{{print len}}')
		bedtools intersect -a <(zcat {input.GenoFiles}) -b <(zcat {input.BED}) | cut -f${{samplelines}} | sed 's/:/\\t/g' | awk -v totalsites=${{totalsites}} '{{if(NR==1){{NS=NF;c=0;for(i=1;i<NS;i++){{for(j=i+1;j<=NS;j++){{c++;a[c]=0}}}}}}comp=0;for(i=1;i<NS;i++){{for(j=i+1;j<=NS;j++){{comp++;if($i!=$j){{a[comp]++}}}}}}}}END{{c=0;for(i=1;i<NS;i++){{for(j=i+1;j<=NS;j++){{c++;sum+=a[c]/totalsites}}}}print sum/c}}' 2> {log.err} > ${{output}}
		"""

# Compute average pairwise differences for each one of a set of BED regions, for a group of samples
# Fix: take regions into account as a whole, not divided in its callable regions
rule Pi_InEachRegion_inBEDregions_PerSubsetOfSamples:
	'''
	Average pairwise differences for each region, for a group of samples
	'''
	input:
		GenoFiles = get_GenoFiles,
		BED = get_bed_file_name,
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		out = "results/Plotting_DNA/{ObsExp}/Pi_InEachRegion_inBEDregions_PerSubsetOfSamples/Pi_InEachRegion_in{BED}Regions_Per{GroupSamples}.txt.gz"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Pi_InEachRegion_inBEDregions_PerSubsetOfSamples/Pi_InEachRegion_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/Pi_InEachRegion_inBEDregions_PerSubsetOfSamples/Pi_InEachRegion_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Pi_InEachRegion_inBEDregions_PerSubsetOfSamples/Pi_InEachRegion_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '20:00:00',
		name = "PIr_{BED}_{GroupSamples}",
		threads = 1,
		mem = 10000,
		samples = get_group_of_samples
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output%.gz}}) 2> {log.err}
		samplelines=$(awk '{{if(NR==FNR){{a[$1]=1;next}}if(a[$1]){{str=str","FNR+4}}}}END{{print str}}' <(echo {params.samples} | sed 's/\\s/\\n/g') <(cat metadata/SamplesOrderInVCF.chr19.txt | sed 's/\\s/\\n/g') | sed 's/^,//g')
		awk '{{if(FNR==NR){{a[$1"_"$2"_"$3]=$4;next}}if(a[$1"_"$2"_"$3]){{print a[$1"_"$2"_"$3]}}else{{print 0}}}}'  <(bedtools intersect -a <(zcat {input.GenoFiles} | cut -f1,2,3,${{samplelines}}) -b <(zcat {input.BED} | cut -f1,2,3) -wo | awk '{{chr=NF-3; st=NF-2; end=NF-1; str=$chr"\\t"$st"\\t"$end; for(i=4; i<chr; i++){{str=str"\\t"$i}} print str}}' | sed 's/:/\\t/g' | awk '{{if(NR==1){{ns=NF-3;}}if(chr"_"st"_"end != $1"_"$2"_"$3){{sum=0;for(i=4;i<ns+3;i++){{for(j=i+1;j<=ns+3;j++){{sum+=a[i"_"j]; a[i"_"j]=0}}}} if(NR>1){{if(sum>0){{sum=sum/length(a)/(end-st+1)}} print chr"\\t"st"\\t"end"\\t"sum}} chr=$1; st=$2; end=$3}}else{{for(i=4;i<ns+3;i++){{for(j=i+1;j<=ns+3;j++){{if($i!=$j){{a[i"_"j]++}}}}}}}}}}') <(zcat {input.BED}) | gzip 2> {log.err} > ${{output}} 
		"""
# Compute average pairwise differences for a set of windows with a fixed amount of callable sequence on them, for a group of samples
# to rerun
rule Pi_InEachRegion_FixedLengthRegions_PerSubsetOfSamples:
	input:
		GenoFiles = get_GenoFiles,
		BED = expand(rules.Windows_AccumulatedCallableSize_BEDs.output.WindowsBED, windowsize=config["windowsize"]),
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		out = "results/Plotting_DNA/{ObsExp}/Pi_InEachRegion_FixedLengthRegions_PerSubsetOfSamples/Pi_InEachRegion_FixedLengthRegions_Per{GroupSamples}.txt.gz"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Pi_InEachRegion_FixedLengthRegions_PerSubsetOfSamples/Pi_InEachRegion_FixedLengthRegions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/Pi_InEachRegion_FixedLengthRegions_PerSubsetOfSamples/Pi_InEachRegion_FixedLengthRegions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Pi_InEachRegion_FixedLengthRegions_PerSubsetOfSamples/Pi_InEachRegion_FixedLengthRegions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '20:00:00',
		name = "PI_FixLen_{GroupSamples}",
		threads = 1,
		mem = 10000,
		samples = get_group_of_samples,
		windowsize=config["windowsize"]
	shell:
		"""
		output={output.out}
		mkdir -p $(dirname ${{output%.gz}}) 2> {log.err}
		samplelines=$(awk '{{if(NR==FNR){{a[$1]=1;next}}if(a[$1]){{str=str","FNR+4}}}}END{{print str}}' <(echo {params.samples} | sed 's/\\s/\\n/g') <(cat metadata/SamplesOrderInVCF.chr19.txt | sed 's/\\s/\\n/g') | sed 's/^,//g')
		awk '{{if(FNR==NR){{a[$1"_"$2"_"$3]=$4;next}}if(a[$1"_"$2"_"$3]){{print a[$1"_"$2"_"$3]}}else{{print "0"}}}}' <(bedtools intersect -a <(zcat {input.GenoFiles} | cut -f1,2,3,${{samplelines}}) -b <(zcat {input.BED} | cut -f1,2,3) -wo | awk '{{chr=NF-3; st=NF-2; end=NF-1; str=$chr"\\t"$st"\\t"$end; for(i=4; i<chr; i++){{str=str"\\t"$i}} print str}}' | sed 's/:/\\t/g' | awk -v reglen={params.windowsize} '{{if(NR==1){{ns=NF-3;}}if(chr"_"st"_"end != $1"_"$2"_"$3){{sum=0;for(i=4;i<ns+3;i++){{for(j=i+1;j<=ns+3;j++){{sum+=a[i"_"j]; a[i"_"j]=0}}}} if(NR>1){{if(sum>0){{sum=sum/length(a)/reglen}} print chr"\\t"st"\\t"end"\\t"sum}} chr=$1; st=$2; end=$3}}else{{for(i=4;i<ns+3;i++){{for(j=i+1;j<=ns+3;j++){{if($i!=$j){{a[i"_"j]++}}}}}}}}}}') <(zcat {input.BED}) | gzip 2> {log.err} > ${{output}} 
		"""
# Build files with a list of samples corresponding to the wildcard RandObsSim (including simulated samples and randomized samples)
rule Prepare_PopulationSamples_Files:
	'''
	Calculate the frequencey of each allele in each posistion for a subset of regions (bed file) and a subset of samples (AllSamples, AtlSamples or MedSamples)
	'''
	input:
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		outAll = "results/Plotting_DNA/{ObsExp}/Prepare_PopulationSamples_Files/{RandObsSim}_AllSamples_{boots}.txt",
		outAtl = "results/Plotting_DNA/{ObsExp}/Prepare_PopulationSamples_Files/{RandObsSim}_AtlSamples_{boots}.txt",
		outMed = "results/Plotting_DNA/{ObsExp}/Prepare_PopulationSamples_Files/{RandObsSim}_MedSamples_{boots}.txt"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Prepare_PopulationSamples_Files/{RandObsSim}_PopulationSamples_{boots}.err",
		out = "logs/Plotting_DNA/{ObsExp}/Prepare_PopulationSamples_Files/{RandObsSim}_PopulationSamples_{boots}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Prepare_PopulationSamples_Files/{RandObsSim}_PopulationSamples_{boots}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "Prepare_PopFile_{RandObsSim}_{boots}",
		threads = 1,
		mem = 10000,
		samplesAtl = config["AtlSamples"],
		samplesMed = config["MedSamples"],
		SimSSize = config["SimulationsSampleSize"]
	shell:
		"""
		mkdir -p $(dirname {output.outAll}) 2> {log.err}
		RandObsSim={wildcards.RandObsSim}
		if [[ ${{RandObsSim}} == "Rand" ]]
		then
			numA=$(echo {params.samplesAtl} | sed 's/ /\\n/g' | awk '{{}}END{{print NR}}' )
			numM=$(echo {params.samplesMed} | sed 's/ /\\n/g' | awk '{{}}END{{print NR}}' )
			shufled=$(cat {input.SamplesOrderInVCF} | sed 's/\\t/\\n/g' | shuf)
			echo ${{shufled}} | sed 's/ /\\n/g'  > {output.outAll}
			echo ${{shufled}} | sed 's/ /\\n/g' | tac | head -${{numA}} > {output.outAtl}
			echo ${{shufled}} | sed 's/ /\\n/g' | head -${{numM}} > {output.outMed}
		elif [[ ${{RandObsSim}} == "Obs" ]]
		then
			echo {params.samplesMed} {params.samplesAtl} | sed 's/ /\\n/g' > {output.outAll}
			echo {params.samplesAtl} | sed 's/ /\\n/g' > {output.outAtl}
			echo {params.samplesMed} | sed 's/ /\\n/g' > {output.outMed}
		else
			if [[ -s {output.outAll} ]];then rm {output.outAll}; fi
			for i in `seq 0 1 {params.SimSSize}`; do echo $i >> {output.outAll}; done
			touch {output.outAtl} {output.outMed}
		fi
		"""

# Compute total callable spand of exons, introns, promoters and intergenic regions
rule get_FunctionalRegionsCallableSpan:
	'''
	'''
	input:
		ExonsBED=expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Exons"),
		IntronsBED=expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Introns"),
		PromotersBED=expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Promoters"),
		IntergenicBED=expand(rules.Intersection2BEDs.output.Intersection, Name1="Callable", Name2="Intergenic"),
		CallableBED=rules.r7_join_extra_callable_regions.output.ExtraCallableRegions
	output:
		out = "results/Plotting_DNA/get_FunctionalRegionsCallableSpan/FunctionalRegionsCallableSpan.txt"
	log:
		err = "logs/Plotting_DNA/get_FunctionalRegionsCallableSpan/get_FunctionalRegionsCallableSpan.err",
		out = "logs/Plotting_DNA/get_FunctionalRegionsCallableSpan/get_FunctionalRegionsCallableSpan.out"
	benchmark:
		"benchmarks/Plotting_DNA/get_FunctionalRegionsCallableSpan/get_FunctionalRegionsCallableSpan.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "FuncSpan",
		threads = 1,
		mem = 50000
	shell:
		"""
		if [[ -s {output.out} ]]; then rm {output.out}; fi
		mkdir -p $(dirname {output.out}) 2>> {log.err}
		labels=("Exons" "Introns" "Promoters" "Intergenic" "Callable")
		inputs=({input.ExonsBED} {input.IntronsBED} {input.PromotersBED} {input.IntergenicBED} {input.CallableBED})
		for l in ${{!labels[@]}}
		do
			echo ${{labels[$l]}} ${{inputs[$l]}}
			length=$( zcat ${{inputs[$l]}} | awk '{{len=len+$3-$2+1}}END{{print len}}' )
			echo ${{labels[$l]}}"\t"${{length}} >> {output.out}
		done
		"""

# Compute PCA for short variants genotypes for a specific set of variants 
rule PCA_short_variants_inBEDregions_PerSubsetOfSamples:
	'''
	'''
	input:
		GenoFiles = expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPE, chr=config["bralan3chrs"]),
		chrsGENOTYPEMatrices = expand(rules.Prepare_basic_analysis_files_PerChr.output.chrGENOTYPEMatrix, chr=config["bralan3chrs"]),
		BED = get_bed_file_name,
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt" 
	output:
		PCA_PerSample = "results/Plotting_DNA/{ObsExp}/PCA_short_variants_inBEDregions_PerSubsetOfSamples/PCA_PerSample_in{BED}Regions_Per{GroupSamples}.txt.gz",
		PCA_PerVariant = "results/Plotting_DNA/{ObsExp}/PCA_short_variants_inBEDregions_PerSubsetOfSamples/PCA_PerVariant_in{BED}Regions_Per{GroupSamples}.txt.gz",
		PCA_PropVariance = "results/Plotting_DNA/{ObsExp}/PCA_short_variants_inBEDregions_PerSubsetOfSamples/PCA_PropVariance_in{BED}Regions_Per{GroupSamples}.txt"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/PCA_short_variants_inBEDregions_PerSubsetOfSamples/PCA_short_variants_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/PCA_short_variants_inBEDregions_PerSubsetOfSamples/PCA_short_variants_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/PCA_short_variants_inBEDregions_PerSubsetOfSamples/PCA_short_variants_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '10:00:00',
		name = "cPCA_{BED}_{GroupSamples}",
		threads = 1,
		mem = 500000,
		samples = get_group_of_samples,
		chrs = config["bralan3chrs"]
	shell:
		"""
		outputS={output.PCA_PerSample}
		outputV={output.PCA_PerVariant}
		mkdir -p $(dirname ${{outputV%.gz}}) 2> {log.err}
		nummat1=$(echo {input.chrsGENOTYPEMatrices} | awk '{{split($1,a,"chr1"); print a[1]}}')
		nummat2=$(echo {input.chrsGENOTYPEMatrices} | awk '{{split($1,a,"chr1"); print a[2]}}')
		genof1=$(echo {input.GenoFiles} | awk '{{split($1,a,"chr1"); print a[1]}}')
		genof2=$(echo {input.GenoFiles} | awk '{{split($1,a,"chr1"); print a[2]}}')
		for c in {params.chrs}
		do
			echo $c > {log.out}
			awk '{{if(NR==FNR){{a[$1]=1; next}}if(a[$1]){{print $0}}}}' <(bedtools intersect -a <(zcat ${{genof1}}${{c}}${{genof2}} | cut -f1,2,3,4) -b <(zcat {input.BED} | cut -f1,2,3 | grep -w ${{c}}) -wa | cut -f4) <(zcat ${{nummat1}}${{c}}${{nummat2}}) | gzip > ${{outputV%.txt.gz}}_IncludedVariants_${{c}}.tmp.gz 2> {log.err}
		done
		zcat ${{outputV%.txt.gz}}_IncludedVariants* | wc -l
		./scripts/Plotting_DNA/PCA_short_variants_inBEDregions_PerSubsetOfSamples.R \
		${{outputV%.txt.gz}}_IncludedVariants \
		{input.SamplesOrderInVCF} \
		${{outputS%.gz}} \
		${{outputV%.gz}} \
		{output.PCA_PropVariance} \
		\"{params.chrs}\" \
		\"{params.samples}\" > {log.out} 2> {log.err}
		gzip ${{outputS%.gz}} ${{outputV%.gz}}
		rm ${{outputV%.txt.gz}}_IncludedVariants*
		"""



# 
rule Join_PSMC_results:
	'''
	PSMC Pairwise Sequentially Markovian Coalescent
	Li H, Durbin R. Inference of human population history from individual whole-genome sequences. Nature. 2011 Jul 13;475(7357):493-6. doi: 10.1038/nature10231. PMID: 21753753; PMCID: PMC3154645.
	'''
	input:
		psmc=expand(rules.PSMC_PerSample.output.psmc, sample=config["samples"]),
	output:
		jpsmc="results/Plotting_DNA/Join_PSMC_results/joined_PSMC.txt",
		jpsmc_thetas="results/Plotting_DNA/Join_PSMC_results/joined_PSMC_thetas.txt",
	log:
		err = "logs/Plotting_DNA/Join_PSMC_results/Join_PSMC_results.err",
		out = "logs/Plotting_DNA/Join_PSMC_results/Join_PSMC_results.out"
	benchmark:
		"benchmarks/Plotting_DNA/Join_PSMC_results/Join_PSMC_results.txt"
	conda:
		'../envs/hierfstat.yaml'
	params:
		time = '3:00:00',
		name = "jPSMCWG",
		threads = 1,
		mem = 10000,
		samples = expand("{sample}", sample=config["samples"]),
		Blangenerationtime = config["Blangenerationtime"],
		rounditerationsPSMC = config["rounditerationsPSMC"],
		highmu = config["highmuPSMC"],
		lowmu = config["lowmuPSMC"],
		step = config["stepPSMC"],
	shell:
		"./scripts/Plotting_DNA/Join_PSMC_results.sh \"{input.psmc}\" {output.jpsmc}  {output.jpsmc_thetas} {params.Blangenerationtime} \"{params.samples}\" {params.rounditerationsPSMC} {params.highmu} {params.lowmu} {params.step} > {log.out} 2> {log.err}"

#
rule Coverage_inAllRegions_PerSample:
	'''
	Calculate average and standard deviation of coverage for each sample in Callable regions
	'''
	input:
		CoverageFiles = rules.r7_join_coverage_per_site.output.covgenome,
	output:
		out = "results/Plotting_DNA/{ObsExp}/Coverage_inAllRegions_PerSample.txt"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Coverage_inAllRegions_PerSample.err",
		out = "logs/Plotting_DNA/{ObsExp}/Coverage_inAllRegions_PerSample.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Coverage_inAllRegions_PerSample.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '03:00:00',
		name = "cov",
		threads = 1,
		mem = 500000
	shell:
		"""
		zcat < {input.CoverageFiles} | cut -f3- | awk '
			{{
			    for (i = 1; i <= NF; i++) {{
			        sum[i] += $i
			        sumsq[i] += $i * $i
			        n[i]++
			    }}
			}}
			END {{
			    for (i = 1; i <= length(sum); i++) {{
			        mean = sum[i] / n[i]
			        sd = sqrt((sumsq[i] / n[i]) - mean^2)
			        print mean\"\\t\"sd
			    }}
			}}
			' > {output}
		"""


#
# Compute total population Fst
rule Population_Fst_Total:
	'''
	Population Fst
	'''
	input:
		VCF_File = rules.r9_run_all_chr.output.passVCF,
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
	output:
		outPerSite = "results/Plotting_DNA/{ObsExp}/Population_Fst_Total/Population_Fst_PerSite.weir.fst",
		outTotal = "results/Plotting_DNA/{ObsExp}/Population_Fst_Total/Population_Fst_Total.txt"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/Population_Fst_Total/Population_Fst_Total.err",
		out = "logs/Plotting_DNA/{ObsExp}/Population_Fst_Total/Population_Fst_Total.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/Population_Fst_Total/Population_Fst_Total.txt"
	conda:
		'../envs/VariantAnalysis_DNA.yaml'
	params:
		time = '20:00:00',
		name = "Fst_perVar",
		threads = 1,
		mem = 100000
	shell:
		"""
		outputT={output.outTotal}
		outputPS={output.outPerSite}
		mkdir -p $(dirname ${{outputT}}) 2> {log.err}
		tail -n +2 {input.Metadata} | grep -w "Banyuls" | cut -f3 > $(dirname ${{outputT}})/MediterraneanSamples.txt
		tail -n +2 {input.Metadata} | grep -w "Roscoff" | cut -f3 > $(dirname ${{outputT}})/AtlanticSamples.txt
		echo "----> Fst by population per site" # only considers PASS variants
		vcftools --gzvcf {input.VCF_File} --remove-filtered-all --weir-fst-pop $(dirname ${{outputT}})/AtlanticSamples.txt --weir-fst-pop $(dirname ${{outputT}})/MediterraneanSamples.txt --out ${{outputPS%.weir.fst}} 2>> {log.err}
		echo "----> Averaging Fst by population to get total Fst"
		awk 'NR>1 && $3 != "nan" {{sum+=$3; n++}} END {{print sum/n}}' ${{outputPS}} > ${{outputT}} 2>> {log.err}
		"""

################################################
## 1. Estimate and describe genomic diversity
################################################
# add synonymous/non-synonymous
rule plot_Figure1_GenomicDiversity:
	'''

	'''
	input:
		Rconfig = config["Rconfig"],
		PiPerRegFiles = expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, ObsExp=["Observed_Data", "Observed_SNPs", "Observed_INDELs"], BED=["Exons", "Introns", "Promoters", "Intergenic"], GroupSamples=["AtlSamples", "MedSamples", "AllSamples"]),
		PiTotalFiles = expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, ObsExp=["Observed_Data", "Observed_SNPs", "Observed_INDELs"], BED=["Exons", "Introns", "Promoters", "Intergenic", "Callable"], GroupSamples=["AtlSamples", "MedSamples", "AllSamples"]),
		HetPerRegFiles = expand(rules.Heterozygosity_InEachRegion_inBEDregions_PerSample.output.out, ObsExp=["Observed_Data", "Observed_SNPs", "Observed_INDELs"], BED=["Exons", "Introns", "Promoters", "Intergenic"]),
		HetTotalFiles = expand(rules.Heterozygosity_Total_inBEDregions_PerSample.output.out, ObsExp=["Observed_Data", "Observed_SNPs", "Observed_INDELs"], BED=["Exons", "Introns", "Promoters", "Intergenic", "Callable"]),
		HetFixedWindowsFiles = expand(rules.Heterozygosity_InEachRegion_FixedLengthRegions_PerSample.output.out, ObsExp=["Observed_Data", "Observed_SNPs", "Observed_INDELs"], windowsize=[50,100,200]),
		FreqPerSiteFiles = expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, ObsExp=["Observed_Data", "Observed_SNPs", "Observed_INDELs"], BED=["Callable", "SNPs", "INDELs", "Exons", "Introns", "Promoters", "Intergenic"], GroupSamples=["AtlSamples", "MedSamples", "AllSamples"]),
		CoverageFile = expand(rules.Coverage_inAllRegions_PerSample.output.out, ObsExp=["Observed_Data"]),
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19"),
		FunctionalRegionsCallableSpan = rules.get_FunctionalRegionsCallableSpan.output.out,
		PerSiteFeatureType = rules.FunctionalFeatures_BEDs.output.PerSite_FeatureType,
		SpeciesTree = rules.Build_Tree_Species.output.tree,
		MapImage = "metadata/MapSampling.png",
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		Lynch2023 = "data/OtherDiversityEstimates/SuppMat_Lynch2023.txt",
		Leffler2012 = "data/OtherDiversityEstimates/SuppMat_Leffler2012.txt",
		Romiguier2014 = "data/OtherDiversityEstimates/SuppMat_Romiguier2014.txt",
		CorbettDetig2015 = "data/OtherDiversityEstimates/SuppMat_Corbett-Detig2015.txt",
		IndependentEstimates = "data/OtherDiversityEstimates/Independent_GenomicDiversityEstimates.tab"
	output:
		PDF = "results/Plotting_DNA/plot_Figure1_GenomicDiversity.pdf",
		REPORT = "results/Plotting_DNA/plot_Figure1_GenomicDiversity_report.txt"
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
		mem = 10000,
		AtlSamples=config["AtlSamples"],
		MedSamples=config["MedSamples"]
	shell:
		"""
		./scripts/Plotting_DNA/plot_Figure1_GenomicDiversity.R \
		\"{input.PiPerRegFiles}\" \
		\"{input.PiTotalFiles}\" \
		\"{input.HetPerRegFiles}\"  \
		\"{input.HetTotalFiles}\"  \
		\"{input.HetFixedWindowsFiles}\"  \
		\"{input.FreqPerSiteFiles}\"  \
		\"{input.CoverageFile}\"  \
		{input.SamplesOrderInVCF}  \
		{input.FunctionalRegionsCallableSpan}  \
		{input.PerSiteFeatureType}  \
		{input.SpeciesTree}  \
		{input.MapImage}  \
		{input.Metadata}  \
		{input.Lynch2023} \
		{input.Leffler2012} \
		{input.Romiguier2014} \
		{input.CorbettDetig2015} \
		{input.IndependentEstimates} \
		{output.PDF}  \
		{output.REPORT}  \
		\"{params.AtlSamples}\"  \
		\"{params.MedSamples}\"  \
		{input.Rconfig} > {log.out} 2> {log.err}
		"""
#
rule plot_Figure2_PopulationStructure:
	'''
	
	Plotting PSMC Pairwise Sequentially Markovian Coalescent
	Li H, Durbin R. Inference of human population history from individual whole-genome sequences. Nature. 2011 Jul 13;475(7357):493-6. doi: 10.1038/nature10231. PMID: 21753753; PMCID: PMC3154645.
	'''
	input:
		Rconfig = config["Rconfig"],
		PCAValues_files = expand(rules.PCA_short_variants_inBEDregions_PerSubsetOfSamples.output.PCA_PerSample, ObsExp="Observed_Data", BED=("Callable", "SNPs", "INDELs", "Exons", "Introns", "Promoters", "Intergenic", "BiAllelic", "AbsFreq_36-40_AllVariants", "AbsFreq_40-45_AllVariants", "AbsFreq_45-50_AllVariants", "AbsFreq_50-55_AllVariants", "AbsFreq_55-60_AllVariants", "AbsFreq_60-65_AllVariants", "AbsFreq_65-70_AllVariants", "AbsFreq_70-73_AllVariants"), GroupSamples="AllSamples"),
		PCAPropVar_files = expand(rules.PCA_short_variants_inBEDregions_PerSubsetOfSamples.output.PCA_PropVariance, ObsExp="Observed_Data", BED=("Callable", "SNPs", "INDELs", "Exons", "Introns", "Promoters", "Intergenic", "BiAllelic", "AbsFreq_36-40_AllVariants", "AbsFreq_40-45_AllVariants", "AbsFreq_45-50_AllVariants", "AbsFreq_50-55_AllVariants", "AbsFreq_55-60_AllVariants", "AbsFreq_60-65_AllVariants", "AbsFreq_65-70_AllVariants", "AbsFreq_70-73_AllVariants"), GroupSamples="AllSamples"),
		ObsShPriv = expand(rules.SharedPrivatePopulation_BEDs.output.Numbers, ObsOrBoots="Observed"),
		RandShPriv = expand(rules.SharedPrivatePopulation_BEDs.output.Numbers, ObsOrBoots=["Rand_%d" % i for i in range(config["BootsRandPop"])]),
		jpsmc= rules.Join_PSMC_results.output.jpsmc,
		jpsmc_thetas= rules.Join_PSMC_results.output.jpsmc_thetas,
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt",
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		DistTree = rules.Build_DistanceTree_Samples.output.tree,
	output:
		PDF = "results/Plotting_DNA/plot_Figure2_PopulationStructure.pdf",
		REPORT = "results/Plotting_DNA/plot_Figure2_PopulationStructure_report.txt"
	log:
		err = "logs/Plotting_DNA/plot_Figure2_PopulationStructure.err",
		out = "logs/Plotting_DNA/plot_Figure2_PopulationStructure.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_Figure2_PopulationStructure.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pFig2",
		threads = 1,
		mem = 50000,
		chrs = config["bralan3chrs"],
		PSMChighmu = config["highmuPSMC"],
		PSMClowmu = config["lowmuPSMC"],
		Blangenerationtime = config["Blangenerationtime"],
	shell:
		"""
		./scripts/Plotting_DNA/plot_Figure2_PopulationStructure.R \
		\"{input.PCAValues_files}\" \
		\"{input.PCAPropVar_files}\" \
		\"{input.ObsShPriv}\" \
		\"{input.RandShPriv}\" \
		{input.jpsmc} \
		{input.jpsmc_thetas} \
		{input.SamplesOrderInVCF} \
		{input.Metadata} \
		{input.DistTree} \
		{output.PDF} \
		{output.REPORT} \
		\"{params.chrs}\" \
		{params.PSMChighmu} \
		{params.PSMClowmu} \
		{params.Blangenerationtime} \
		{input.Rconfig} > {log.out} 2> {log.err}
		"""
# 
rule plot_Figure3_vsSimulations:
	'''

	'''
	input:
		Rconfig = config["Rconfig"],
		ObsHet = expand(rules.Heterozygosity_InEachRegion_inBEDregions_PerSample.output.out, BED="SelectedWindowsSimChrSize", GroupSamples="AllSamples", ObsExp="Observed_Data"),
		SimHet = expand(rules.Heterozygosity_Total_inBEDregions_PerSample.output.out, BED=["Simulated_%d_FiniteAlleles" % i for i in config["SimulationsReplicates"]], GroupSamples="SimSamples", ObsExp=["Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_%s" % i for i in config["SimulationsRangeMuNe"]]),
		ObsSFSAtl = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutMajorFreq, BED="SelectedWindowsSimChrSize", GroupSamples="AtlSamples", ObsExp="Observed_SNPs"),
		ObsSFSMed = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutMajorFreq, BED="SelectedWindowsSimChrSize", GroupSamples="MedSamples", ObsExp="Observed_SNPs"),
		SimSFS = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutMajorFreq, BED=["Simulated_%d_FiniteAlleles" % i for i in config["SimulationsReplicates"]], GroupSamples="SimSamples", ObsExp=["Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_%s" % i for i in config["SimulationsRangeMuNe"]]),
		ObsNAllelesAtl = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutNumAlleles, BED="SelectedWindowsSimChrSize", GroupSamples="AtlSamples", ObsExp="Observed_SNPs"),
		ObsNAllelesMed = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutNumAlleles, BED="SelectedWindowsSimChrSize", GroupSamples="MedSamples", ObsExp="Observed_SNPs"),
		SimNAlleles = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutNumAlleles, BED=["Simulated_%d_FiniteAlleles" % i for i in config["SimulationsReplicates"]], GroupSamples="SimSamples", ObsExp=["Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_%s" % i for i in config["SimulationsRangeMuNe"]]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt",
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
	output:
		PDF = "results/Plotting_DNA/plot_Figure3_vsSimulations.pdf",
		REPORT = "results/Plotting_DNA/plot_Figure3_vsSimulations_report.txt"
	log:
		err = "logs/Plotting_DNA/plot_Figure3_vsSimulations.err",
		out = "logs/Plotting_DNA/plot_Figure3_vsSimulations.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_Figure3_vsSimulations.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pFig3",
		threads = 1,
		mem = 50000,
		SimulationsSampleSize = config["SimulationsSampleSize"]
	shell:
		"""
		./scripts/Plotting_DNA/plot_Figure3_vsSimulations.R \
		{input.ObsHet} \
		\"{input.SimHet}\" \
		{input.ObsSFSAtl} \
		{input.ObsSFSMed} \
		\"{input.SimSFS}\" \
		{input.ObsNAllelesAtl} \
		{input.ObsNAllelesMed} \
		\"{input.SimNAlleles}\" \
		{input.SamplesOrderInVCF} \
		{input.Metadata} \
		{output.PDF} \
		{output.REPORT} \
		{params.SimulationsSampleSize} \
		{input.Rconfig} > {log.out} 2> {log.err}
		"""
#
rule ENA_Submision:
	'''
	Following instructions from:
		- https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html
		- https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html
	Fields in ManifestFile:
		STUDY
		SAMPLE
		NAME
		PLATFORM
		INSTRUMENT
		INSERT_SIZE
		LIBRARY_NAME (optional)
		LIBRARY_SOURCE
		LIBRARY_SELECTION
		LIBRARY_STRATEGY
		DESCRIPTION (optional)
		FASTQ
	java -jar webin-cli-<version>.jar
		-context reads 
		-userName: the Webin submission account name.
		-password: the Webin submission account password.
		-centerName: the center name of the submitter (mandatory for broker accounts).
		-manifest: the manifest file name.
		-outputDir: directory for output files.
		-inputDir: input directory for files declared in manifest file.
		-validate: validates the files defined in the manifest file.
		-submit: validates and submits the files defined in the manifest file.
		-test: use Webin test service
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		SamplesENA = "results/Plotting_DNA/ENA_Submision/samples-ENA.csv",
		sampleFASTQ1 = "data/DNAseqFASTQ/{sample}_{lane}_{platform}_R1.fastq.gz",
		sampleFASTQ2 = "data/DNAseqFASTQ/{sample}_{lane}_{platform}_R2.fastq.gz",
	output:
		ManifestFile = "results/Plotting_DNA/ENA_Submision/ManifestFile_{sample}_{lane}_{platform}.txt",
		SubmissionReport = "results/Plotting_DNA/ENA_Submision/SubmissionReport_{sample}_{lane}_{platform}.txt",
	log:
		err = "logs/Plotting_DNA/ENA_Submision_{sample}_{lane}_{platform}.err",
		out = "logs/Plotting_DNA/ENA_Submision_{sample}_{lane}_{platform}.out"
	benchmark:
		"benchmarks/Plotting_DNA/ENA_Submision_{sample}_{lane}_{platform}.txt"
	conda:
		'../envs/ENAsubmission.yaml'
	params:
		time = '1:00:00',
		name = "ENA{sample}_{lane}_{platform}",
		threads = 1,
		mem = 5000,
		ENAjavafile = config["ENAjavafile"],
		ENAuserName = config["ENAuserName"],
		ENApassword = config["ENApassword"],
	shell:
		"""
		pwd=$(pwd)
		# Build Manifest file per sample, lane & platform
		# STUDY
		echo "STUDY\tPRJEB106885" > {output.ManifestFile} 2> {log.err}
		# SAMPLE
		ENAsample=$(cat {input.SamplesENA} | grep ,{wildcards.sample}, | cut -f1 -d',') 2> {log.err}
		echo "SAMPLE\t$ENAsample" >> {output.ManifestFile} 2> {log.err}
		# NAME
		echo "NAME\tEXP_{wildcards.sample}_{wildcards.lane}_{wildcards.platform}" >> {output.ManifestFile} 2> {log.err}
		# INSTRUMENT
		awk -v platform={wildcards.platform} 'BEGIN {{if (platform == \"HiSeq4000\") print \"INSTRUMENT\tIllumina HiSeq 4000\"; else if (platform == \"NovaSeq6000\") print \"INSTRUMENT\tIllumina NovaSeq 6000\"}}' >> {output.ManifestFile} 2> {log.err}
		# INSERT_SIZE
		echo "INSERT_SIZE\t300" >> {output.ManifestFile} 2> {log.err}
		# LIBRARY_SOURCE
		echo "LIBRARY_SOURCE\tGENOMIC" >> {output.ManifestFile} 2> {log.err}
		# LIBRARY_SELECTION
		echo "LIBRARY_SELECTION\tRANDOM" >> {output.ManifestFile} 2> {log.err}
		# LIBRARY_STRATEGY
		echo "LIBRARY_STRATEGY\tWGS" >> {output.ManifestFile} 2> {log.err}
		# FASTQ
		echo "FASTQ\t$(basename {input.sampleFASTQ1})" >> {output.ManifestFile} 2> {log.err}
		echo "FASTQ\t$(basename {input.sampleFASTQ2})" >> {output.ManifestFile} 2> {log.err}
		cp {output.ManifestFile} {output.SubmissionReport}


		# Uploading fastq files
		echo "Uploading fastq files..." >> {log.out}
		cd $(dirname {input.sampleFASTQ1})
		echo ascp -QT -l300M -L- $(basename {input.sampleFASTQ1}) $(basename {input.sampleFASTQ2}) {params.ENAuserName}@webin.ebi.ac.uk:.
		cd $pwd
		echo "Done uploading fastq files." >> {log.out}

		"""
#
rule ENA_Submision_AllSamplesLanesPlatforms:
	'''
	'''
	input:
		ManifestFiles_HiSeq4000 = expand("results/Plotting_DNA/ENA_Submision/ManifestFile_{samplelane}_HiSeq4000.txt", samplelane=config["samplelanesHiSeq4000"]),
		ManifestFiles_NovaSeq6000_1 = expand("results/Plotting_DNA/ENA_Submision/ManifestFile_{samplelane}_NovaSeq6000.txt", samplelane=config["samplelanesNovaSeq6000_1"]),
		ManifestFiles_NovaSeq6000_2 = expand("results/Plotting_DNA/ENA_Submision/ManifestFile_{samplelane}_NovaSeq6000.txt", samplelane=config["samplelanesNovaSeq6000_2"]),
	output:
		OutputFile = "results/Plotting_DNA/ENA_Submision/ManifestFile_all_samplelaneplatform.txt",
	log:
		err = "logs/Plotting_DNA/ENA_Submision_all.err",
		out = "logs/Plotting_DNA/ENA_Submision_all.out"
	benchmark:
		"benchmarks/Plotting_DNA/ENA_Submision_all.txt"
	conda:
		'../envs/ENAsubmission.yaml'
	params:
		time = '1:00:00',
		name = "ENA_all",
		threads = 1,
		mem = 5000,
	shell:
		"""
		echo "ManifestFile\tExperimentAccession\tRunAccession" > {output.OutputFile} 2> {log.err}
		for run in {input.ManifestFiles_HiSeq4000} {input.ManifestFiles_NovaSeq6000_1} {input.ManifestFiles_NovaSeq6000_2}; do
			samplelaneplatform=$(basename $run | cut -f3,4,5 -d'_' | sed 's/.txt//g')
			expid=$(cat $run | grep 'The following experiment accession was assigned to the submission:' | rev | cut -f1 -d' ' | rev) 2> {log.err}
			runid=$(cat $run | grep 'The following run accession was assigned to the submission:' | rev | cut -f1 -d' ' | rev) 2> {log.err}
			echo "$run\t$expid\t$runid" >> {output.OutputFile} 2> {log.err}
		done
		"""
#
rule scp_fastq_data_from_NAS:
	'''
	Copy fastq data from NAS (cluster) to local data folder (meant only to upload data to ENA)
	'''
	output:
		fastqfile = "data/DNAseqFASTQ/{idfile}.fastq.gz"
	log:
		err = "logs/Plotting_DNA/scp_fastq_data_from_NAS/{idfile}.err",
		out = "logs/Plotting_DNA/scp_fastq_data_from_NAS/{idfile}.out"
	benchmark:
		"benchmarks/Plotting_DNA/scp_fastq_data_from_NAS/{idfile}.txt"
	shell:
		"scp mbrasovi@curnagl.dcsr.unil.ch:/nas/FAC/FBM/DEE/mrobinso/default/D2c/mbrasovi/Banyuls_Roscoff/DNAseqFASTQ/{wildcards.idfile}.fastq.gz {output.fastqfile}"
#




############################
#### Deprecated plots
rule plot_PCA_short_variants_perGenomicFeature:
	'''
	Plot PCA dividing variants per genomic feature (exonic, intronic, promoter, intergenic) for all chrs
	'''
	input:
		Rconfig = config["Rconfig"],
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		chrsGENOTYPEMatrices = expand("results/VariantAnalysis_DNA/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.genotype.numericmatrix", chr=config["bralan3chrs"]),
		chrsDividedListExon = expand("results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.classified.{chr}_exon.multiallelic.txt", chr=config["bralan3chrs"]),
		chrsDividedListBiallelic = expand("results/VariantAnalysis_DNA/Divide_variants_perGenomicFeature_PerChr/ShortVariants.filtered.classified.{chr}.biallelic.txt", chr=config["bralan3chrs"]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt" 
	output:
		PCA = "results/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature/plot_PCA_short_variants_perGenomicFeature.pdf",
		PCnumericmatrix = "results/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature/Numericmatrix_PCvalues_PerVariant.txt"
	log:
		err = "logs/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature/plot_PCA_short_variants_perGenomicFeature.err",
		out = "logs/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature/plot_PCA_short_variants_perGenomicFeature.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature/plot_PCA_short_variants_perGenomicFeature.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '10:00:00',
		name = "PCA",
		threads = 1,
		mem = 500000,
		chrs = expand("{chr}", chr=config["bralan3chrs"])
	shell:
		"./scripts/Plotting_DNA/plot_PCA_short_variants_perGenomicFeature.sh {input.Metadata} \"{input.chrsGENOTYPEMatrices}\" \"{input.chrsDividedListExon}\" {input.SamplesOrderInVCF} {output.PCA} \"{params.chrs}\" {input.Rconfig} > {log.out} 2> {log.err}"

rule plot_join_PSMC:
	'''
	Plotting PSMC Pairwise Sequentially Markovian Coalescent
	Li H, Durbin R. Inference of human population history from individual whole-genome sequences. Nature. 2011 Jul 13;475(7357):493-6. doi: 10.1038/nature10231. PMID: 21753753; PMCID: PMC3154645.
	'''
	input:
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
		psmc=expand("results/VariantAnalysis_DNA/PSMC_PerSample/{sample}.psmc", sample=config["samples"]),
		SeaLevel = "data/Earth_GlacialCycles/SupplementaryTable1_Miller2020.tbl"
	output:
		psmc_joinTXT="results/Plotting_DNA/plot_join_PSMC/joined_PSMC.txt",
		psmc_joinPDF="results/Plotting_DNA/plot_join_PSMC/joined_PSMC.pdf"
	log:
		err = "logs/Plotting_DNA/plot_join_PSMC/plot_join_PSMC.err",
		out = "logs/Plotting_DNA/plot_join_PSMC/plot_join_PSMC.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_join_PSMC/plot_join_PSMC.txt"
	conda:
		'../envs/hierfstat.yaml'
	params:
		time = '3:00:00',
		name = "pPSMCWG",
		threads = 1,
		mem = 10000,
		samples = expand("{sample}", sample=config["samples"]),
		Blangenerationtime = config["Blangenerationtime"]
	shell:
		"./scripts/Plotting_DNA/plot_join_PSMC.sh {input.Metadata} \"{input.psmc}\" {input.SeaLevel} {output.psmc_joinTXT} {output.psmc_joinPDF} {params.Blangenerationtime} \"{params.samples}\" > {log.out} 2> {log.err}"


rule plot_PerSite_inBEDregions_PerSubsetOfSamples:
	'''
	Plot:
	- Site Frequency Spectrum
	- Number of alleles histogram
	'''
	input:
		Rconfig = config["Rconfig"],
		Freqfile = rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out
	output:
		PDF = "results/Plotting_DNA/{ObsExp}/plot_PerSite_inBEDregions_PerSubsetOfSamples/plot_PerSite_in{BED}Regions_Per{GroupSamples}.pdf"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/plot_PerSite_inBEDregions_PerSubsetOfSamples/plot_PerSite_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/plot_PerSite_inBEDregions_PerSubsetOfSamples/plot_PerSite_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/plot_PerSite_inBEDregions_PerSubsetOfSamples/plot_PerSite_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pSite{BED}_{GroupSamples}",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_PerSite_inBEDregions_PerSubsetOfSamples.R {input.Freqfile} {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"

# in dev
# plot het distribution across samples
rule plot_Heterozygosity_inBEDregions_PerSubsetOfSamples:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		Freqfile = rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out
	output:
		PDF = "results/Plotting_DNA/{ObsExp}/plot_Heterozygosity_inBEDregions_PerSubsetOfSamples/plot_Heterozygosity_in{BED}Regions_Per{GroupSamples}.pdf"
	log:
		err = "logs/Plotting_DNA/{ObsExp}/plot_Heterozygosity_inBEDregions_PerSubsetOfSamples/plot_Heterozygosity_in{BED}Regions_Per{GroupSamples}.err",
		out = "logs/Plotting_DNA/{ObsExp}/plot_Heterozygosity_inBEDregions_PerSubsetOfSamples/plot_Heterozygosity_in{BED}Regions_Per{GroupSamples}.out"
	benchmark:
		"benchmarks/Plotting_DNA/{ObsExp}/plot_Heterozygosity_inBEDregions_PerSubsetOfSamples/plot_Heterozygosity_in{BED}Regions_Per{GroupSamples}.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pSite{BED}_{GroupSamples}",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_Heterozygosity_in.R {input.Freqfile} {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"


############################
####
rule plot_HeterozygosityPi_InPop:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		HetAll=expand(rules.Heterozygosity_Total_inBEDregions_PerSample.output.out, BED="Callable", ObsExp="Observed_Data"),
		PiTObs=expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		SamplesOrderInVCF = expand(rules.r8_filter_VCF_for_callable_regions_PerChr.output.SamplesOrderInVCF, chr="chr19")
	output:
		PDF = "results/Plotting_DNA/plot_HeterozygosityPi_InPop/plot_HeterozygosityPi_InPop.pdf"
	log:
		err = "logs/Plotting_DNA/plot_HeterozygosityPi_InPop/plot_HeterozygosityPi_InPop.err",
		out = "logs/Plotting_DNA/plot_HeterozygosityPi_InPop/plot_HeterozygosityPi_InPop.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_HeterozygosityPi_InPop/plot_HeterozygosityPi_InPop.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pHetSim",
		threads = 1,
		mem = 50000,
		AtlSamples=config["AtlSamples"],
		MedSamples=config["MedSamples"]
	shell:
		"./scripts/Plotting_DNA/plot_HeterozygosityPi_InPop.R {input.HetAll} \"{input.PiTObs}\" {input.SamplesOrderInVCF} {output.PDF} \"{params.AtlSamples}\" \"{params.MedSamples}\" {input.Rconfig} > {log.out} 2> {log.err}"

# in dev!
rule plot_PiPerFixedWindows:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		PiWObs=expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, BED="FixedLengthRegions", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		BED = expand(rules.Windows_AccumulatedCallableSize_BEDs.output.WindowsBED, windowsize=config["windowsize"]),
		NonExtraCallableRegions=rules.r7_join_extra_callable_regions.output.NonExtraCallableRegions,
		ChrLen = "data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt",
	output:
		PDF = "results/Plotting_DNA/plot_PiPerFixedWindows/plot_PiPerFixedWindows.pdf"
	log:
		err = "logs/Plotting_DNA/plot_PiPerFixedWindows/plot_PiPerFixedWindows.err",
		out = "logs/Plotting_DNA/plot_PiPerFixedWindows/plot_PiPerFixedWindows.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_PiPerFixedWindows/plot_PiPerFixedWindows.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pPiWind",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_PiPerFixedWindows.R \"{input.PiWObs}\" {input.BED} {input.NonExtraCallableRegions} {input.ChrLen} {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"


# add significance test
rule plot_SharedPrivate_Populations:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		ObsShPriv = expand(rules.SharedPrivatePopulation_BEDs.output.Numbers, ObsOrBoots="Observed"),
		RandShPriv = expand(rules.SharedPrivatePopulation_BEDs.output.Numbers, ObsOrBoots=["Rand_%d" % i for i in range(config["BootsRandPop"])]),
	output:
		PDF = "results/Plotting_DNA/plot_SharedPrivate_Populations/plot_SharedPrivate_Populations.pdf"
	log:
		err = "logs/Plotting_DNA/plot_SharedPrivate_Populations/plot_SharedPrivate_Populations.err",
		out = "logs/Plotting_DNA/plot_SharedPrivate_Populations/plot_SharedPrivate_Populations.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_SharedPrivate_Populations/plot_SharedPrivate_Populations.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pShPriv",
		threads = 1,
		mem = 2000
	shell:
		"./scripts/Plotting_DNA/plot_SharedPrivate_Populations.R \"{input.ObsShPriv}\" \"{input.RandShPriv}\" {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"

# debug!
rule plot_PerSite_VarTypes_FuncRegions_Populations:
	'''
	Plot:
	- Site Frequency Spectrum
	- Number of alleles histogram
	'''
	input:
		Rconfig = config["Rconfig"],
		FreqObs=expand(rules.FrequencyPerSite_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples="AllSamples", ObsExp="Observed_Data"),
		TypesOfVariants=expand(rules.VariantType_BEDs_PerSubsetOfSamples.output.PerSite_VariantType, GroupSamples="AllSamples"),
		FunctionalFeatOfVariants=rules.FunctionalFeatures_BEDs.output.PerSite_FeatureType,
		#SharedPrivatePop=expand(rules.SharedPrivatePopulation_BEDs.output.AllVarBothPop, ObsOrBoots="Observed"),
		SynNonSyn=rules.SynonymousNonSynonymous_BEDs.output.AllVar_SynNonSyn,
		SpanFuncRegions=rules.get_FunctionalRegionsCallableSpan.output.out
	output:
		PDF = "results/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations/plot_PerSite_VarTypes_FuncRegions_Populations.pdf"
	log:
		err = "logs/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations/plot_PerSite_VarTypes_FuncRegions_Populations.err",
		out = "logs/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations/plot_PerSite_VarTypes_FuncRegions_Populations.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations/plot_PerSite_VarTypes_FuncRegions_Populations.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "pSite",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_PerSite_VarTypes_FuncRegions_Populations.R {input.FreqObs} {input.TypesOfVariants} {input.FunctionalFeatOfVariants} {input.SharedPrivatePop} {input.SynNonSyn} {input.SpanFuncRegions} {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"

# in dev!
rule plot_PerRegion_Diversity:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		PiTObs=expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, BED="Callable", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiTObsExons=expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, BED="Exons", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiTObsIntrons=expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, BED="Introns", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiTObsPromoters=expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, BED="Promoters", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiTObsIntergenic=expand(rules.Pi_Total_inBEDregions_PerSubsetOfSamples.output.out, BED="Intergenic", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiPerRObsExons=expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, BED="Exons", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiPerRObsIntrons=expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, BED="Introns", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiPerRObsPromoters=expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, BED="Promoters", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
		PiPerRObsIntergenic=expand(rules.Pi_InEachRegion_inBEDregions_PerSubsetOfSamples.output.out, BED="Intergenic", GroupSamples=["AllSamples", "AtlSamples", "MedSamples"], ObsExp="Observed_Data"),
	output:
		PDF = "results/Plotting_DNA/plot_PerRegion_Diversity/plot_PerRegion_Diversity.pdf"
	log:
		err = "logs/Plotting_DNA/plot_PerRegion_Diversity/plot_PerRegion_Diversity.err",
		out = "logs/Plotting_DNA/plot_PerRegion_Diversity/plot_PerRegion_Diversity.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_PerRegion_Diversity/plot_PerRegion_Diversity.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "pPerRegion_Div",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_PerRegion_Diversity.R \"{input.PiTObs}\" \"{input.PiTObsExons}\" \"{input.PiTObsIntrons}\" \"{input.PiTObsPromoters}\" \"{input.PiTObsIntergenic}\" \"{input.PiPerRObsExons}\" \"{input.PiPerRObsIntrons}\" \"{input.PiPerRObsPromoters}\" \"{input.PiPerRObsIntergenic}\" {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"


rule plot_GOEnrichment:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		FeatMatrices = expand("results/VariantAnalysis_DNA/PerFeature_StatsMatrix_PerChr/PerFeatureStats.{chr}.tab.gz", chr=config["bralan3chrs"]),
		GOAnnotations = "data/BraLan3Gene_2_GOterm.txt"
	output:
		PDF = "results/Plotting_DNA/plot_GOEnrichment/plot_GOEnrichment.pdf"
	log:
		err = "logs/Plotting_DNA/plot_GOEnrichment/plot_GOEnrichment.err",
		out = "logs/Plotting_DNA/plot_GOEnrichment/plot_GOEnrichment.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_GOEnrichment/plot_GOEnrichment.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "GOEnrich",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_GOEnrichment.R \"{input.FeatMatrices}\" {input.GOAnnotations} {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"
####

rule plot_HeterozygosityDistribution_vsSimulations:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		ObsHet = expand(rules.Heterozygosity_InEachRegion_inBEDregions_PerSample.output.out, BED="SelectedWindowsSimChrSize", GroupSamples="AllSamples", ObsExp="Observed_Data"),
		SimHet = expand(rules.Heterozygosity_Total_inBEDregions_PerSample.output.out, BED=["Simulated_%d_FiniteAlleles" % i for i in config["SimulationsReplicates"]], GroupSamples="SimSamples", ObsExp=["Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_%s" % i for i in config["SimulationsRangeMuNe"]]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt",
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
	output:
		PDF = "results/Plotting_DNA/plot_HeterozygosityDistribution_vsSimulations/plot_HeterozygosityDistribution_vsSimulations.pdf"
	log:
		err = "logs/Plotting_DNA/plot_HeterozygosityDistribution_vsSimulations/plot_HeterozygosityDistribution_vsSimulations.err",
		out = "logs/Plotting_DNA/plot_HeterozygosityDistribution_vsSimulations/plot_HeterozygosityDistribution_vsSimulations.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_HeterozygosityDistribution_vsSimulations/plot_HeterozygosityDistribution_vsSimulations.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pHetSim",
		threads = 1,
		mem = 50000
	shell:
		"./scripts/Plotting_DNA/plot_HeterozygosityDistribution_vsSimulations.R {input.ObsHet} \"{input.SimHet}\" {input.SamplesOrderInVCF} {input.Metadata} {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"


rule plot_SFS_NumAlleles_vsSimulations:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		ObsSFSAtl = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutMajorFreq, BED="SelectedWindowsSimChrSize", GroupSamples="AtlSamples", ObsExp="Observed_SNPs"),
		ObsSFSMed = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutMajorFreq, BED="SelectedWindowsSimChrSize", GroupSamples="MedSamples", ObsExp="Observed_SNPs"),
		SimSFS = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutMajorFreq, BED=["Simulated_%d_FiniteAlleles" % i for i in config["SimulationsReplicates"]], GroupSamples="SimSamples", ObsExp=["Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_%s" % i for i in config["SimulationsRangeMuNe"]]),
		ObsNAllelesAtl = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutNumAlleles, BED="SelectedWindowsSimChrSize", GroupSamples="AtlSamples", ObsExp="Observed_SNPs"),
		ObsNAllelesMed = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutNumAlleles, BED="SelectedWindowsSimChrSize", GroupSamples="MedSamples", ObsExp="Observed_SNPs"),
		SimNAlleles = expand(rules.Calc_NumAlleles_MajorFreq_InEachRegion_inBEDregions.output.OutNumAlleles, BED=["Simulated_%d_FiniteAlleles" % i for i in config["SimulationsReplicates"]], GroupSamples="SimSamples", ObsExp=["Expected_AsInSimulations/SLiM_simulation_singlePopulation/SimParam_1000_%s" % i for i in config["SimulationsRangeMuNe"]]),
		SamplesOrderInVCF = "metadata/SamplesOrderInVCF.chr19.txt",
		Metadata = "metadata/DNA_ReadGroups_Metadata_HiSeq4000_NovaSeq6000.txt",
	output:
		PDF = "results/Plotting_DNA/plot_SFS_NumAlleles_vsSimulations/plot_SFS_NumAlleles_vsSimulations.pdf"
	log:
		err = "logs/Plotting_DNA/plot_SFS_NumAlleles_vsSimulations/plot_SFS_NumAlleles_vsSimulations.err",
		out = "logs/Plotting_DNA/plot_SFS_NumAlleles_vsSimulations/plot_SFS_NumAlleles_vsSimulations.out"
	benchmark:
		"benchmarks/Plotting_DNA/plot_SFS_NumAlleles_vsSimulations/plot_SFS_NumAlleles_vsSimulations.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '1:00:00',
		name = "pSFSSim",
		threads = 1,
		mem = 50000,
		SimulationsSampleSize = config["SimulationsSampleSize"]
	shell:
		"./scripts/Plotting_DNA/plot_SFS_NumAlleles_vsSimulations.R {input.ObsSFSAtl} {input.ObsSFSMed} \"{input.SimSFS}\" {input.ObsNAllelesAtl} {input.ObsNAllelesMed} \"{input.SimNAlleles}\" {input.SamplesOrderInVCF} {input.Metadata} {output.PDF} {params.SimulationsSampleSize} {input.Rconfig} > {log.out} 2> {log.err}"






























