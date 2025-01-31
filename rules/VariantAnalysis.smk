
################################################
## Prepare basic files for variant analysis
################################################

# Read VCF and build genotype, numericmatrix and 0hom1het files
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
		chrGENOTYPE = "results/VariantAnalysis/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.genotype.bed.gz",
		chrGENOTYPEMatrix = "results/VariantAnalysis/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.genotype.numericmatrix.gz",
		HeterozygosityFile = "results/VariantAnalysis/Prepare_basic_analysis_files_PerChr/ShortVariants.filtered.{chr}.0hom1het.bed.gz"
	log:
		err = "logs/VariantAnalysis/Prepare_basic_analysis_files_PerChr/Prepare_basic_analysis_files_PerChr.{chr}.err",
		out = "logs/VariantAnalysis/Prepare_basic_analysis_files_PerChr/Prepare_basic_analysis_files_PerChr.{chr}.out"
	benchmark:
		"benchmarks/VariantAnalysis/Prepare_basic_analysis_files_PerChr/Prepare_basic_analysis_files_PerChr.{chr}.txt"
	conda:
		'../envs/VariantAnalysis.yaml'
	params:
		time = '5:00:00',
		name = "{chr}Genotype",
		threads = 1,
		mem = 25000
	shell:
		"./scripts/VariantAnalysis/Prepare_basic_analysis_files_PerChr.sh {input.chrVCF} {output.chrGENOTYPE} > {log.out} 2> {log.err}"
#
#
#








##################################################
# Rules from plotting_DNA





##########################
# Functions
# Provide bed file name according to the wildcard BED
def get_bed_file_name(wildcards):
	if wildcards.BED == "Callable":
		return rules.r7_join_extra_callable_regions.output.ExtraCallableRegions
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
	elif wildcards.BED == "PrivateA":
		return rules.SharedPrivatePopulation_BEDs.output.PrivateA
	elif wildcards.BED == "PrivateM":
		return rules.SharedPrivatePopulation_BEDs.output.PrivateM
	elif wildcards.BED == "Shared":
		return rules.SharedPrivatePopulation_BEDs.output.Shared
	elif wildcards.BED == "WindowsSimChrSize":
		return rules.SimualtedChrSizeWindows_BEDs.output.WindowsBED
	elif wildcards.BED == "SelectedWindowsSimChrSize":
		return rules.ChooseTreeRegions.output.SelectedBED
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
#
#
#
#
#


# Compute total average pairwise differences in a set of BED regions and a group of samples
rule Pi_Total_inBEDregions_PerSubsetOfSamples:
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

# add printing last region!
# Compute average pairwise differences for each one of a set of BED regions, for a group of samples
rule Pi_InEachRegion_inBEDregions_PerSubsetOfSamples:
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

# Build files with a list of samples corresponding to the wildcard RandObsSim (including simulated samples and randomized samples)
rule Prepare_PopulationSamples_Files:
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

############################
#### Rules generating BED files
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
		PerSite_FeatureType = "results/Plotting_DNA/FunctionalFeatures_BEDs/PerSite_FeatureType.txt.gz"
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

# *
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

#


