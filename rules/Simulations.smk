
def get_prevSimInput_ifExists(wildcards):
	i = int(wildcards.iteration)
	if i > 0:
		return "results/Simulations/SLiM_simulation_singlePopulation/SimParam_%s_%s_%s/SLiM_OutputFull_%s_%s_%s_%s_%d.binary" % (wildcards.chrsize, wildcards.mu, wildcards.Ne, wildcards.chrsize, wildcards.mu, wildcards.Ne, wildcards.Replica, i-1)
	elif i == 0:
		return "scripts/Simulations/SLiM_simulation_template_singlePopulation.slim"
	else:
		raise ValueError("Iteration numbers must be 0 or greater: received %s" % wildcards.iteration)

def get_timepoints(wildcards):
	return expand("{timepoint}", timepoint=config["SimulationsOutputTimepoints%s" % wildcards.iteration])

# To test!
rule SLiM_simulation_twoPopulations:
	'''
	'''
	input:
		PrevSimInput=get_prevSimInput_ifExists
	output:
		SlimSampleOutputP1="results/Simulations/SLiM_simulation_twoPopulations/SimParam_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}/SLiM_Sample_P1_SimParam_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.slim.gz",
		SlimSampleOutputP2="results/Simulations/SLiM_simulation_twoPopulations/SimParam_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}/SLiM_Sample_P2_SimParam_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.slim.gz",
		SLiMGeneralOutput="results/Simulations/SLiM_simulation_twoPopulations/SimParam_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}/SLiM_General_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.slim.gz",
		SLiMOutputFull="results/Simulations/SLiM_simulation_twoPopulations/SimParam_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}/SLiM_OutputFull_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.binary"
	wildcard_constraints:
		iteration="[0-99]"
	log:
		err = "logs/Simulations/SLiM_simulation_twoPopulations/SLiM_simulation_twoPopulations_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.err",
		out = "logs/Simulations/SLiM_simulation_twoPopulations/SLiM_simulation_twoPopulations_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.out"
	benchmark:
		"benchmarks/Simulations/SLiM_simulation_twoPopulations/SLiM_simulation_twoPopulations_{chrsize}_{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '50:00:00',
		name = "{mu}_{size1}_{size2}_{mig12}_{mig21}_{Replica}_{iteration}_2p",
		threads = 1,
		mem = 5000,
		timepoints = get_timepoints,
		SampleSize = expand("{SampleSize}", SampleSize=config["SimulationsSampleSize"]),
		ScriptTemplate="scripts/Simulations/SLiM_simulation_template_singlePopulation.slim"
	shell:
		"""
		outputP1={output.SLiMSampleOutputP1}
		outputP2={output.SLiMSampleOutputP2}
		outputG={output.SLiMGeneralOutput}
		outputF={output.SLiMOutputFull}
		for file in ${{outputP1%.gz}}* ${{outputP2%.gz}}* ${{outputG%.gz}}* ${{outputF%.gz}}* {log.out} {log.err}
		do
			if [[ -s ${{file}} ]]
			then
				rm ${{file}}
			fi
		done
		mkdir -p $(dirname ${{outputP1}}) 2>> {log.err}

		###############################
		# Prepare script for this SLiM run
		## Print initialize simulation
		Runscript=$(echo ${{output%.gz}} | sed 's/.slim/.slimsc/g' | sed 's/_Sample//g')
		cat {params.ScriptTemplate} | sed '/END initialize simulation/q' | sed '$ d' > ${{Runscript}}

		## Print first early
		echo {input.PrevSimInput} >> {log.out}
		if [[ {input.PrevSimInput} == "scripts/Simulations/SLiM_simulation_template_singlePopulation.slim" ]]
		then
			lineS=$(cat {params.ScriptTemplate} | grep -n 'START Create two subpopulations with the corresponding Ne' | cut -f1 -d':')
			lineP=$(cat {params.ScriptTemplate} | tail -n +${{lineS}} | grep -n 'END Create two subpopulations with the corresponding Ne' | cut -f1 -d':')
			cat {params.ScriptTemplate} | tail -n +${{lineS}} | head -${{lineP}} >> ${{Runscript}}
		else
			lineS=$(cat {params.ScriptTemplate} | grep -n 'START Read previous iteration results' | cut -f1 -d':')
			lineP=$(cat {params.ScriptTemplate} | tail -n +${{lineS}} | grep -n 'END Read previous iteration results' | cut -f1 -d':')
			cat {params.ScriptTemplate} | tail -n +${{lineS}} | head -${{lineP}} >> ${{Runscript}}
		fi
		## Print the per timepoint printing part
		lineP=$(cat {params.ScriptTemplate} | grep -n 'START Print on a range of time points' | cut -f1 -d':')
		lineE=$(cat {params.ScriptTemplate} | tail -n +${{lineP}} | grep -n 'END Print on a range of time points' | cut -f1 -d':')
		for timepoint in {params.timepoints}
		do
			cat {params.ScriptTemplate} | tail -n +${{lineP}} | head -${{lineE}} | sed "s/XXX/$timepoint/g" >> ${{Runscript}}
		done

		## Print fixed mutations at the end
		lineS=$(cat {params.ScriptTemplate} | grep -n 'START of print fixed muations' | cut -f1 -d':')
		lineP=$(cat {params.ScriptTemplate} | tail -n +${{lineS}} | grep -n 'END of print fixed muations' | cut -f1 -d':')
		cat {params.ScriptTemplate} | tail -n +${{lineS}} | head -${{lineP}} | sed "s/XXX/$timepoint/g" >> ${{Runscript}}

		###############################
		# Run SLiM
		slim -d "output1='${{output1%.gz}}'" -d "output2='${{output2%.gz}}'" -d "outputG='${{outputG%.gz}}'" -d "outputF='${{outputF%.gz}}'" -d "input='{input.PrevSimInput}'" -d SampleSize={params.SampleSize} -d chrsize={wildcards.chrsize} -d mu={wildcards.mu} -d size1={wildcards.size1} -d size2={wildcards.size2} -d mig12={wildcards.mig12} -d mig21={wildcards.mig21} ${{Runscript}}  >> {log.out} 2>> {log.err}
		gzip ${{outputP1%.gz}} >> {log.out} 2>> {log.err}
		gzip ${{outputP2%.gz}} >> {log.out} 2>> {log.err}
		gzip ${{outputG%.gz}} >> {log.out} 2>> {log.err}
		"""

rule SLiM_simulation_singlePopulation:
	'''
	'''
	input:
		PrevSimInput=get_prevSimInput_ifExists
	output:
		SLiMSampleOutput="results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{mu}_{Ne}/SLiM_Sample_{chrsize}_{mu}_{Ne}_{Replica}_{iteration}.slim.gz",
		SLiMGeneralOutput="results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{mu}_{Ne}/SLiM_General_{chrsize}_{mu}_{Ne}_{Replica}_{iteration}.slim.gz",
		SLiMOutputFull="results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{mu}_{Ne}/SLiM_OutputFull_{chrsize}_{mu}_{Ne}_{Replica}_{iteration}.binary"
	wildcard_constraints:
		iteration="[0-9]*"
	log:
		err = "logs/Simulations/SLiM_simulation_singlePopulation/SLiM_simulation_singlePopulation_{chrsize}_{mu}_{Ne}_{Replica}_{iteration}.err",
		out = "logs/Simulations/SLiM_simulation_singlePopulation/SLiM_simulation_singlePopulation_{chrsize}_{mu}_{Ne}_{Replica}_{iteration}.out"
	benchmark:
		"benchmarks/Simulations/SLiM_simulation_singlePopulation/SLiM_simulation_singlePopulation_{chrsize}_{mu}_{Ne}_{Replica}_{iteration}.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '50:00:00',
		name = "{mu}_{Ne}_{Replica}_1p",
		threads = 1,
		mem = 5000,
		timepoints = get_timepoints,
		SampleSize = expand("{SampleSize}", SampleSize=config["SimulationsSampleSize"]),
		ScriptTemplate="scripts/Simulations/SLiM_simulation_template_singlePopulation.slim"
	shell:
		"""
		output={output.SLiMSampleOutput}
		outputG={output.SLiMGeneralOutput}
		outputF={output.SLiMOutputFull}
		for file in ${{output%.gz}}* ${{outputG%.gz}}* ${{outputF%.gz}}* {log.out} {log.err}
		do
			if [[ -s ${{file}} ]]
			then
				rm ${{file}}
			fi
		done
		mkdir -p $(dirname ${{output}}) 2>> {log.err}

		###############################
		# Prepare script for this SLiM run
		## Print initialize simulation
		Runscript=$(echo ${{output%.gz}} | sed 's/.slim/.slimsc/g' | sed 's/_Sample//g')
		cat {params.ScriptTemplate} | sed '/END initialize simulation/q' | sed '$ d' > ${{Runscript}}

		## Print first early
		echo {input.PrevSimInput} >> {log.out}
		if [[ {input.PrevSimInput} == "scripts/Simulations/SLiM_simulation_template_singlePopulation.slim" ]]
		then
			lineS=$(cat {params.ScriptTemplate} | grep -n 'START Create one subpopulation with the corresponding Ne' | cut -f1 -d':')
			lineP=$(cat {params.ScriptTemplate} | tail -n +${{lineS}} | grep -n 'END Create one subpopulation with the corresponding Ne' | cut -f1 -d':')
			cat {params.ScriptTemplate} | tail -n +${{lineS}} | head -${{lineP}} >> ${{Runscript}}
		else
			lineS=$(cat {params.ScriptTemplate} | grep -n 'START Read previous iteration results' | cut -f1 -d':')
			lineP=$(cat {params.ScriptTemplate} | tail -n +${{lineS}} | grep -n 'END Read previous iteration results' | cut -f1 -d':')
			cat {params.ScriptTemplate} | tail -n +${{lineS}} | head -${{lineP}} >> ${{Runscript}}
		fi
		## Print the per timepoint printing part
		lineP=$(cat {params.ScriptTemplate} | grep -n 'START Print on a range of time points' | cut -f1 -d':')
		lineE=$(cat {params.ScriptTemplate} | tail -n +${{lineP}} | grep -n 'END Print on a range of time points' | cut -f1 -d':')
		for timepoint in {params.timepoints}
		do
			cat {params.ScriptTemplate} | tail -n +${{lineP}} | head -${{lineE}} | sed "s/XXX/$timepoint/g" >> ${{Runscript}}
		done

		## Print fixed mutations at the end
		lineS=$(cat {params.ScriptTemplate} | grep -n 'START of print fixed muations' | cut -f1 -d':')
		lineP=$(cat {params.ScriptTemplate} | tail -n +${{lineS}} | grep -n 'END of print fixed muations' | cut -f1 -d':')
		cat {params.ScriptTemplate} | tail -n +${{lineS}} | head -${{lineP}} | sed "s/XXX/$timepoint/g" >> ${{Runscript}}

		###############################
		# Run SLiM
		slim -d "output='${{output%.gz}}'" -d "outputG='${{outputG%.gz}}'" -d "outputF='${{outputF%.gz}}'" -d "input='{input.PrevSimInput}'" -d SampleSize={params.SampleSize} -d chrsize={wildcards.chrsize} -d mu={wildcards.mu} -d Ne={wildcards.Ne} ${{Runscript}}  >> {log.out} 2>> {log.err}
		gzip ${{output%.gz}} >> {log.out} 2>> {log.err}
		gzip ${{outputG%.gz}} >> {log.out} 2>> {log.err}
		"""

rule SLiM_output_to_chrgenotypes_InfiniteAlleles:
	'''
	Infinite Alleles
	SimParam = {chrsize}_{mu}_{Ne}_{Replica}
	'''
	input:
		SLiMOutputs=expand("results/Simulations/{{SimFolder}}/SLiM_Sample_{{SimParam}}_{it}.slim.gz", it=config["SimulationsIterations"])
	output:
		ChrGenotypesFile="results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_InfiniteAlleles.chrgenotypes.gz"
	log:
		err = "logs/Simulations/{SimFolder}/SLiM_output_to_chrgenotypes_{SimParam}_InfiniteAlleles.err",
		out = "logs/Simulations/{SimFolder}/SLiM_output_to_chrgenotypes_{SimParam}_InfiniteAlleles.out"
	benchmark:
		"benchmarks/Simulations/{SimFolder}/SLiM_output_to_chrgenotypes_{SimParam}_InfiniteAlleles.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '3:00:00',
		name = "2GenoIA{SimParam}",
		threads = 1,
		mem = 2000,
	shell:
		"scripts/Simulations/SLiM_output_to_chrgenotypes_InfiniteAlleles.py \"{input.SLiMOutputs}\" {output.ChrGenotypesFile} > {log.out} 2> {log.err};"

rule SLiM_output_to_chrgenotypes_FiniteAlleles:
	'''
	Finite Alleles (accounting for back mutations)

	For single population simulatons:
		SimParam = {chrsize}_{mu}_{Ne}_{Replica}
	'''
	input:
		SLiMOutputs=expand("results/Simulations/{{SimFolder}}/SLiM_Sample_{{SimParam}}_{it}.slim.gz", it=config["SimulationsIterations"])
	output:
		ChrGenotypesFile="results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_FiniteAlleles.chrgenotypes.gz"
	log:
		err = "logs/Simulations/{SimFolder}/SLiM_output_to_chrgenotypes_{SimParam}_FiniteAlleles.err",
		out = "logs/Simulations/{SimFolder}/SLiM_output_to_chrgenotypes_{SimParam}_FiniteAlleles.out"
	benchmark:
		"benchmarks/Simulations/{SimFolder}/SLiM_output_to_chrgenotypes_{SimParam}_FiniteAlleles.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '3:00:00',
		name = "2GenoFA{SimParam}",
		threads = 1,
		mem = 2000,
	shell:
		"scripts/Simulations/SLiM_output_to_chrgenotypes_FiniteAlleles.py \"{input.SLiMOutputs}\" {output.ChrGenotypesFile} > {log.out} 2> {log.err};"

def get_chrgenotype_file_name(wildcards):
	if "InfiniteAlleles" in wildcards.FiniteInfiniteAlleles:
		return rules.SLiM_output_to_chrgenotypes_InfiniteAlleles.output.ChrGenotypesFile
	elif "FiniteAlleles" in wildcards.FiniteInfiniteAlleles:
		return rules.SLiM_output_to_chrgenotypes_FiniteAlleles.output.ChrGenotypesFile
	else:
		raise ValueError("Accepted values: InfiniteAlleles or FiniteAlleles; not %s" % wildcards.FiniteInfiniteAlleles)

rule Formating_SimulationFiles:
	'''
	Only last generation
	'''
	input:
		ChrGenotypesFile=get_chrgenotype_file_name
	output:
		HetFile = "results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_{FiniteInfiniteAlleles}.0hom1het.bed.gz",
		GenoFile = "results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_{FiniteInfiniteAlleles}.genotype.bed.gz",
		BED = "results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_{FiniteInfiniteAlleles}.bed.gz"
	log:
		err = "logs/Simulations/{SimFolder}/Formating_SimulationFiles_{SimParam}_{FiniteInfiniteAlleles}.err",
		out = "logs/Simulations/{SimFolder}/Formating_SimulationFiles_{SimParam}_{FiniteInfiniteAlleles}.out"
	benchmark:
		"benchmarks/Simulations/{SimFolder}/Formating_SimulationFiles_{SimParam}_{FiniteInfiniteAlleles}.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '3:00:00',
		name = "FormSim",
		threads = 1,
		mem = 2000,
		chrsize = expand("{chrsize}", chrsize=config["SimulationsChrSize"])
	shell:
		"scripts/Simulations/Formating_SimulationFiles.py {input.ChrGenotypesFile} {output.HetFile} {output.GenoFile} {output.BED} {params.chrsize} > {log.out} 2> {log.err};"

rule Calc_Pi_from_GenotypesFile:
	'''
	For single population simulatons:
		SimParam = {chrsize}_{mu}_{Ne}_{Replica}_{InfiniteAlleles/FiniteAlleles}
	'''
	input:
		ChrGenotypesFile=get_chrgenotype_file_name
	output:
		PIfile="results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_{FiniteInfiniteAlleles}.pi.gz"
	log:
		err = "logs/Simulations/{SimFolder}/Calc_Pi_from_GenotypesFile_{SimParam}_{FiniteInfiniteAlleles}.err",
		out = "logs/Simulations/{SimFolder}/Calc_Pi_from_GenotypesFile_{SimParam}_{FiniteInfiniteAlleles}.out"
	benchmark:
		"benchmarks/Simulations/{SimFolder}/Calc_Pi_from_GenotypesFile_{SimParam}_{FiniteInfiniteAlleles}.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '3:00:00',
		name = "piSim",
		threads = 1,
		mem = 2000,
		chrsize = expand("{chrsize}", chrsize=config["SimulationsChrSize"])
	shell:
		"scripts/Simulations/Calc_Pi_from_GenotypesFile.py {input.ChrGenotypesFile} {output.PIfile} {params.chrsize} > {log.out} 2> {log.err};"

rule Calc_Het_from_GenotypesFile:
	'''
	For single population simulatons:
		SimParam = {chrsize}_{mu}_{Ne}_{Replica}_{InfiniteAlleles/FiniteAlleles}
	'''
	input:
		ChrGenotypesFile=get_chrgenotype_file_name
	output:
		Hetfile="results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_{FiniteInfiniteAlleles}.het.gz"
	log:
		err = "logs/Simulations/{SimFolder}/Calc_Het_from_GenotypesFile_{SimParam}_{FiniteInfiniteAlleles}.err",
		out = "logs/Simulations/{SimFolder}/Calc_Het_from_GenotypesFile_{SimParam}_{FiniteInfiniteAlleles}.out"
	benchmark:
		"benchmarks/Simulations/{SimFolder}/Calc_Het_from_GenotypesFile_{SimParam}_{FiniteInfiniteAlleles}.txt"
	conda:
		'../envs/Simulations.yaml'
	params:
		time = '3:00:00',
		name = "hetSim",
		threads = 1,
		mem = 2000,
		chrsize = expand("{chrsize}", chrsize=config["SimulationsChrSize"])
	shell:
		"scripts/Simulations/Calc_Het_from_GenotypesFile.py {input.ChrGenotypesFile} {output.Hetfile} {params.chrsize} > {log.out} 2> {log.err};"

rule plot_Pi_RangeNe_RangeMu_singlePopulation:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		PIfilesI=expand("results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{muNe}/SLiM_Sample_{chrsize}_{muNe}_{Replica}_InfiniteAlleles.pi.gz", chrsize=config["SimulationsChrSize"], muNe=config["SimulationsRangeMuNe"], Replica=config["SimulationsReplicates"]),
		PIfilesF=expand("results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{muNe}/SLiM_Sample_{chrsize}_{muNe}_{Replica}_FiniteAlleles.pi.gz", chrsize=config["SimulationsChrSize"], muNe=config["SimulationsRangeMuNe"], Replica=config["SimulationsReplicates"]),
	output:
		PDF="results/Simulations/plot_Pi_RangeNe_RangeMu_singlePopulation/plot_Pi_RangeNe_RangeMu_singlePopulation.pdf"
	log:
		err = "logs/Simulations/plot_Pi_RangeNe_RangeMu_singlePopulation/plot_Pi_RangeNe_RangeMu_singlePopulation.err",
		out = "logs/Simulations/plot_Pi_RangeNe_RangeMu_singlePopulation/plot_Pi_RangeNe_RangeMu_singlePopulation.out"
	benchmark:
		"benchmarks/Simulations/plot_Pi_RangeNe_RangeMu_singlePopulation/plot_Pi_RangeNe_RangeMu_singlePopulation.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "pPi1",
		threads = 1,
		mem = 2000
	shell:
		"scripts/Simulations/plot_Pi_RangeNe_RangeMu_singlePopulation.R \"{input.PIfilesI} {input.PIfilesF}\" {output.PDF} {input.Rconfig} > {log.out} 2> {log.err}"

rule plot_Het_RangeNe_RangeMu_singlePopulation:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		HetfilesI=expand("results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{muNe}/SLiM_Sample_{chrsize}_{muNe}_{Replica}_InfiniteAlleles.het.gz", chrsize=config["SimulationsChrSize"], muNe=config["SimulationsRangeMuNe"], Replica=config["SimulationsReplicates"]),
		HetfilesF=expand("results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{muNe}/SLiM_Sample_{chrsize}_{muNe}_{Replica}_FiniteAlleles.het.gz", chrsize=config["SimulationsChrSize"], muNe=config["SimulationsRangeMuNe"], Replica=config["SimulationsReplicates"]),
	output:
		PDF="results/Simulations/plot_Het_RangeNe_RangeMu_singlePopulation/plot_Het_RangeNe_RangeMu_singlePopulation.pdf"
	log:
		err = "logs/Simulations/plot_Het_RangeNe_RangeMu_singlePopulation/plot_Het_RangeNe_RangeMu_singlePopulation.err",
		out = "logs/Simulations/plot_Het_RangeNe_RangeMu_singlePopulation/plot_Het_RangeNe_RangeMu_singlePopulation.out"
	benchmark:
		"benchmarks/Simulations/plot_Het_RangeNe_RangeMu_singlePopulation/plot_Het_RangeNe_RangeMu_singlePopulation.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "pHet1",
		threads = 1,
		mem = 2000,
		samplesize = expand("{samplesize}", samplesize=config["SimulationsSampleSize"])
	shell:
		"scripts/Simulations/plot_Het_RangeNe_RangeMu_singlePopulation.R \"{input.HetfilesI} {input.HetfilesF}\" {output.PDF} {input.Rconfig} {params.samplesize} > {log.out} 2> {log.err}"


rule Build_DistanceTree_LastGeneration_from_GenotypesFile:
	'''
	'''
	input:
		ChrGenotypesFile=rules.SLiM_output_to_chrgenotypes_FiniteAlleles.output.ChrGenotypesFile
	output:
		tree = "results/Simulations/{SimFolder}/SLiM_Sample_{SimParam}_FiniteAlleles.treefile"
	log:
		err = "logs/Simulations/{SimFolder}/Build_DistanceTree_LastGeneration_from_GenotypesFile_{SimParam}_FiniteAlleles.err",
		out = "logs/Simulations/{SimFolder}/Build_DistanceTree_LastGeneration_from_GenotypesFile_{SimParam}_FiniteAlleles.out"
	benchmark:
		"benchmarks/Simulations/{SimFolder}/Build_DistanceTree_LastGeneration_from_GenotypesFile_{SimParam}_FiniteAlleles.txt"
	conda:
		'../envs/Phylogenetics.yaml'
	params:
		time = '02:00:00',
		name = "treeSim",
		threads = 1,
		mem = 2000
	shell:
		"scripts/Simulations/Build_DistanceTree_LastGeneration_from_GenotypesFile.sh {input.ChrGenotypesFile} {output.tree} {params.threads} > {log.out} 2> {log.err}"

rule plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation:
	'''
	'''
	input:
		Rconfig = config["Rconfig"],
		trees=expand("results/Simulations/SLiM_simulation_singlePopulation/SimParam_{chrsize}_{muNe}/SLiM_Sample_{chrsize}_{muNe}_{Replica}_FiniteAlleles.treefile", chrsize=config["SimulationsChrSize"], muNe=config["SimulationsRangeMuNe"], Replica=config["SimulationsReplicates"])
	output:
		BranchLengths="results/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation/ExternalVsInternalBranchLengths_RangeNe_RangeMu_singlePopulation.tab",
		PDF="results/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation.pdf"
	log:
		err = "logs/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation.err",
		out = "logs/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation.out"
	benchmark:
		"benchmarks/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation.txt"
	conda:
		'../envs/Plotting_DNA.yaml'
	params:
		time = '3:00:00',
		name = "pTree1",
		threads = 1,
		mem = 2000,
		muNe = expand("{muNe}", muNe=config["SimulationsRangeMuNe"])
	shell:
		"""
		scripts/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation.sh \"{input.trees}\" {output.BranchLengths} > {log.out} 2> {log.err}
		scripts/Simulations/plot_DistanceTree_ExternalVsInternalBranchLengthsRatio_RangeNe_RangeMu_singlePopulation.R {output.BranchLengths} {output.PDF} {input.Rconfig} >> {log.out} 2>> {log.err}
		"""

















