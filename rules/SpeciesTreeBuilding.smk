################################################
## 121 ortholog analysis to build species tree
################################################


rule DownloadProteinSeq_from_Ensembl:
	'''
	Given a list of species and their genome reference version (in config file), download their protein sequences from ENSEMBL
	'''
	output:
		ProteinSequencesFiles = expand("data/Sequences/{version}_proteins.fa.gz", version=config["SpTreeGenomeVersions"]),
		TranscriptSequencesFiles = expand("data/Sequences/{version}_transcript.fa.gz", version=config["SpTreeGenomeVersions"])
	log:
		err = "logs/SpeciesTreeBuilding/DownloadProteinSeq_from_Ensembl/DownloadProteinSeq_from_Ensembl.err",
		out = "logs/SpeciesTreeBuilding/DownloadProteinSeq_from_Ensembl/DownloadProteinSeq_from_Ensembl.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/DownloadProteinSeq_from_Ensembl/DownloadProteinSeq_from_Ensembl.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '50:00:00',
		name = "DownloadProt",
		threads = 40,
		mem = 30000,
		species = expand("{sp}", sp=config["SpTreeSpecies"]),
		GenomeVersions = expand("{sp}", sp=config["SpTreeGenomeVersions"]),
		EnsemblRelease = expand("{rel}", rel=config["EnsemblRelease"]),
		EnsMetazoaRelease = expand("{rel}", rel=config["EnsMetazoaRelease"])
	shell:
		"./scripts/SpeciesTreeBuilding/DownloadProteinSeq_from_Ensembl.sh \"{output.ProteinSequencesFiles}\" \"{output.TranscriptSequencesFiles}\" \"{params.species}\" \"{params.GenomeVersions}\" {params.EnsemblRelease} {params.EnsMetazoaRelease} > {log.out} 2> {log.err}"

rule Filter_ProteinSequences_PerSpecies:
	'''
	Filter peptides of less than minlen in length
	'''
	input:
		ProteinSequences = "data/Sequences/{version}_proteins.fa.gz",
		TranscriptSequences = "data/Sequences/{version}_transcript.fa.gz"
	output:
		FilteredProteinSequences = "results/SpeciesTreeBuilding/Filter_ProteinSequences_PerSpecies/{version}_proteins.fa.gz",
		FilteredTranscriptSequences = "results/SpeciesTreeBuilding/Filter_ProteinSequences_PerSpecies/{version}_transcript.fa.gz"
	log:
		err = "logs/SpeciesTreeBuilding/Filter_ProteinSequences_PerSpecies/Filter_ProteinSequences_{version}.err",
		out = "logs/SpeciesTreeBuilding/Filter_ProteinSequences_PerSpecies/Filter_ProteinSequences_{version}.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/Filter_ProteinSequences_PerSpecies/Filter_ProteinSequences_{version}.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '00:20:00',
		name = "FilterSeq{version}",
		threads = 1,
		mem = 1000
	shell:
		"./scripts/SpeciesTreeBuilding/Filter_ProteinSequences_PerSpecies.sh {input.ProteinSequences} {input.TranscriptSequences} {output.FilteredProteinSequences} {output.FilteredTranscriptSequences} {wildcards.version} > {log.out} 2> {log.err}"

# rule to be run locally in the cluster in interactive (-X) mode
rule BuildOGs_withOrthoLoger:
	'''
	Given a list of species, build the orthologous groups with OrthoLoger
	Automates Orthologer setup and execution using run_loger_prompt.sh to handle interactive prompts.
	Author: Marina from Sagane code made with code from Giulia Campli 
	'''
	input:
		ProteinSequencesFiles = expand(rules.Filter_ProteinSequences_PerSpecies.output.FilteredProteinSequences, version=config["SpTreeGenomeVersions"]),
	output:
		OGs = "results/SpeciesTreeBuilding/BuildOGs_withOrthoLoger/OGs.tab.gz"
	log:
		err = "logs/SpeciesTreeBuilding/BuildOGs_withOrthoLoger/BuildOGs_withOrthoLoger.err",
		out = "logs/SpeciesTreeBuilding/BuildOGs_withOrthoLoger/BuildOGs_withOrthoLoger.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/BuildOGs_withOrthoLoger/BuildOGs_withOrthoLoger.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '10:10:00',
		name = "OGsOrthoLoger",
		threads = 40,
		mem = 40000,
		species = expand("{sp}", sp=config["SpTreeSpecies"]),
		orthologer_dir = "./scripts/SpeciesTreeBuilding/orthologer_3.0.5/ORTHOLOGER-3.0.5/bin",
		brhclus_dir = "./scripts/SpeciesTreeBuilding/orthologer_3.0.5/BRHCLUS-5.1.7/bin",
		script_prompt = "./scripts/SpeciesTreeBuilding/run_loger_prompt.sh"
	shell:
		"./scripts/SpeciesTreeBuilding/BuildOGs_withOrthoLoger.sh \"{input.ProteinSequencesFiles}\" {output.OGs} {params.orthologer_dir} {params.brhclus_dir} {params.script_prompt} \"{params.species}\" > {log.out} 2> {log.err}"

checkpoint Filter_nonOne2One_Orthogroups:
	'''
	Get One2One orthogroups only
	'''
	input:
		OGs = rules.BuildOGs_withOrthoLoger.output.OGs
	output:
		OGsDIR=directory("results/SpeciesTreeBuilding/Filter_nonOne2One_Orthogroups")
	log:
		err = "logs/SpeciesTreeBuilding/Filter_nonOne2One_Orthogroups/Filter_nonOne2One_Orthogroups.err",
		out = "logs/SpeciesTreeBuilding/Filter_nonOne2One_Orthogroups/Filter_nonOne2One_Orthogroups.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/Filter_nonOne2One_Orthogroups/Filter_nonOne2One_Orthogroups.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '20:00:00',
		name = "FilterOG",
		threads = 1,
		mem = 5000,
		species = expand("{sp}", sp=config["SpTreeSpecies"]),
		SpeciesGenePrefix = expand("{sp}", sp=config["SpTreeGenePrefix"])
	shell:
		"./scripts/SpeciesTreeBuilding/Filter_nonOne2One_Orthogroups.sh {input.OGs} {output.OGsDIR} \"{params.species}\" \"{params.SpeciesGenePrefix}\" > {log.out} 2> {log.err}"

rule ExtractProtSeq_MSAwMAFFT_perOrthogroup:
	'''
	Per orthogroup:
	- Get the protein sequences of all genes
	- Do a multiple sequence alignment with MAFFT
	- Trim the alignment with ClipKIT
	'''
	input:
		ProteinSequencesFiles = expand(rules.Filter_ProteinSequences_PerSpecies.output.FilteredProteinSequences, version=config["SpTreeGenomeVersions"]),
		TranscriptSequencesFiles = expand(rules.Filter_ProteinSequences_PerSpecies.output.FilteredTranscriptSequences, version=config["SpTreeGenomeVersions"]),
		OGGeneList = "results/SpeciesTreeBuilding/Filter_nonOne2One_Orthogroups/{OG}.tbl.gz"
	output:
		OG_PSeqs = "results/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/{OG}_Pseq.fa.gz",
		MSA_P="results/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/{OG}_aln_Pseq.fa.gz",
		TrimmedMSA_P="results/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/{OG}_aln_trimmed_Pseq.fa.gz",
		OG_TSeqs = "results/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/{OG}_Tseq.fa.gz",
		MSA_T="results/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/{OG}_aln_Tseq.fa.gz",
		TrimmedMSA_T="results/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/{OG}_aln_trimmed_Tseq.fa.gz"
	log:
		err = "logs/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/ExtractProtSeq_MSAwMAFFT_{OG}.err",
		out = "logs/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/ExtractProtSeq_MSAwMAFFT_{OG}.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup/ExtractProtSeq_MSAwMAFFT_{OG}.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '05:00:00',
		name = "MSA{OG}",
		threads = 1,
		mem = 1000,
		species = expand("{sp}", sp=config["SpTreeSpecies"])
	shell:
		"./scripts/SpeciesTreeBuilding/ExtractProtSeq_MSAwMAFFT_perOrthogroup.sh \"{input.ProteinSequencesFiles}\" \"{input.TranscriptSequencesFiles}\" {input.OGGeneList} {output.OG_PSeqs}  {output.OG_TSeqs} {output.MSA_P} {output.MSA_T} {output.TrimmedMSA_P} {output.TrimmedMSA_T} \"{params.species}\" > {log.out} 2> {log.err}"

def get_OG_file_names(wildcards):
	'''
	Aggregate the file names of the OGs after filtering and MSA
	'''
	checkpoint_output = checkpoints.Filter_nonOne2One_Orthogroups.get(**wildcards).output.OGsDIR
	if wildcards.MSAType == "Protein":
		return expand(rules.ExtractProtSeq_MSAwMAFFT_perOrthogroup.output.TrimmedMSA_P, OG=glob_wildcards(os.path.join(checkpoint_output, '{OG}.tbl.gz')).OG)
	elif wildcards.MSAType == "Transcript":
		return expand(rules.ExtractProtSeq_MSAwMAFFT_perOrthogroup.output.TrimmedMSA_T, OG=glob_wildcards(os.path.join(checkpoint_output, '{OG}.tbl.gz')).OG)
	else:
		raise ValueError("Unknown values for MSAType in get_OG_file_names: %s" % (wildcards.MSAType))

rule Concatenate_MSA_One2OneOR:
	input:
		OG_MSAs=get_OG_file_names
	output:
		ConcatenatedMSA="results/SpeciesTreeBuilding/Concatenate_MSA_One2OneOR/Concatenate_{MSAType}_MSA_One2OneOR.fa.gz"
	log:
		err = "logs/SpeciesTreeBuilding/Concatenate_MSA_One2OneOR/Concatenate_{MSAType}_MSA_One2OneOR.err",
		out = "logs/SpeciesTreeBuilding/Concatenate_MSA_One2OneOR/Concatenate_{MSAType}_MSA_One2OneOR.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/Concatenate_MSA_One2OneOR/Concatenate_{MSAType}_MSA_One2OneOR.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '20:00:00',
		name = "ConcatMSA",
		threads = 1,
		mem = 5000,
		species = expand("{sp}", sp=config["SpTreeSpecies"]),
		SpeciesGenePrefix = expand("{sp}", sp=config["SpTreeGenePrefix"])
	shell:
		"./scripts/SpeciesTreeBuilding/Concatenate_MSA_One2OneOR.sh \"{input.OG_MSAs}\" {output.ConcatenatedMSA} \"{params.species}\" {wildcards.MSAType} > {log.out} 2> {log.err}"

rule Build_Tree_Species:
	'''
	Using the protein concatenated alignment
	'''
	input:
		Rconfig = config["Rconfig"],
		ConcatenatedMSA = expand(rules.Concatenate_MSA_One2OneOR.output.ConcatenatedMSA, MSAType="Protein")
	output:
		tree = "results/SpeciesTreeBuilding/Build_Tree_Species/Build_Tree_Species.treefile"
	log:
		err = "logs/SpeciesTreeBuilding/Build_Tree_Species/Build_Tree_Species.err",
		out = "logs/SpeciesTreeBuilding/Build_Tree_Species/Build_Tree_Species.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/Build_Tree_Species/Build_Tree_Species.txt"
	conda:
		'../envs/Phylogenetics.yaml'
	params:
		time = '15:00:00',
		name = "treeSp",
		threads = 10,
		mem = 100000,
		SpTreeConstrainTree = expand("{tree}", tree=config["SpTreeConstrainTree"])
	shell:
		"scripts/SpeciesTreeBuilding/Build_Tree_Species.sh {input.ConcatenatedMSA} {output.tree} \"{params.SpTreeConstrainTree}\" {params.threads} > {log.out} 2> {log.err}"

rule dNdS_SpeciesTree:
	'''
	Using the concatenation of all the nucleotide backtranslation of the aminoacid alignments
	'''
	input:
		ConcatenatedMSA_FASTA = expand(rules.Concatenate_MSA_One2OneOR.output.ConcatenatedMSA, MSAType="Transcript"),
	output:
		ConcatenatedMSA_PHYLIP = "results/SpeciesTreeBuilding/dNdS_SpeciesTree/ConcatenatedMSA.phy.gz",
		CodemlOut = "results/SpeciesTreeBuilding/dNdS_SpeciesTree/dNdS_SpeciesTree.out"
	log:
		err = "logs/SpeciesTreeBuilding/dNdS_SpeciesTree/dNdS_SpeciesTree.err",
		out = "logs/SpeciesTreeBuilding/dNdS_SpeciesTree/dNdS_SpeciesTree.out"
	benchmark:
		"benchmarks/SpeciesTreeBuilding/dNdS_SpeciesTree/dNdS_SpeciesTree.txt"
	conda:
		'../envs/SpeciesTreeBuilding.yaml'
	params:
		time = '15:00:00',
		name = "dNdSSp",
		threads = 10,
		mem = 100000,
		SpTreeConstrainTree = expand("{tree}", tree=config["SpTreeConstrainTree_unroted"])
	shell:
		"""
		scripts/SpeciesTreeBuilding/dNdS_SpeciesTree.sh {input.ConcatenatedMSA_FASTA} {output.ConcatenatedMSA_PHYLIP} {output.CodemlOut} \"{params.SpTreeConstrainTree}\" > {log.out} 2> {log.err}
		"""











