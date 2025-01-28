#!/usr/bin/python3

################
## Debug & develop
#./scripts/Plotting_DNA/plot_Circos.py data/genomes/BraLan3/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt "results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr1.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr2.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr3.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr4.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr5.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr6.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr7.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr8.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr9.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr10.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr11.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr12.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr13.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr14.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr15.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr16.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr17.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr18.tab.gz results/VariantAnalysis_DNA/PerWindow_StatsMatrix_PerChr/PerWindowStats_100000_10000.chr19.tab.gz" "results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr1.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr2.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr3.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr4.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr5.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr6.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr7.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr8.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr9.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr10.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr11.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr12.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr13.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr14.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr15.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr16.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr17.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr18.tab.gz results/VariantAnalysis_DNA/PerSite_StatsMatrix_PerChr/PerSiteStats.chr19.tab.gz" results/Plotting_DNA/plot_Circos/karyotype_Callable_Varsites_Pi_PC1_GC_ColapsedPropAlt.txt results/Plotting_DNA/plot_Circos/circos_Callable_Varsites_Pi_PC1_GC_ColapsedPropAlt.conf results/Plotting_DNA/plot_Circos/plot_Circos_Callable_Varsites_Pi_PC1_GC_ColapsedPropAlt.png Callable_Varsites_Pi_PC1_GC_ColapsedPropAlt

import sys
import os
import subprocess

chrfile = sys.argv[1]
WindMatFiles = sys.argv[2].split(" ")
SiteMatFiles = sys.argv[3].split(" ")
strSiteMatFiles = sys.argv[3]
karyotype = sys.argv[4]
configfile = sys.argv[5]
out = sys.argv[6]
housekeepingconf = sys.argv[7]
fileID = sys.argv[8]

###############
## Initiate + general parameters
os.system("mkdir -p $(dirname "+ karyotype + ")")
outdir = "/".join(karyotype.split("/")[0:-1])
tmpfile = outdir + "/tmpfile" 
outbasename = out.split("/")[-1]
chrcolor = "black"
guide = fileID.split("_")



###############
## Karyotype file
os.system("cat " + chrfile + " | grep '^>chr' | sed 's/>//g' > " + tmpfile)
chrfile_r = open(tmpfile, 'r')
chrinfo = []
for l in chrfile_r:
	line=l.split("\n")[0].split("\t")
	chrinfo.append(line)
chrfile_r.close()

karyotype_w = open(karyotype, 'w') #chr - ID LABEL START END COLOR
for c in chrinfo:
	line = ["chr", "-", c[0], c[2], "0", c[1], chrcolor]
	karyotype_w.write("\t".join(line) + "\n")
karyotype_w.close()

###############
## Plot files

# Line plots 
print("--> Line plots")
print('--> Line plots', file=sys.stderr)
lp_files = []
lp_rwidth = 0.05
lp_rstarts = [0.9]
for i in range(1, len(guide)):
	lp_rstarts.append(lp_rstarts[-1]-lp_rwidth)

lp_rlim = [lp_rstarts, [x+lp_rwidth for x in lp_rstarts]]
lp_ylim = [	[],[] ]

## define tmp files and axis limits
for i in guide:
	if "Callable" == i:
		lp_files.append("line_propcallable_for_" + fileID + ".tmp")
		lp_ylim[0].append(0)
		lp_ylim[1].append(1)
	if "Varsites" == i:
		lp_files.append("line_varsites_for_" + fileID + ".tmp")
		lp_ylim[0].append(0)
		lp_ylim[1].append(1)
	if "Pi" == i:
		lp_files.append("line_pi_for_" + fileID + ".tmp")
		lp_ylim[0].append(0)
		lp_ylim[1].append(.2)
	if "PC1" == i:
		lp_files.append("line_PC1_for_" + fileID + ".tmp")
		lp_ylim[0].append(-2)
		lp_ylim[1].append(2)
	if "pFst" == i:
		lp_files.append("line_pfst_for_" + fileID + ".tmp")
		lp_ylim[0].append(0)
		lp_ylim[1].append(0.1)
	if "sFst" == i:
		lp_files.append("line_sfst_for_" + fileID + ".tmp")
		lp_ylim[0].append(0)
		lp_ylim[1].append(0.02)
	if "GC" == i:
		lp_files.append("line_GC_for_" + fileID + ".tmp")
		lp_ylim[0].append(0.3)
		lp_ylim[1].append(0.6)


if "Callable" in guide:
	print("    --> Callable")
	print('    --> Callable', file=sys.stderr)
	os.system("rm " + outdir + "/line_propcallable_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"callablelength\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colcallable = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3," + colcallable + " | tail -n +2 | awk '{w=$3-$2; val=$4/w; print $1\"\t\"$2\"\t\"$3\"\t\"val}' >> " + outdir + "/line_propcallable_for_" + fileID + ".tmp")


if "Varsites" in guide:
	print("    --> Varsites")
	print('    --> Varsites', file=sys.stderr)
	os.system("rm " + outdir + "/line_varsites_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"varsites\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colvarsites = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3," + colvarsites + "," + colcallable + " | tail -n +2 | awk '{if($5>0){val=$4/$5}else{val=0};print $1\"\t\"$2\"\t\"$3\"\t\"val}' >> " + outdir + "/line_varsites_for_" + fileID + ".tmp")

if "Pi" in guide:
	print("    --> Pi")
	print('    --> Pi', file=sys.stderr)
	os.system("rm " + outdir + "/line_pi_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"PI\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colpi = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3," + colpi + " | tail -n +2 >> " + outdir + "/line_pi_for_" + fileID + ".tmp")

if "PC1" in guide:
	print("    --> PC1")
	print('    --> PC1', file=sys.stderr)
	os.system("rm " + outdir + "/line_PC1_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"PC1all\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colPC1 = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3," + colPC1 + " | tail -n +2 >> " + outdir + "/line_PC1_for_" + fileID + ".tmp")

if "pfst" in guide:
	print("    --> pfst")
	print('    --> pfst', file=sys.stderr)
	os.system("rm " + outdir + "/line_pfst_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"mean_pop_FST\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colpfst = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3," + colpfst + " | tail -n +2 >> " + outdir + "/line_pfst_for_" + fileID + ".tmp")

if "sfst" in guide:
	print("    --> sfst")
	print('    --> sfst', file=sys.stderr)
	os.system("rm " + outdir + "/line_sfst_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"mean_sex_FST\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colsfst = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3," + colsfst + " | tail -n +2 >> " + outdir + "/line_sfst_for_" + fileID + ".tmp")

if "GC" in guide:
	print("    --> GC")
	print('    --> GC', file=sys.stderr)
	os.system("rm " + outdir + "/line_GC_for_" + fileID + ".tmp 2> ~/null")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{str=\"\";for(i=1;i<=NF;i++){if($i==\"ATcontent\" || $i==\"GCcontent\"){str=str\",\"i}}print str}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colsATGC = result.stdout.split("\n")[0]
	for i in WindMatFiles:
		os.system("zcat " + i + " | cut -f1,2,3" + colsATGC + " | awk '{t=$4+$5; if(t>0){v=$5/t}else{v=0}; print $1\"\t\"$2\"\t\"$3\"\t\"v}' | tail -n +2 >> " + outdir + "/line_GC_for_" + fileID + ".tmp")


###############
## Config file
conf_w = open(configfile, 'w')
conf_text="karyotype = " + karyotype + "\n \
<ideogram>\n \
	<spacing>\n \
		default = 0.005r\n \
	</spacing>\n \
	radius = 0.9r\n \
	thickness = 20p\n \
	fill = yes  \n \
	stroke_color = dgrey\n \
	stroke_thickness = 2p   \n \
	show_label       = yes\n \
	label_font       = default\n \
	label_radius     = 1.01r\n \
	label_with_tag   = yes\n \
	label_size       = 36\n \
	label_parallel   = yes\n \
	label_case       = lower\n \
	#label_format     = eval(sprintf(\"chr%s\",var(label)))\n \
</ideogram>\n \
<image>\n \
	dir = " + outdir + "\n \
	file  = " + outbasename + "\n \
	radius         = 1500p\n \
	background     = white\n \
	angle_offset   = -90\n \
</image>\n \
<plots>\n "

##################################
## Line plots
for lp in range(len(lp_files)):
	conf_text = conf_text + " \n \
	<plot>\n \
		type    = line\n \
		max_gap = 1u\n \
		file    = " + outdir + "/" + lp_files[lp] + "\n \
		color   = black\n \
		min     = " + str(lp_ylim[0][lp]) + "\n \
		max     = " + str(lp_ylim[1][lp]) + "\n \
		r0      = " + str(lp_rlim[0][lp]) + "r\n \
		r1      = " + str(lp_rlim[1][lp]) + "r\n \
		thickness = 1\n \
		fill_color = black\n \
		<axes>\n \
			<axis>\n \
				color     = lgreen\n \
				thickness = 1\n \
				spacing   = 0.1r\n \
			</axis>\n \
			</axes>\n \
	</plot>\n \
	"
conf_text = conf_text + " \n \
</plots>\n \
<<include " + housekeepingconf + ">>\n \
<colors>\n \
<<include etc/colors.conf>>\n \
</colors>\n \
<fonts>\n \
<<include etc/fonts.conf>>\n \
</fonts>\n \
"
conf_w.write(conf_text)
conf_w.close()






os.system("rm " + tmpfile + " 2> ~/null")
#os.system("rm " + outdir + "/line_*_for_" + fileID + ".tmp 2> ~/null")





















