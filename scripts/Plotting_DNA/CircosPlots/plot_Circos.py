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
fileID = sys.argv[7]

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


if "ColapsedPropAlt" in guide:
	# Collapsed line plots 
	print("--> Collapsed line plots")
	print('--> Collapsed line plots', file=sys.stderr)
	clp_filebasename = "line_propalt"
	clp_rwidth = 0.05
	clp_rstart = lp_rstarts[-1]
	clp_rlim = [clp_rstart, clp_rstart+clp_rwidth]
	clp_ylim = [0, .2]

	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"callablelength\"){print i}}}'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colcallable = result.stdout.split("\n")[0]
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{str=\"\"; for(i=1;i<=NF;i++){if($i ~ /altAl/){str=str\",\"i}}print str}' | sed 's/^,//g'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colAlt = result.stdout.split("\n")[0].split(",")
	bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{str=\"\"; for(i=1;i<=NF;i++){if($i ~ /altAl/){str=str\",\"$i}}print str}' | sed 's/^,//g'"
	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	colnamesAlt = [x[5:-1] for x in result.stdout.split("\n")[0].split(",")]
	SampleColors = ["orange" if x[0]=="R" else "purple" for x in colnamesAlt]
	
	for j in range(len(colAlt)):
		os.system("rm " + outdir + "/" + clp_filebasename + "_" + colnamesAlt[j] + "_for_" + fileID + ".tmp 2> ~/null")
	
	for i in WindMatFiles:
		for j in range(len(colAlt)):
			os.system("zcat " + i + " | cut -f1,2,3," + colcallable + "," + colAlt[j] + " | tail -n +2 | awk '{if($4>0){val=$5/$4;print $1\"\t\"$2\"\t\"$3\"\t\"val}}' >> " + outdir + "/" + clp_filebasename + "_" + colnamesAlt[j] + "_for_" + fileID + ".tmp")


# Scatter plots
#print("--> Scatter plots")
#sp_files = ["scatter_PC1"]
#sp_rstarts = [0.48]
#sp_rwidth = 0.1
#sp_rlim = [sp_rstarts, [x+sp_rwidth/2 for x in sp_rstarts], [x+sp_rwidth for x in sp_rstarts]]
#sp_ylim = [	[-7], 
#			[7.75]]
#sp_midrange = [	[-6], 
#				[5]]
#sp_axismargin = [0.01]
#bash_command = "zcat " + SiteMatFiles[0] + " | head -1 | awk '{for(i=0;i<=NF;i++){if($i==\"PC1all\"){print i}}}'"
#result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#colPC1 = result.stdout.split("\n")[0]
#for i in SiteMatFiles:
#	print(i)
#	bash_command = "zcat " + SiteMatFiles[0] + " | head -1 | awk '{for(i=0;i<=NF;i++){if($i==\"PC1all\"){print i}}}'"
#	result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#	colPC1 = result.stdout.split("\n")[0]
#	os.system("zcat " + i + " | cut -f1,2,3," + colPC1 + " | tail -n +2 | grep -v 'NA$' | awk '{if($4 <= " + str(sp_midrange[0][0]) + "){print $1\"\t\"$2\"\t\"$3\"\t\"$4}}' >> " + outdir + "/scatter_PC1.tmp")
#	os.system("zcat " + i + " | cut -f1,2,3," + colPC1 + " | tail -n +2 | grep -v 'NA$' | awk '{if($4 >= " + str(sp_midrange[1][0]) + "){print $1\"\t\"$2\"\t\"$3\"\t\"$4}}' >> " + outdir + "/scatter_PC1.tmp")
#print("       Bottom")
#os.system("cat " + outdir + "/scatter_PC1.tmp | sort -k4,4g | tail -50000 >> " + outdir + "/scatter_PC1bottom.tmp")
#bash_command = "cat " + outdir + "/scatter_PC1bottom.tmp | awk '{max=-10000; min=10000; if($4>max){max=$4}; if($4<min){min=$4}}END{print min,max}'"
#result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#minmax = result.stdout.split("\n")[0].split(" ")
#print(minmax)
#sp_ylim[0][0] = float(minmax[0]) - sp_axismargin[0]
#sp_midrange[0][0] = float(minmax[1]) + sp_axismargin[0]
#print("       Top")
#os.system("cat " + outdir + "/scatter_PC1.tmp | sort -k4,4g | awk '{if(NR <= 50000){print $0}}' >> " + outdir + "/scatter_PC1top.tmp")
#bash_command = "cat " + outdir + "/scatter_PC1top.tmp | awk '{max=-10000; min=10000; if($4>max){max=$4}; if($4<min){min=$4}}END{print min,max}'"
#result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#minmax = result.stdout.split("\n")[0].split(" ")
#print(minmax)
#sp_midrange[1][0] = float(minmax[0]) - sp_axismargin[0]
#sp_ylim[1][0] = float(minmax[1]) + sp_axismargin[0]





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
##################################
## Collapsed line plots
if "ColapsedPropAlt" in guide:
	for j in range(len(colAlt)):
		conf_text = conf_text + " \n \
		<plot>\n \
			type    = line\n \
			max_gap = 1u\n \
			file    = " + outdir + "/" + clp_filebasename + "_" + colnamesAlt[j] + "_for_" + fileID + ".tmp\n \
			color   = " + SampleColors[j] + "\n \
			min     = " + str(clp_ylim[0]) + "\n \
			max     = " + str(clp_ylim[1]) + "\n \
			r0      = " + str(clp_rlim[0]) + "r\n \
			r1      = " + str(clp_rlim[1]) + "r\n \
			thickness = 1\n \
		"
		if j == 0:
			conf_text = conf_text + " \n \
				<axes>\n \
					<axis>\n \
						color     = vvlgrey\n \
						thickness = 1\n \
						spacing   = 0.1r\n \
					</axis>\n \
					</axes>\n \
			"
		conf_text = conf_text + " \n \
		</plot>\n \
		"


###################################
### Scatter plots
#for sp in range(len(sp_files)):
#	conf_text = conf_text + " \n \
#	<plot>\n \
#		type    = scatter\n \
#		stroke_thickness = 1\n \
#		file    = " + outdir + "/" + sp_files[sp] + "bottom.tmp\n \
#		color   = black\n \
#		fill_color       = dred\n \
#		stroke_color     = dred\n \
#		glyph            = circle\n \
#		glyph_size       = 10\n \
#		min     = " + str(sp_ylim[0][sp]) + "\n \
#		max     = " + str(sp_midrange[0][sp]) + "\n \
#		r0      = " + str(sp_rlim[0][sp]) + "r\n \
#		r1      = " + str(sp_rlim[1][sp]) + "r\n \
#		<backgrounds>\n \
#			<background>\n \
#				color     = vvlred\n \
#			</background>\n \
#		</backgrounds>\n \
#		<axes>\n \
#			<axis>\n \
#				color     = black\n \
#				thickness = 1\n \
#				spacing   = 0.1r\n \
#			</axis>\n \
#			</axes>\n \
#	</plot>\n \
#	<plot>\n \
#		type    = scatter\n \
#		stroke_thickness = 1\n \
#		file    = " + outdir + "/" + sp_files[sp] + "top.tmp\n \
#		color   = black\n \
#		fill_color       = dgreen\n \
#		stroke_color     = dgreen\n \
#		glyph            = circle\n \
#		glyph_size       = 10\n \
#		min     = " + str(sp_midrange[1][sp]) + "\n \
#		max     = " + str(sp_ylim[1][sp]) + "\n \
#		r0      = " + str(sp_rlim[1][sp]) + "r\n \
#		r1      = " + str(sp_rlim[2][sp]) + "r\n \
#		<backgrounds>\n \
#			<background>\n \
#				color     = vvlgreen\n \
#			</background>\n \
#		</backgrounds>\n \
#		<axes>\n \
#			<axis>\n \
#				color     = black\n \
#				thickness = 1\n \
#				spacing   = 0.1r\n \
#			</axis>\n \
#			</axes>\n \
#	</plot>\n \
#	"

conf_text = conf_text + " \n \
</plots>\n \
<<include scripts/Plotting_DNA/Circos_housekeeping_modif.conf>>\n \
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





















