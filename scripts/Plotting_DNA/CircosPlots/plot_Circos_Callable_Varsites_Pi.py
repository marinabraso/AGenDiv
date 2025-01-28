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
noncall = sys.argv[4]
karyotype = sys.argv[5]
configfile = sys.argv[6]
out = sys.argv[7]
housekeepingconf = sys.argv[8]


###############
## Initiate + general parameters
os.system("mkdir -p $(dirname "+ karyotype + ")")
outdir = "/".join(karyotype.split("/")[0:-1])
tmpfile = outdir + "/tmpfile" 
outbasename = out.split("/")[-1]
chrcolor = "black"




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
lp_rwidth = 0.15
lp_rstarts = [.8]
lp_rstarts.append(lp_rstarts[-1]-lp_rwidth)
#lp_rstarts.append(lp_rstarts[-1]-lp_rwidth)

lp_rlim = [lp_rstarts, [x+lp_rwidth for x in lp_rstarts]]
lp_ylim = [	[],[] ]

## define tmp files and axis limits
lp_files.append("line_varsites_for_Callable_Varsites_Pi.tmp")
lp_ylim[0].append(0)
lp_ylim[1].append(1)
#lp_files.append("line_pi_for_Callable_Varsites_Pi.tmp")
#lp_ylim[0].append(0)
#lp_ylim[1].append(.2)
lp_files.append("line_exoncontent_for_Callable_Varsites_Pi.tmp")
lp_ylim[0].append(0)
lp_ylim[1].append(1)


bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"callablelength\"){print i}}}'"
result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
colcallable = result.stdout.split("\n")[0]


print("    --> Varsites")
print('    --> Varsites', file=sys.stderr)
os.system("rm " + outdir + "/line_varsites_for_Callable_Varsites_Pi.tmp 2> ~/null")
bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"varsites\"){print i}}}'"
result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
colvarsites = result.stdout.split("\n")[0]
for i in WindMatFiles:
	os.system("zcat " + i + " | cut -f1,2,3," + colvarsites + "," + colcallable + " | tail -n +2 | awk '{if($5>0){val=$4/$5}else{val=0};print $1\"\t\"$2\"\t\"$3\"\t\"val}' >> " + outdir + "/line_varsites_for_Callable_Varsites_Pi.tmp")

print("    --> Pi")
print('    --> Pi', file=sys.stderr)
os.system("rm " + outdir + "/line_pi_for_Callable_Varsites_Pi.tmp 2> ~/null")
bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"PI\"){print i}}}'"
result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
colpi = result.stdout.split("\n")[0]
for i in WindMatFiles:
	os.system("zcat " + i + " | cut -f1,2,3," + colpi + " | tail -n +2 >> " + outdir + "/line_pi_for_Callable_Varsites_Pi.tmp")

print("    --> Exon content")
print('    --> Exon content', file=sys.stderr)
os.system("rm " + outdir + "/line_exoncontent_for_Callable_Varsites_Pi.tmp 2> ~/null")
bash_command = "zcat " + WindMatFiles[0] + " | head -1 | awk '{for(i=1;i<=NF;i++){if($i==\"exoncontent\"){print i}}}'"
result = subprocess.run(bash_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
colexoncontent = result.stdout.split("\n")[0]
for i in WindMatFiles:
	os.system("zcat " + i + " | cut -f1,2,3," + colcallable + "," + colexoncontent + " | tail -n +2 | awk '{if($4>0){val=$5/$4}else{val=0};print $1\"\t\"$2\"\t\"$3\"\t\"val}' >> " + outdir + "/line_exoncontent_for_Callable_Varsites_Pi.tmp")

print("    --> Noncallable")
print('    --> Noncallable', file=sys.stderr)
os.system("rm " + outdir + "/line_Noncallable_for_Callable_Varsites_Pi.tmp 2> ~/null")
os.system("zcat " + noncall + " | cut -f1,2,3 | awk '{if($3-$2 > 100000){print $0}}' >> " + outdir + "/line_Noncallable_for_Callable_Varsites_Pi.tmp")


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
"
conf_text = conf_text + "\n \
<highlights>\n \
	<highlight>\n \
		z = 0\n \
		fill_color = vlgrey\n \
		file       = " + outdir + "/line_Noncallable_for_Callable_Varsites_Pi.tmp\n \
		r0         = " + str(lp_rlim[0][-1]) + "r\n \
		r1         = " + str(lp_rlim[1][0]) + "r\n \
	</highlight>\n \
</highlights>\n \
"
conf_text = conf_text + "\n \
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
"



conf_text = conf_text + " \n \
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
#os.system("rm " + outdir + "/line_*_for_Callable_Varsites_Pi.tmp 2> ~/null")





















