# LincOFinder
LincOFinder is a program designed to identify shared microsyntenic clusters surrounding lincRNAs between two species. 

The main script usage is as follows "python3 find_synteny.py -s Input_order_orthologies.txt -o Output.out

-s sorted orthologs file 
-o name of the output file
**FileModeling.py and cluster_syntheny.py must be in the same folder to run**

-----INPUT-----
The sorted orthologs file must have the following format:
Ref_GeneID  Cluster1(Int_GeneID|virtual_coordinate|Chromosome|strand) Cluster2  Cluster3...ClusterN for coding genes
Ref_GeneID "Hypnc_"Name_of_the_lincRNA for lincRNAs

As in the following example:
BL24937	ENSG00000122565|141|7|+	ENSG00000108468|776|17|-	ENSG00000094916|482|12|-
BL06743	ENSG00000141293|778|17|-	ENSG00000005020|143|7|-
BL24569	ENSG00000128654|846|2|+
BL20528	Hypnc_BL20528
BL12289	ENSG00000105991|144|7|-	ENSG00000120094|779|17|-	ENSG00000128645|845|2|+
BL38782	Hypnc_BL38782
BL01409	ENSG00000105996|145|7|-	ENSG00000173917|780|17|-

Chromosomes must be numbers, we advise to change chromosome X to "0" and Y to "1000".

----OUTPUT----
The output must be read as follows:

lincRNA (number of genes in the best cluster) >"standard deviation of the distance between the further genes" <"standard deviation of the distance between the closer genes"  #BEST CLUSTER ->Int_GeneID_1,virtual_coordinate_1,Chromosome_1,Ref_GeneID_1 Int_GeneID_2,virtual_coordinate_2,Chromosome_2,Ref_GeneID_2 Int_GeneID_N,virtual_coordinate_N,Chromosome_N,Ref_GeneID_N

  Other possible clusters (one line for cluster)

An example of the output:
Hypnc_BL48403	(2)	>0.7	<0.7	#ENSG00000106443,Pos:83,Strand:+,Chr:7,BL:BL05756	ENSG00000189043,Pos:82,Strand:-,Chr:7,BL:BL11653

 	ENSG00000106443,Pos:83,Strand:+,Chr:7,BL:BL05756	ENSG00000189043,Pos:82,Strand:-,Chr:7,BL:BL11653	
	ENSG00000185633,Pos:579,Strand:-,Chr:12,BL:BL11653	ENSG00000139372,Pos:786,Strand:+,Chr:12,BL:BL24944	



A perl script is also included to help generating the INPUT file, usage available in preparation_README.md
