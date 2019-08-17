# LincOFinder  
LincOFinder is a program designed to identify shared microsyntenic clusters surrounding lincRNAs between two species. 

**Introduction**  
Lon-Non coding RNAs (lnc-RNAs) take part in a wide range of cellular and developmental processes. Unfortunately, week sequence conservation of lnc-RNAs between two species makes comparisons by sequence homology difficult. This program aims to find lnc-RNAs that are conserved by their surrounding protein coding genes between two species. This conservation of a genomic neighborhood is called “synteny”.  
  
**lincOfinder**  
In a first step, an ordered gene file (see input section) has to be generate, where every gene of Ref species is sorted by its genomic position. A second column lists the corresponding orthologs of Int species. Thus, each line contains a gene of Ref species and its corresponding, potential orthologs of Int species.  
  
In a first step of finding microsynteny, lincOfinder searches for linc-RNAs in the ordered gene file. Since microsynteny involves genes in close vicinity, lincOfinder then takes 3 genes upstream and 3 genes downstream of a linc-RNA in Ref species, together with all the orthologs of Int species. Finally, a UPGMA clustering on the genomic positions is performed to find combinations of orthologs that are close together. In the best case, one would find a cluster of Int species genes, where all the genes are their direct genomic neighbors. Given the fact, that the clustering was performed on Int species gene positions but the selection of genes was based on the 6 lincRNA surrounding genes of Ref species, such a cluster could suggest microsynteny. In case there are multiple clusters, which are equal in their distance, lincOfinder reports the best cluster as the one with the orthologs closest to the lincRNA. LincOfinder additionally reports the standard deviation of the best cluster. This is not related to the clustering process, but facilitates to assess cluster distances by eye. Furthermore, all other possible clusters are part of the output.  
  
The main script usage is as follows **python3 find_synteny.py -s Input_order_orthologies.txt -o Output.out**  

-s sorted orthologs file  
-o name of the output file  
**FileModeling.py and cluster_syntheny.py must be in the same folder to run**  

**-----INPUT-----**  
The sorted orthologs file must have the following format:  
Ref_GeneID  Ortholog1(Int_GeneID|virtual_coordinate|Chromosome|strand) Ortholog2  Ortholog3...OrthologN for coding genes  
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
  
**----OUTPUT----**  
The output must be read as follows:  
  
lincRNA (number of genes in the best cluster) >"standard deviation of distance of all genes" <"standard deviation of distance of closest genes"  #BEST CLUSTER ->Int_GeneID_1,virtual_coordinate_1,Chromosome_1,Ref_GeneID_1 Int_GeneID_2,virtual_coordinate_2,Chromosome_2,Ref_GeneID_2 Int_GeneID_N,virtual_coordinate_N,Chromosome_N,Ref_GeneID_N  
  
  Other possible clusters (one line for cluster)  
  
An example of the output:  
Hypnc_BL48403	(2)	>0.7	<0.7	#ENSG00000106443,Pos:83,Strand:+,Chr:7,BL:BL05756  	  ENSG00000189043,Pos:82,Strand:-,Chr:7,BL:BL11653  

 ENSG00000106443,Pos:83,Strand:+,Chr:7,BL:BL05756	ENSG00000189043,Pos:82,Strand:-,Chr:7,BL:BL11653	
 ENSG00000185633,Pos:579,Strand:-,Chr:12,BL:BL11653	ENSG00000139372,Pos:786,Strand:+,Chr:12,BL:BL24944	



A perl script is also included to help generating the INPUT file, usage available in preparation_README.md
  
**Required packages**  
  
- python 3.4 or higher  
- operator  
- abc  
- pickle  
- optparse  
- itertools  
- collections  
- Bio.Phylo.TreeConstruction  
- numpy  
- statistics  
  
