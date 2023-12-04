The preparation_input_tool.pl usage is the following:  
**perl preparation_input_tool.pl file_orthologies file_int_species file_ref_species output**  
file_orthologies: A file including all the orthology relathionships between the two species  
file_ref_species: A file including the coding and lincRNA genes for the Reference species ordered  
file_int_species: A file including the coding genes for the Interrogated species ordered  
output: Output file name (input for find_synteny.py)  
  
Some considerations for the format of the Input files below  
  
**----INPUT----  
-file_orthologies**  
the format must be:  
Orthologic_family_code(number)  Species GeneID  (Optional)Orthologic_subfamily_code(number)  
With the Int species always first for every orthologyc family  
For more info see example below:  
22	Hsa	ENSG00000127951	22  
22	Bla	BL74412	22  
29	Hsa	ENSG00000116194	29  
29	Hsa	ENSG00000136859	29  
29	Hsa	ENSG00000130812	3623  
29	Bla	BL05848	29  
29	Bla	BL12211	3623  
31	Hsa	ENSG00000101280	31  
31	Hsa	ENSG00000091879	31  
31	Hsa	ENSG00000154188	32  
31	Bla	BL13639	31  
31	Bla	BL24935	31  
31	Bla	BL00184	31  
34	Hsa	ENSG00000129083	34  
34	Bla	BL09202	34  
35	Hsa	ENSG00000166548	35  
35	Bla	BL09208	35  
  
**-file_ref_species**  
Format must be a comma separated file with:  
Chromosome,start,end,GeneID,TranscriptID,strand,class  
Ordered by start possition for every chromosome  
Must be collapsed to the canonical isoform for each gene (or the longest isoform if no other info is available)  
For more info see example below:  
Sc0000000,372950,374384,BL73495,BL73495_cuf0,+,coding  
Sc0000000,543362,545073,BL60128,BL60128_cuf0,+,linc  
Sc0000000,548933,607377,BL04777,BL04777_cuf6,+,coding  
Sc0000000,610969,614494,BL07583,BL07583_evm0,-,coding  
  
**-file_int_species**  
GTF file of the Interrogated species with only the coding genes  
Must be ordered by start possition for every chromosome  
The gene_id must be comparable with the one in file_orthologies  
For more info see example below:  
1	ensembl_havana	gene	65419	71585	.	+	.	gene_id "ENSG00000186092"; gene_version "6"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding";  
1	ensembl_havana	gene	450703	451697	.	-	.	gene_id "ENSG00000284733"; gene_version "1"; gene_name "OR4F29"; gene_source "ensembl_havana"; gene_biotype "protein_coding";  
1	ensembl_havana	gene	685679	686673	.	-	.	gene_id "ENSG00000284662"; gene_version "1"; gene_name "OR4F16"; gene_source "ensembl_havana"; gene_biotype "protein_coding";  
  
**----OUTPUT----**  
The output will be like the input described in README.md if everything went well  
