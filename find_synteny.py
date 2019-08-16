from operator import itemgetter
import pickle
import FileModeling
import re
import sys
import cluster_syntheny
import os


def mapping_files(ordering, ghash, dict_blast):
    '''
    mapping for clustering
    :param ordering: Bla_orden_file.spliF(",")
    :param ghash: ghash database
    :param dict_blast: Blast database
    :return: BLidentifier mapped with best blast hit + data about hit position/strand/chromosome
    '''
    for row in ordering:
        # ['Sc0000240', '75054', '90098', 'BL38170', 'Hypnc5241', '+']
        BL_id = row[3]
        if BL_id in dict_blast.keys() and not re.match("Hypnc", row[4]):
            sorted_hits = sorted(dict_blast[BL_id], key=itemgetter(1), reverse=False)[:5]
            gene_list = []
            for x in sorted_hits:
                if x[0] in ghash.keys():
                    # x[0] = gene_name;
                    # ghash[x[0]] = gene|chr|..
                    gene = ghash[x[0]]
                # incomplete entries
                elif x[0] in ghash.keys() and len(ghash[x[0]]) != 4:
                    continue
                # missing
                else:
                    continue
                    # alternative: x[0]
                gene_list.append(gene)
            ordering[ordering.index(row)] = [row[0],BL_id, "\t".join(gene_list)]

        elif re.match("Hypnc", row[4]):
            ordering[ordering.index(row)] = [row[0],BL_id, row[4]]
            # to adjust hypnc lines to gene lines
        else:
            # if nothing in BLAST
            ordering[ordering.index(row)] = [row[0],BL_id, "none"]
    return ordering



def first_file(sorting_file, ghash, dict_blast, output):
    '''
    :param ordering: Bla-orden.splitF(,)
    :param ghash: hashed genome file (second spec. )
    :param dict_blast: hashed blast file
    :param output: file for clustering
    :return: written output file
    '''
    # loading, preparing files
    sorting_file = FileModeling.casual(sorting_file)
    ordering= sorting_file.splitF(delim=",")
    # header
    print("amphiID\thuman_name1\thuman_name2"
          "\thuman_name3\thuman_name4\thuman_name5", file=output)

    # merging all to one -> scaf89 BL3434 GENE|chr|..\tGENE2
    mapping_files(ordering,ghash,dict_blast)
    # 1 find hypnc cluster positions
    cluster_pos = sorting_file.compareLines(ordering,"Hypnc")
    # 1.1 merge hypnc clusters to single entry before poping excess
    for x in sorted(cluster_pos, reverse=True):
        ordering[x-1][-1] = "{}&{}".format(ordering[x][-1],ordering[x-1][-1])
    # 2 find identical entries in a row
    double_pos = sorting_file.compareLines(ordering,2)
    # 3 find blank and "none" lines
    blank_pos = [ordering.index(line) for line in ordering if line[2] == "" or line[2] == "none"]
    # 4 remove cluster, double entries, blanks and nones
    pos_2_pop = sorted(set(cluster_pos)|set(double_pos)|set(blank_pos), reverse=True)
    [ordering.pop(x) if not x == len(ordering) else ordering.pop(x-1) for x in pos_2_pop]
    # 5 counts scaffolds
    word_abund = sorting_file.wordCount(ordering, show="n")

    # print to new file if: !
    for i in ordering:
        index = ordering.index(i)
        # no crappy scaffolds
        if word_abund[i[0]] > 7:
            # Hypn only if 3 up, down same scaff
            if re.match("Hypnc",i[2]):
                if index < len(ordering)-3 \
                        and i[0] == ordering[index+3][0] == ordering[index-3][0]:
                    output.write("\t".join(i[1:])+"\n")
            # Non-hypn in any case
            else:
                output.write("\t".join(i[1:])+"\n")
        # omit crappy scaffs
        else:
            continue

    return output


### first functions to generate the first output-file
### following functions used to go on with this output
### actual clustering takes place


def split_line(line, line_index, name_length):
    '''
    :param line: BL GENE|.. GENE|
    :param line_index: line numer in file
    :param name_length: how many genes already collected
    :return: names, chr., etc vectors
    '''
    name,chrom, genomepos, strand=[],[],[],[]
    #amphi ids, same length
    amphi= [line[0]] * len(line[1:])
    for z in line[1:]:
        # z for single genes
        z =z.split("|")
        # name vector stores: name of human, position of BL in dict, position of gene in nam_vector
        # TPH#6.10 = gene: TPH, 6 = position in dict_file (the key) . 10= position 10 in name_vector
        name.append(z[0]+"#"+str(line_index)+","+str(len(name)+name_length))
        chrom.append(int(z[2]))
        genomepos.append(int(z[1]))
        strand.append(z[-1])
    return amphi,name,genomepos,chrom,strand


def threestream(file_list, hypnc, dir):
    '''
    checks 3 upstream or downstream if there are other hypnc
    in the range of possible synteny close to a hypnc
    :param file_list: .readline()
    :param hypnc: line number of hypnc
    :param dir: upstream = +1; downstream = -1
    :return: line number where to start downstream / end upstream
    '''
    if type(dir) != int:
        raise ("dir is an integer")

    if not file_list[hypnc + (1 * dir)][1].startswith("Hypnc"):

        if not file_list[hypnc + (2 * dir)][1].startswith("Hypnc"):

            if not file_list[hypnc + (3 * dir)][1].startswith("Hypnc"):
                line_number = hypnc + (3 * dir)
            else:
                line_number = hypnc + (2 * dir)
        else:
            line_number = hypnc + (1 * dir)

    else:
        line_number = hypnc

    return line_number

def coding_neighbors(file_list, hypnc, dir):
    '''
    searches for 3 protein coding genes upstream/downstream of a hypnc
    other hypnc or empty entries (no human ortholog) will be omitted
    :param file_list: file_list: .readline()
    :param hypnc: line number of hypnc
    :param dir: upstream = +1; downstream = -1
    :return: a list of line numbers of protein coding genes
    '''
    if type(dir) != int:
        raise ("dir is an integer")
    line_index = 1
    gene_index_list = []
    while len(gene_index_list) < 3:
        if hypnc + (line_index * dir) == 0:
            break
        if hypnc + (line_index * dir) == len(file_list):
            break
        if not len(file_list[hypnc + (line_index * dir)]) == 1:
            gene = file_list[hypnc + (line_index * dir)][1]
            if not gene.startswith("Hypnc"):
                gene_index_list.append(hypnc + (line_index * dir))
        line_index += 1
    return gene_index_list



def final_file(first_file_list, final_output):
    '''
    applies the clustering and prints the final synteny results
    :param first_file_list: an ordered orthologs files
    :param final_output: synteny results
    :return:
    '''
    gen_translation = FileModeling.casual(first_file_list).splitF()
    for line in gen_translation:
        # line = BL GENE|asdf Gene|sdfd
        if len(line) == 1:
            continue
        if line[1].startswith("Hypnc"):
            hypnc_name = line[1]
            hypnc = gen_translation.index(line)
            # generating vectors
            name,chrom, genomepos, strand, amphi=[],[],[],[],[]
            #print(line)
            # look 3 up and down
            downstream_indices = coding_neighbors(gen_translation, hypnc, dir = -1)
            upstream_indices = coding_neighbors(gen_translation, hypnc, dir = 1)
            list_of_neighbors = sorted(downstream_indices + upstream_indices)

            # rna position with respect to all downstream genes (eg. 3 lines * 5 genes per line)
            rna_pos = sum([len(gen_translation[y][1:]) for y in downstream_indices]) - 0.5
            for line_index in [x for x in list_of_neighbors]:
                amphi.extend(split_line(gen_translation[line_index], line_index, len(name))[0])
                name.extend(split_line(gen_translation[line_index], line_index, len(name))[1])
                genomepos.extend(split_line(gen_translation[line_index], line_index, len(name))[2])
                chrom.extend(split_line(gen_translation[line_index], line_index, len(name))[3])
                strand.extend(split_line(gen_translation[line_index], line_index, len(name))[4])
            goodness,seed,path,path_len, all_path=cluster_syntheny.search_syntheny(name,chrom,genomepos,strand,amphi,rna_pos)
            print("{}\t({})\t>{}\t<{}\t#{}\n\n{}\n____".format(hypnc_name,str(path_len),str(goodness),str(seed), path, all_path),
                  file=final_output)

    return final_output



# controling the commands

'''
def main():
    # parse command line input
    blastdb,genomedb,sorted_file,outfile=FileModeling.flag_import(script="find_synteny")
    # loading dictionary of a Blastfile: key=amphiID value=[[gene,evalue],[gene,evalue]]
    dict_blast = pickle.load(open(blastdb,"rb"))
    # loading dictionary of SORTED genome file
    ghash = pickle.load(open(genomedb, "rb"))
    # used to sort the amphi id according to their position on the chr. (Bla_orden_extended file)
    order_file = open(sorted_file)
    FileModeling.casual(order_file).control_input(4,1,0,5)
    order_file.seek(0)
    # first output of the mapping
    output = open(outfile + ".txt", "w")

    # final outputfile (Synteny)
    finPath, outName = os.path.split(outfile)[0:2]
    finName = "Synteny_"+outName.split(".")[0] + ".txt"
    final_results = open(os.path.join(finPath,finName), "w")

    # template for
    print("\n...merging data\n")
    first_file(order_file,ghash,dict_blast,output)
    output.close()
    order_file.close()
    print("# merged file:\t" + os.getcwd() +"/"+ str(finPath+outName))
    # load the first outputfile again
    first_output=open(outfile + ".txt")
    print("\n...searching for synteny\n")

    final_file(first_output,final_results)
    first_output.close()
    final_results.close()
    print("# synteny results:\t" + os.getcwd() +"/"+ str(os.path.join(finPath,finName)))
    print ("\n\n###### Ya está! ######\n")
'''


def main():
    '''
    execute only as shell script. In contrast to main(), main_second()
    uses an already generated, ordered file (sorted file). This file contains
    the odered amphi genes and the human orthologs. The original main function generates this file
    as intermediate output. It requires output from the the makeHash function which itself requires
    the species genomes and a blast output file.
    :return: final synteny file
    '''
    # parse command line input
    sorted_file, outfile = FileModeling.flag_import(script="find_synteny")

    # final outputfile (Synteny)
    outPath, outName = os.path.split(outfile)[0:2]
    resultsName = "Synteny_"+outName.split(".")[0] + ".txt"
    synteny_results = open(os.path.join(outPath,resultsName), "w")

    # load the first outputfile again
    ordered_orthologs=open(sorted_file)
    print("\n...searching for synteny\n")

    final_file(ordered_orthologs,synteny_results)
    ordered_orthologs.close()
    synteny_results.close()
    print("# synteny results:\t" + os.getcwd() +"/"+ str(os.path.join(outPath,resultsName)))
    print ("\n\n###### Ya está! ######\n")


if __name__ == "__main__":
    # execute only form shell
    main()





#### TO DO
# Add the possiblity to let the user decide if he wants
# use an pre-ordered orthologs file between the two species
# or if he wants this script to generate this pre-ordered file
# basted on the genomes of the 2 species. This is already possible
# ... see the out-commented main function. It just requires to run the make hash script and the files
# must have a specified format.