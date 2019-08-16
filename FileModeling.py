from collections import defaultdict
from abc import ABCMeta, abstractmethod
from collections import Counter
import itertools
import re
import optparse
import sys

#blaorden = open("/home/janosch/PycharmProjects/lncRNA/Syntheny_search/syntenizer/Bla_orden_genes_extendedI")
#genomefile = open("/home/janosch/PycharmProjects/lncRNA/Syntheny_search/hgenome_adj.txt")

# at first is a function file to produce databases for a __main__Script later on

class prepare_file(object):
    ''' Metaclass - a Parent class just for inheritance, to derive properties to child classes
    can not be used to create own object instances! For this use the children
    '''
    __metaclass__ = ABCMeta

    def __init__(self, file):
        self.file= file


    def splitF (self, delim="\t", to_strip ="\n"):
        '''
        just to make it more comfortable to create a spriped, splitted, file
        :param delim: default tab
        :param to_strip: default \n
        :return: new readlines() file
        '''
        content = self.file.readlines()
        for line in content:
            split_line = line.strip(to_strip).split(delim)
            content[content.index(line)] = split_line
        return content



    def control_input (self, name, position, chrom, strand):
        '''
        control function to check if the file are of the correct style (comma or tab separated)
        :param name: column in the file for gene name
        :param position: column in the file for gene position (integer)
        :param chrom: chromosome column in the file
        :param strand: pos of strand in the file
        :return: OK or Not Okay (error messages)
        '''
        name_list = []
        order_dict ={}
        counter = 1
        strand_sym = ["+","-","1","+1","-1"]
        # remove header
        #self.file.readlines(0)
        for line in self.file:
            line = line.strip("\n")
            line = re.split("[,\t]",line)
            try:
                gene_post = int(line[position])
            except ValueError:
                print(self.file.name)
                print("\n!ERROR IN GENE POSTITION COLUMN: ’{}’\n".format(line[position]))
                print(line)
                raise ValueError
            if not line[strand] in strand_sym:
                print(self.file.name)
                print("\n!ERROR IN GENE STRAND COlUMN: ’{}’\n".format(line[strand]))
                print(line)
                raise ValueError
            order_dict[counter] = (line[chrom],gene_post)
            name_list.append(line[name].upper())
            counter +=1
        if not len([x for x in name_list if "HOX" in x]) >= 0:
            print("\n!ERROR IN GENE NAME COLUMN \n{}\n- "
                  "lack of hox in column ’{}’ may indicate wrong column"
                  " or lack of data\n".format(self.file.name,name))
            raise EnvironmentError

        for gene in range(1, len(order_dict)):
            if order_dict[gene][0] == order_dict[gene+1][0]:
                if order_dict[gene][1] > order_dict[gene+1][1]:
                    print(self.file.name, order_dict[gene],order_dict[gene+1], " THESE TWO ARE NOT PROPERLY SORTED")

        return print("\nFILE SEEMS OK: {}\n".format(str(self.file.name)))








    def wordCount(self,file_list, show = "n"):
        ''' counts the abundance of single words in an file.readline()
        :param file_list = file.readline()
        :param show: print the results to the screen? "n" or "y".. default = n
        :return: a dictionary of all single words and their counts
        '''
        reading = list(itertools.chain.from_iterable(file_list))
        wordcount = Counter(reading)
        if show == "y":
            for item in wordcount.items(): print("{}\t{}".format(*item))
        return wordcount

    def compareLines(self,list, pattern):
        ''' Compares a line with the next following line in a list_file.
        Entries/lines can be strings [asdf, ...] or lists itself [[],[]]
        :param pattern: either string to search for in line or List_pos. if line is a list
        :return: List with positions of doubling
        '''
        content = list
        line1 = content[0]
        line2 = content[1]
        double_list = []
        x = 0
        while line1 and line2:
            if type(line1) == str:
                if pattern in line1 and pattern in line2:
                    double_list.append(content.index(line2.split("\t")))

            elif type(pattern) == int:
                if line1[pattern] == line2[pattern]:
                    double_list.append(content.index(line2))

            else:
                if pattern in "\t".join(line1) and pattern in "\t".join(line2):
                    double_list.append(content.index(line2))
            try:
                x += 1
                line1 = content[x]
                line2 = content[x+1]
            except IndexError:
                break

        return double_list

    @abstractmethod
    def file_type(self):
        """"Return a string representing the type of the file."""
        pass


###### Children clase to create particular object instances
###### They derive their properties form the prepare_file class

class casual (prepare_file):
    '''without any further specification / just normal text files'''

    def has_type(selfs):
        '''Return a string representing the type of vehicle this is'''

        return "text file"


class hash_file (prepare_file):
    """hash a file."""

    def hashBlast(self):
        """Return a dictionary with blast query as keys and values = [[gene_name_hit, evalue],[x,y],...]
        filters unique combination of query and hit (query isoforms)"""
        unique_dict = {}
        dict_blast = defaultdict(list)

        for line in self.file:
            if not line.startswith("#"):
                line = line.strip("\n").split("\t")
                # should always be in the same position
                query, gene_name, evalue = line[0].split("_")[0], line[1].split("|")[-1], float(line[-2])
                amphi_ensembl = gene_name + query
                if amphi_ensembl not in unique_dict.keys():
                    # unique query/name-combination (isoforms)
                    unique_dict[amphi_ensembl] = "once"
                    ranking = [gene_name, evalue]
                    dict_blast[query].append(ranking)
        return dict_blast

    def hashGenome(self):
        '''
        A ORDERED gene file with information of gene_name, chrom, strand
        :return: hash in style of key = geneName; value = geneName|position in genome|chrom|strand
        '''
        g_hash = {}
        content = self.file.readlines()
        for anno in content:
            if content.index(anno) != 0:
                pos_counter = content.index(anno)                  # genome pos.
                anno = anno.strip("\n").split("\t")
                gene, strand= anno[1], anno[3]

                try:
                    chrom = int(anno[2].upper().replace("Y","100").
                                replace("X","0").strip("CHR"))
                except ValueError:
                    print("#ERROR: invalid literal for int()-->",anno[2])
                    print("#check chromosome in line {1}\t{0}"
                          .format(anno, pos_counter))
                    break

                g_hash[gene] = "{0}|{1}|" \
                               "{2}|{3}".format(gene,str(pos_counter),chrom,str(strand))
        return g_hash



    def hash_type(self):
        """"Return a string representing the type of file."""

        return 'file to hash'




##### import file function

def flag_import(script):
    sys.argv[1:]
    parser = optparse.OptionParser()
    parser.add_option('-o', '--output',
                    dest="output_filename",
                    default="default.output",
                    )
    parser.add_option('-s', '--sortedFile',
                    dest="sorted_file",
                    default="orderedFile path error"
                    )
    parser.add_option('-d', '--blastdb',
                    dest="blast_hash",
                    default="blastp path error"
                    )
    parser.add_option('-g', '--genomedb',
                      dest="genome_hash",
                      default="genomep path error"
                      )
    parser.add_option('-t', '--toHashGenome',
                      dest="Genome_to_hash",
                      default="GENOME TO HASH path error"
                      )
    parser.add_option('-b', '--toHashBlast',
                      dest="Blast_to_hash",
                      default="BLAST to HASH path error"
                      )
    options, remainder = parser.parse_args()
    if script == "find_synteny":
        return options.sorted_file, options.output_filename
    elif script == "makeHash":
        return options.Genome_to_hash, options.Blast_to_hash
    else:
        return remainder