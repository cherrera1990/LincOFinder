# clustering script to identify syntheny
# the output is a vector of position values
# these values indicate the positions of the genes, in the input list, which make up the best hit for syntheny
import Bio.Phylo.TreeConstruction
#from Bio.Phylo.TreeConstruction import Matrix
import numpy
import re
import statistics


def generate_matrix(names):
    '''
    generate a list of lists: the rows grow with colum number -> col=1 row=1 etc
    afterwards shaped into a triangular matrix
    :param names: list of gene names (also used to define the size)
    :return: a triangular matrix with 10^10 as diagonal
    '''
    hole_list = []
    for x in range(1, len(names) + 1):
        line_list = []
        for y in range(0, x):
            line_list.append(10e10)
        hole_list.append(line_list)
    return Bio.Phylo.TreeConstruction._Matrix(names, hole_list)


def Initialize(names, chrom, position, RNApos):
    '''
    Function which conducts the first steps 1) first_mat 2) seed_cluster 2) first iteration 3) minimum control
    :param names: orig. names
    :param chrom: orig. chroms
    :param position: orig. vector of pos.
    :param RNApos: orig. RNApos
    :return: the minimum of the first iteration step, the best_cluster and
    the first iteration-matrix for further clustering. If no minmum, return:No_minimum, No_minimum, ini_mat
    ----> means "STOP" for all further clustering
    '''
    m = generate_matrix(names)
    ini_mat = distance_matrix(m, chrom, position, names)
    # import!!! ini is different from iter, has it uses a distinct function to find the minimum
    minimum = ini_minimum(ini_mat, RNApos)
    if minimum != "No_minimum":
        # best cluster a the first clustered name combi and acts a seed to find the best path
        best_cluster = names[minimum[0]] + "|" + names[minimum[1]]
    else:
        best_cluster = "No_minimum"
    try:
        # try iteration, else use ini_mat and stop downstream compluting with "No_minimum"
        iter_mat = collapse_matrix(ini_mat, names, position, chrom, type="ini", RNApos=RNApos)
        minimum = find_minimum(iter_mat)
    except ValueError:
        minimum = "No_minimum"
        iter_mat=ini_mat
    return minimum, best_cluster, iter_mat


def ini_minimum(dist_matrix, lncRNA):
    '''
    USE ONLY FOR INI-MATRIX because in case of multiple minimums, the lncRNA position is used evaluate the best min.
    :param dist_matrix: initial distance matrix
    :param value: theoretical position of the lncRNA in the gene names vector
    :return: the 2 surrounding values of the lncRNA position (e.g. lnc=14 -> return = 12,16 (the min closest to 14)
    '''
    num_mat = numpy.array(dist_matrix)
    # filter matrix which has no entries ... line1=hypno, line2=hyponc, line3=hyponc -> no values in the matrix
    if num_mat.shape[0] > 1:
        array = numpy.where(num_mat == num_mat.min())
    else:
        return "No_minimum"
    # if more than one minimum
    if len(array[0]) > 2:
        min_list = []
        # calculate which minimum coordinate combinations are the closest to the lncRNA pos.
        # array1[11,6,5] array2[6,7,7]  -> entry one of array1 belongs to entry 1 of array2...they form the min.
        for x in range(0, len(array[0])):
            min_list.append((array[0][x] + array[1][x]) / 2)  ## min.-combi.-list [5.5, 6.5,6]
        # take the one which is closest to lncRNA...if lnceRNA= 5 take [5.5]
        closest = min(min_list, key=lambda x: abs(x - lncRNA))
        # get the position in min list ...pos 0 is [5.5]
        close_idx = [i for i, j in enumerate(min_list) if j == closest][0]
        # convert min list position in minimum coordinates .. [11] [6]
        cordi1, cordi2 = array[0][close_idx], array[1][close_idx]
        # return them unless!!!! those coordinates yield 10e10
        if num_mat[cordi1, cordi2] != 10e10:
            return cordi1, cordi2
        else:
            return "No_minimum"
            # take the easy way if there is just one minimum
    else:
        cordi1, cordi2 = numpy.where(num_mat == num_mat.min())[0]
        if num_mat[cordi1, cordi2] != 10e10:
            return cordi1, cordi2
        else:
            return "No_minimum"


def find_minimum(dist_matrix):
    '''
    searches the minimum in a non-ini distance matrix ... it takes the minimum with the highest coordinates!
    This is, because clusters are appended to the running lists and have therefore the highest coordinates in the
    lists. Thus, by this, it is more likely to choose a minimum of an already generated
    cluster if more than one minimum exists.
    :param dist_matrix: distance matrix
    :return: the minimum coordinates
    '''
    num_mat = numpy.array(dist_matrix)
    array = numpy.where(num_mat == num_mat.min())
    if len(array[0]) > 2:
        cordi1 = array[0][-1]  # take minimum which appears last in the array (likely to be the min of an appended cluster)
        cordi2 = array[1][-1]
        if num_mat[cordi1, cordi2] != 10e10:
            return cordi1, cordi2
        else:
            return "No_minimum"
    else:
        cordi1, cordi2 = numpy.where(num_mat == num_mat.min())[0]
        if num_mat[cordi1, cordi2] != 10e10:
            return cordi1, cordi2
        else:
            return "No_minimum"


def collapse_list(list_input, coordHigh, coordLow):
    '''
    fucntion with merges 2 entries in a list. List sizes reduces by 1. THIS FUNCTION DOES THE ACTUAL CLUSTERING
    :param list_input: list to cluster
    :param coordHigh: highest coordinate (list position) of the matrix minima
    :param coordLow:  lowest coordinate (list position) of the matrix minima
    :return: a list where the list entries of the 2 coordinates are merged into one entry. This entry is added to the
    old rest of the list at the first position
    '''
    new_list = []
    if type(list_input[0]) == int or type(list_input[0]) == float:
        # if a number like chrom number or position, take the average of it
        # new position
        # chrom number should remain the same after deviation since they were the same before
        cluster = (list_input.pop(coordHigh) + list_input.pop(coordLow)) / 2
        new_list = list_input
        new_list.append(cluster)
    elif type(list_input[0]) == str:
        cluster = list_input.pop(coordHigh) + "|" + list_input.pop(coordLow)
        new_list = list_input
        new_list.append(cluster)
    return new_list


def distance_matrix(m, chrom, position, names):
    '''
    Function to compute distance values and fill with them a raw matrix
    :param m: raw matrix with default values
    :param chrom: list with chromosome numbers in order
    :param position: list with positions of the genes in same order
    :return: the distance matrix
    '''
    for i in range(len(chrom)):
        for j in range(i + 1, len(chrom)):
            if chrom[i] != chrom[j]:
                m[j, i] = m[i, j] = 10e10
            # splits names as is it stores the info about BL8778
            # clustered names store this info in each second element after double split starting from [1]
            elif bool(set(re.split("[#,]", names[i])[1::2]) & set(re.split("[#,]", names[j])[1::2])) == True:
                m[j, i] = m[i, j] = 10e10
            # forbid the same names... ever second starting from 0 (good agains cluster)
            elif bool(set(re.split("[#|]", names[i])[0::2]) & set(re.split("[#|]", names[j])[0::2])) == True:
                m[j, i] = m[i, j] = 10e10
            elif position[i] == position[j]:
                m[j, i] = m[i, j] = 10e10
            else:
                m[j, i] = m[i, j] = abs(position[i] - position[j])
    return m


def collapse_matrix(m, names, position, chrom, type=None, RNApos=None):
    '''
    :param m: old distance matrix
    :param names: names of the old distance matrix
    :param position: postion of the old distance matrix
    :param chrom: chrom numbers of the old distance matrix
    :return: new distance matrix (clustered and collapsed)
    '''
    # finding the minimum coordinates in the matrix
    # they have to be sorted by their size to pop out list elements properly
    if type == "ini":
        # the initial matrix needs ini_minimum
        high, low = sorted(ini_minimum(m, RNApos), reverse=True)
    else:
        high, low = sorted(find_minimum(m), reverse=True)
    # use minimum coordinates to cluster the lists
    new_names = collapse_list(names, high, low)
    new_pos = collapse_list(position, high, low)
    new_chrom = collapse_list(chrom, high, low)
    # generate default matrix
    raw_mat = generate_matrix(new_names)
    # compute collapsed, new, clustered distance matrix...overwritting the raw_mat
    collapsed_mat = distance_matrix(raw_mat, new_chrom, new_pos, new_names)
    return collapsed_mat


def find_path(matrix_names, best_cluster):
    '''
    Function to find the best path in a matrix, starting from the best initial cluster
    :param matrix_names: distance matrix
    :return: the position values of how the genes of the best path appear in the original "names"-list
    '''
    path = ""
    gene_value_list = list()
    for cluster in matrix_names:
        cluster_size = len(re.split("[,|]", cluster)[1::2])
        if cluster_size > 1:
            gene_value_list.extend(re.split("[,|]", cluster)[1::2])
        if re.search(best_cluster, cluster):
            path = cluster
    # genes are stored like "SLC17A5#8.14|PLG#5.9|HIST1H3F#11.22", wherebey the number after the dot holds the position
    path_position = re.split("[,|]", path)[1::2]
    return path_position


#def get_all_possible_path(matrix_names):
#    '''
#    The function is similar to the find_path function. But instead of retruning only the best cluster
#    it retuns the position of all genes that ended up in clusters of 2 or more genes.
#    :param matrix_names: distance matrix
#    :return: position values of all the genes that formed a cluster >= 2 genes. (position of how they appear in the names-list)
#    '''
#    gene_pos_list = list()
#    for cluster in matrix_names:
#        cluster = re.split("[,|]", cluster)[1::2]
#        if len(cluster) > 1:
#            gene_pos_list.extend(cluster)
#    return gene_pos_list

def get_all_possible_path(matrix_names, org_n,org_pos,org_str,org_chr,org_amp):
    '''
    The function is similar to the find_path function. But instead of retruning only the best cluster
    it retuns the position of all genes that ended up in clusters of 2 or more genes.
    :param matrix_names: distance matrix
    :return: position values of all the genes that formed a cluster >= 2 genes. (position of how they appear in the names-list)
    '''
    string_list = list(" ")
    for cluster in matrix_names:
        cluster = re.split("[,|]", cluster)[1::2]
        if len(cluster) > 1:
            for pos in cluster:
                gene = int(pos)
                gene_string = org_n[gene].split("#")[0]+",Pos:"+str(org_pos[gene])+\
                      ",Strand:"+org_str[gene]+",Chr:"+str(org_chr[gene])+",BL:"+org_amp[gene]
                string_list.append(gene_string)
            string_list.append("\n")
    string_path = "\t".join(string_list)

    return string_path

def path_goodness(path_values, org_n,org_pos,org_str,org_chr,org_amp):
    string_list=[]
    distance_list=[]
    if len(path_values) > 0:
        for val in path_values:
            gene=int(val)
            gene_string = org_n[gene].split("#")[0]+",Pos:"+str(org_pos[gene])+\
                      ",Strand:"+org_str[gene]+",Chr:"+str(org_chr[gene])+",BL:"+org_amp[gene]
            string_list.append(gene_string)
            distance_list.append(int(org_pos[gene]))
        goodness=statistics.stdev(distance_list)
        string_path= "\t".join(string_list)
    else:
        goodness= 10e10
        string_path="No Syntheny found!"
    return round(goodness,1), string_path


####################################################################################################################

# EVOKER FUNCTION !!!!!!
# allows to just import the this function from the script into another script for doing all the clustering

def search_syntheny(names, chrom, position, strand, amphiIdenti, RNApos):
    '''
    starting function which evokes all other function and returns the final clustering output
    :param names: original names vector
    :param chrom: original chromosome vector (as int)
    :param position: original position vector (as int)
    :return: position values of the longest-best path through the final distance matrix (those can then be used to
    parse the original list according to ones needs !!!!!!
    '''
    # store important lists from the beginning to keep them from changing
    orig_str = tuple(strand)
    orig_amp = tuple(amphiIdenti)
    orig_n = tuple(names)
    orig_chr = tuple(chrom)
    orig_pos = tuple(position)
    # initialize the first dist. matrix with the input lists
    minimum, best_cluster, iter_mat = Initialize(names, chrom, position, RNApos)
    # output is the first minimum, the best name-cluser and the first iteration-matrix
    # than calculate the goodness for the best_cluster
    seed_path= find_path(iter_mat.names, best_cluster)
    seed_goodness,first_pathway= path_goodness(seed_path,orig_n,orig_pos,orig_str,orig_chr,orig_amp)
    while minimum != "No_minimum":
        iter_mat = collapse_matrix(iter_mat, names, position, chrom)
        minimum = find_minimum(iter_mat)
    # take the names of the last matrix .. as a list
    # find_path returns a list of position values to access the original lists
    try:
        path_values = find_path(iter_mat.names, best_cluster)
        goodness,path=path_goodness(path_values,orig_n,orig_pos,orig_str,orig_chr,orig_amp)
        all_path = get_all_possible_path(iter_mat.names,orig_n,orig_pos,orig_str,orig_chr,orig_amp)

        return goodness,seed_goodness, path, len(path_values), all_path
    except IndexError:
        return ["No_unique_path"]


