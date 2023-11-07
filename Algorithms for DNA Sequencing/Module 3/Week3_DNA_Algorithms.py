
from itertools import permutations
from collections import defaultdict
# dynamic programming for edit distance


# initial version which is slower and does not involve dynamic programming
from datetime import time


def edit_dist_recursive(x,y):
    # slow implementation --> not to use
    if len(x) == 0:
        return len(y)
    elif len(y) == 0:
        return len(x)
    else:
        dist_horiz = edit_dist_recursive(x[:-1], y) + 1
        dist_vert = edit_dist_recursive(x, y[:-1]) + 1
        if x[-1] == y[-1]:
            dist_diag = edit_dist_recursive(x[:-1], y[:-1])
        else:
            dist_diag = edit_dist_recursive(x[:-1], y[:-1]) + 1

        return min(dist_horiz, dist_vert, dist_diag)


def edit_dist(x,y):
    matrix_d = []
    # initializing the matrix
    for i in range(len(x) +1):
        matrix_d.append([0]* (len(y) + 1))
    #maximum distances is one string is 0 and the other on is i to the length of non 0 string/ initialize the first crow and column
    for i in range(len(x) + 1):
        matrix_d [i][0] = i
    for i in range(len(y) + 1):
        matrix_d [0][i] = i
    #fill in the rest of the matrix
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            dist_horiz = matrix_d[i][j-1] + 1
            dist_vert = matrix_d[i-1][j] + 1
            if x[i-1] == y[j-1]:
                dist_diag = matrix_d[i-1][j-1]
            else:
                dist_diag = matrix_d[i-1][j-1] + 1

            matrix_d [i][j] = min(dist_horiz, dist_vert, dist_diag)

    return matrix_d[-1][-1], min(matrix_d[-1])

#testing and checking how much faster the edit_dist function is
#x = 'shake spea'
#y = 'shakespear'
#start = timer()
#print(edit_dist_recursive(x,y))
#end = timer()
#print('The time of execution for recursive algorithm ', (end-start)*10**3, 'ms')
#start = timer()
#print(edit_dist(x,y))
#end = timer()
#print('The time of execution for edit_dist algorithm ', (end-start)*10**3, 'ms')


## global alignment, different mismatches are penalized differently


alphabet = ['A', 'C', 'G', 'T']
score = [[0, 4, 2, 4, 8],
         [4, 0, 4, 2, 8],
         [2, 4, 0, 4, 8],
         [4, 2, 4, 0, 8],
         [8, 8, 8, 8, 8]]

def global_alignment(x,y):
    matrix_d = []
    # initializing the matrix
    for i in range(len(x) +1):
        matrix_d.append([0]* (len(y) + 1))

    for i in range(1, len(x) + 1):
        matrix_d [i][0] = matrix_d[i-1][0] + score[alphabet.index(x[i-1])][-1]
    for i in range(1, len(y) + 1):
        matrix_d [0][i] = matrix_d[0][i-1] + score[-1][alphabet.index(y[i-1])]

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            dist_horiz = matrix_d[i][j-1] + score[-1][alphabet.index(y[j-1])]
            dist_vert = matrix_d[i-1][j] + score[alphabet.index(x[j-1])][-1]
            if x[i-1] == y[j-1]:
                dist_diag = matrix_d[i-1][j-1]
            else:
                dist_diag = matrix_d[i-1][j-1] + score[alphabet.index(x[j-1])][alphabet.index(y[j-1])]

            matrix_d [i][j] = min(dist_horiz, dist_vert, dist_diag)

    return matrix_d[-1][-1]


x = 'TACCAGATTCGA'
y = 'TACCAGATTCG'

#print(global_alignment(x, y))

def overlap(a, b, min_length = 3):
    '''Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0'''
    start = 0  #start wll the way at the left
    while True:
        start = a.find(b[:min_length], start) #look for b sufix in a
        if start == -1: # no more occurrences to the right
            return 0
        #found occurrences: check full suffiz/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1 #move just past the previous match

print(overlap('TTACGT', 'CGTACCGT'))


c = list(permutations([1,2,3],3))



def naive_overlap_map(reads, k):
    olaps = {}
    for a,b in permutations(reads,2):
        olen = overlap(a,b, min_length=k)
        if olen > 0:
            olaps[(a,b)]  = olen
    return olaps


reads = ['ACGGATGATC', 'GATCAAGT', 'TTCACGGA']
#print(naive_overlap_map(reads, 3))

test_p = 'GCGTATGC'
test_t = 'TATTGGCTATACGGTT'

def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

chrm1_human = read_genome('../Module 2/chr1.GRCh38.excerpt.fasta')
#print(chrm1_human[:5])

def approx_match(x,y):
    matrix_d = []
    # initializing the matrix
    for i in range(len(x) +1):
        matrix_d.append([0]* (len(y) + 1))
    #maximum distances is one string is 0 and the other on is i to the length of non 0 string/ initialize the first crow and column
    for i in range(len(x) + 1):
        matrix_d [i][0] = i
    for i in range(len(y) + 1):
        matrix_d [0][i] = 0
    #fill in the rest of the matrix
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            dist_horiz = matrix_d[i][j-1] + 1
            dist_vert = matrix_d[i-1][j] + 1
            if x[i-1] == y[j-1]:
                dist_diag = matrix_d[i-1][j-1]
            else:
                dist_diag = matrix_d[i-1][j-1] + 1

            matrix_d [i][j] = min(dist_horiz, dist_vert, dist_diag)

    return matrix_d[-1][-1], min(matrix_d[-1])





q1 = approx_match('GCTGATCGATCGTACG',chrm1_human)
q2 = approx_match('GATTTACCAGATTGAG',chrm1_human)
test_q1 = approx_match('GCGTATGC','TATTGGCTATACGGTT')
print(q1,q2)
#print(test_q1)

def read_fastq(filename):
    sequence = []
    qualities = []
    with open(filename) as fh:
        while True:
            # will read every 4 lines until the end of the file because each 4 lines corresponds to one read
            fh.readline()  # we read it but dont strore is because we don't need it, is the first line that has information on the read
            seq = fh.readline().rstrip()  # this one contains the actual read so we add it to seq
            fh.readline()  # the '+' that we dont  need
            qual = fh.readline().rstrip()  # quality information that we want to keep
            if len(seq) == 0:
                break
            sequence.append(seq)
            qualities.append(qual)
    return sequence, qualities




def overlap_map_q3(reads, k):
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i + k]].add(read)

    # Make graph
    graph = defaultdict(set)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                if overlap(r, o, k):
                    graph[r].add(o)
    edges = 0
    for read in graph:
        edges += len(graph[read])
    return (edges, len(graph))




reads = read_fastq('ERR266411_1.for_asm.fastq')[0]
reads_2 = read_fastq('ERR266411_1.for_asm_v2.fastq')[0]
reads_subset = reads[:10]

print(overlap_map_q3(reads_subset, 30))
print(overlap_map_q3(reads, 30))
