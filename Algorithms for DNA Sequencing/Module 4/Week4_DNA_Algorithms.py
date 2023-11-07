import itertools
from itertools import permutations
from collections import defaultdict


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

#shortest common super force by brute force
def scs_original(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

#new scs version where all the shortest common superstrings are kept variable shortest_sup_list
#there can be more than one shortest supperstring this way we can get a list of those
def scs(ss):
    """ Returns shortest common superstring of given
            strings, which must be the same length """
    shortest_sup = None
    sup_list = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]   # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length = 1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) <= len(shortest_sup):
            sup_list.append(sup)
            if shortest_sup is None or len(sup) < len(shortest_sup):
                shortest_sup = sup  # found shorter superstring

    min_length = min(len(sup_) for sup_ in sup_list)

    shortest_sup_list = []
    for i in sup_list:
        if len(i) == min_length:
            shortest_sup_list.append(i)


    return shortest_sup, shortest_sup_list # return shortest


#print(scs(['ACGGTACGAGC', 'GAGCTTCGGA', 'GACACGG']))

print('Q1 ', len(scs(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])[0]))

print('Q2 ', scs(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']))

print('Test Q2 ', scs(['ABC', 'BCA', 'CAB']))
print('Test Q2 ', scs(['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']))


#greedy shortest common superstring
#although a lot faster than the scs it does not always report the shortest common string

def pick_maximal_overlap_notoptimal(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length = k) #in case of ties basically the algorith will just pick the first one it encounters

        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen

    return reada, readb, best_olen

def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    # Make index
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i + k]].add(read)

    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                olen = overlap(r, o, k)
                if olen > best_olen:
                    reada, readb = r, o
                    best_olen = olen

    return reada, readb, best_olen

    return reada, readb, best_olen


def greedy_scs(reads,k):
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads,k)
    return ''.join(reads)

#in this example the greedy_scs does not effectly give out the scs
print(greedy_scs(['ABCD', 'CDBC', 'BCDA'], 1))
print(scs(['ABCD', 'CDBC', 'BCDA']))


def de_bruijn_ize (st,k):
    edges = []
    nodes = set()
    for i in range(len(st) - k + 1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1 : i+k])
    return nodes, edges

nodes, edges = de_bruijn_ize('ACGCGTCG', 3)
print(nodes, edges)


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


seqs, _ = read_fastq('ads1_week4_reads.fq')

for k in range(100, 1, -1):
    genome = greedy_scs(seqs, k)
    if len(genome) == 15894:
        print(genome.count('A'))
        print(genome.count('T'))
        print(genome)
        break