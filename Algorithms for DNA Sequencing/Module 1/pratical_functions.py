import matplotlib.pyplot as plt

def longest_common_prefix(s1, s2):
    # identify matching patches of sequences
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i += 1
    return s1[:i]


print(longest_common_prefix('ACCATGT', 'ACCAGAC'))


# get the complement in DNA
# create a dictionary
def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}  #when the confidence on a read is below a certain threshold the program will give an N
    t = ''
    for base in s:
        t = complement[base] + t
    return t


print(reverse_complement('ACCGTCG'))


def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


genome = read_genome('lambda_virus.fa')
print(genome[:100])
print(len(genome))

# what is the frequency of each base? create a dictionary to keep track
counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
for base in genome:
    counts[base] += 1

print(counts)

# this can also be done using a python module
import collections

print(collections.Counter(genome))


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


seqs, quals = read_fastq('SRR835775_1.first1000.fastq')

print(seqs[:5])
print(quals[:5])


def phred33_to_q(qual):
    # conver the ascii encoding to the actual quality value
    # ord() returns the unicode from a given character
    return ord(qual) - 33


print(phred33_to_q('#'))


def create_hist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33_to_q(phred)
            hist[q] += 1
    return hist


h = create_hist(quals)



def find_gc_by_pos(reads):
    gc = [0] * 100
    totals = [0] * 100

    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1

            totals[i] += 1
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc


gc = find_gc_by_pos(seqs)

# plt.plot(range(len(gc)), gc)
# plt.show()

# get the amount of bases of each type
count = collections.Counter()
for seq in seqs:
    count.update(seq)


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters\
                match = False  # mismatch; reject alignment
                break
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences


genome = read_genome('phix.fa')

import random


# generate artificial reads from a genome
def generate_reads(genome, numreads, readlen):
    '''Generate reads from randomm positions in a given genome.'''

    reads = []
    for _ in range(numreads):
        start = random.randint(0, len(genome) - readlen) - 1
        reads.append(genome[start: start + readlen])
    return reads

reads = generate_reads(genome, 100, 100)
num_matches = 0
for r in reads:
    matches = naive(r, genome)
    if len(matches) > 0:
        num_matches += 1

print('%d / %d reads matched exactly!' % (num_matches, len(reads)))

phix_reads, _ = read_fastq('ERR266411_1.first1000.fastq')
num_matches_1 = 0
for r in phix_reads:
    matches = naive(r, genome)
    if len(matches) > 0:
        num_matches_1 += 1

print('%d / %d reads matched exactly!' % (num_matches_1, len(reads)))
# not those many reads match, it could be a problem with the quality of the reads (we can read just for instance the first 30 bases)
# but it can also be due to the fact that genomes are double-stranded and we are
# only reading on direction (also try to match the complement

num_matches_2 = 0
for r in phix_reads:
    r = r[:30]  ##only consider the beginner of the read
    matches = naive(r, genome)
    matches.extend(naive(reverse_complement(r), genome)) ##adds to the new list to list
    if len(matches) > 0:
        num_matches_2 += 1
print('%d / %d reads matched exactly!' % (num_matches_2, len(reads))) #better but still not 100 out of 100. Maybe exact match are not the best option

def naive_w_complement(p, t):
    p_reverse = reverse_complement(p)
    occurrences = []

    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters
                match = False

        if match:
            occurrences.append(i)  # all chars matched; record

        if p != p_reverse:
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i + j] != p_reverse[j]:
                    match = False
                    break
            if match:
                occurrences.append(i)

    return occurrences


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        missmatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters\
                missmatches +=1
                if missmatches >2:
                    break

        if missmatches <= 2:
            occurrences.append(i)

    return occurrences


p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
test1 = naive_w_complement(p, t)

phix_genome = read_genome('phix.fa')
lambda_genome = read_genome('lambda_virus.fa')
test2 = naive_w_complement('ATTA', phix_genome)
test3 = naive_w_complement('AGTCGA', lambda_genome)
test4 = naive('TTAA', lambda_genome)
test5 = naive_2mm('TTCAAGCC', lambda_genome)
test6 = naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT')
test7 = naive_2mm('AGGAGGTT', lambda_genome)

seqs_test, quals_test = read_fastq('ERR037900_1.first1000.fastq')

overall_qual = create_hist(quals_test)

#plt.bar(range(len(overall_qual)), overall_qual)
#plt.show()

def lowest_quality_base(qs):
    total = [0] * len(qs[0])

    for q in qs:
        for i, phred in enumerate(q):
            total[i] += phred33_to_q(phred)


    return total.index(min(total))


