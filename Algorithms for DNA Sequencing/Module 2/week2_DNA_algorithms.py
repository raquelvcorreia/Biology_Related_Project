import bisect
from bm_preproc import BoyerMoore


#implementing Boyer-Moore





def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
           p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i <len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if not p[j] == t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences


class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1): # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index): # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits



def query_index(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i+k : i+len(p)]:
            offsets.append(i)
    return offsets


t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'

index = Index(t, 2)
print(query_index(p, t, index))


def naive_hamming(p, t, max_dist):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        nmm = 0
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters
                nmm += 1
                if nmm > max_dist:
                    break
        if nmm <= max_dist:
            occurrences.append(i)  # approx match
    return occurrences


def approximate_match(p, t, n):
        segment_length = round(len(p)/(n+1))
        all_matches = set()
        for i in range(n+1):
            start = i*segment_length
            end = min((i+1)*segment_length, len(p))
            p_bm = boyer_moore(p[start:end], p_bm, t)

            for m in matches:
                if m < start or m-start+len(p) > len(t):
                    continue

                mismatches = 0
                for j in range(0, start):
                    if not p[j] == t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break
                for j in range(end, len(p)):
                    if not p[j] == t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break

                if mismatches <= n:
                    all_matches.add(m-start)

        return list(all_matches)




def naive_with_counts(p, t):
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        num_alignments += 1
        match = True
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i + j] != p[j]:  # compare characters\
                match = False  # mismatch; reject alignment
                break
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences, num_alignments, num_character_comparisons

## test 1
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print('Test 1 naive with counts ', occurrences, num_alignments, num_character_comparisons)


## test 2
p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print('Test 2 naive with counts ', occurrences, num_alignments, num_character_comparisons)



def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences




def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(p)-1, -1, -1):
            num_character_comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_alignments, num_character_comparisons


def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

chrm1_human = read_genome('chr1.GRCh38.excerpt.fasta')


def match_approximate_ph(p, t, n):
    """ Find all approximate occurrences of p in t with up to n mismatches.
            p=pattern, t=text, n = number of maximum mismatches
            (insertions and deletions are not allowed ) using class Index"""

    segment_length = round((len(p) / (n + 1)))
    all_matches = set()
    p_idx = Index(t, 8)
    idx_hits = 0
    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))
        matches = p_idx.query(p[start:end])

        # check if the whole p matches
        for m in matches:
            idx_hits +=1
            if m < start or m - start + len(p) > len(t):
                continue

            mismatches = 0

            for j in range(0, start):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m - start)

    return list(all_matches), idx_hits




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


class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i + self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits



def match_approximate_ph_subseq(p, t, n, ival):
    """ Find all approximate occurrences of p in t with up to n mismatches.
            using subsequences (SubseqIndex) - taking only the nth character
            p=pattern, t=text, n = number of maximum mismatches
            (insertions and deletions are not allowed ) using class Index"""

    segment_length = round((len(p) / (n + 1)))
    all_matches = set()
    p_idx = SubseqIndex(t, 8, ival)
    idx_hits = 0
    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))
        matches = p_idx.query(p[start:end])

        # check if the whole p matches
        for m in matches:
            idx_hits +=1
            if m < start or m - start + len(p) > len(t):
                continue

            mismatches = 0

            for j in range(0, start):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m - start)

    return list(all_matches), idx_hits





## Test 1
p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print('Test 1 BM with counts ', occurrences, num_alignments, num_character_comparisons)

# Test 2
p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print('Test 2 BM with counts ', occurrences, num_alignments, num_character_comparisons)


## Question one, two
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = chrm1_human
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print('Q1, Q2 naive with counts ', occurrences, num_alignments, num_character_comparisons)

## Question three
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
t = chrm1_human
p_bm = BoyerMoore(p)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print('Q3 BM with counts ', occurrences, num_alignments, num_character_comparisons)


## Question 4, 5
p_4 = 'GGCGCGGTGGCTCACGCCTGTAAT'
print('Q4 ', len((match_approximate_ph(p_4, t, 2)[0])))
print('Q5 ', match_approximate_ph(p_4, t, 2)[1])

print('Occurrences using the naive_2mm function ', len(naive_2mm(p_4, t)))
print('Number of hits using the Boyer-Moore function ', len(boyer_moore(p_4, BoyerMoore(p_4), t)))
print('Number of hits using the naive function ', len(naive(p_4, t)))


## Question 6
print('Q6 ', match_approximate_ph_subseq(p_4, t, 2, 3)[1])