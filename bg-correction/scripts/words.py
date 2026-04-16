from itertools import product

kmer2index = {} # map kmer to index
kmers = []  # map index to kmer
kmer_rcs = {} # map kmer to its reverse complement kmer
n_words = 256

comp_lut = {"A":["T"],"C":["G"],"G":["C","T"],"T":["A","G"]}
def get_rc(s):
    last_comps = [""]
    for c in s:
        comps = []
        for comp in last_comps:
            for cc in comp_lut[c]:
                comps.append(comp + cc)
        last_comps = comps
    return [comp[::-1] for comp in comps]

def init(k=4):
    for i,kmer in enumerate(list(product(*["ACGT"]*k))):
        kmer = "".join(kmer)
        kmer2index[kmer] = i
        kmers.append(kmer)

    for kmer in kmers:
        kmer_rcs[kmer] = get_rc(kmer)
    global n_words
    n_words = len(kmers)
