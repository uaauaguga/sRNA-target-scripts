#!/usr/bin/env python
from itertools import product
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('get 6mer energy')

def get_rc(s):
    comp_lut = {"A":["T"],"C":["G"],"G":["C","T"],"T":["A","G"]}
    last_comps = [""]
    for c in s:
        comps = []
        for comp in last_comps:
            for cc in comp_lut[c]:
                comps.append(comp + cc)
        last_comps = comps
    return [comp[::-1] for comp in comps]


def main():
    energies = {}
    k = 6
    cmd = ["RNAcofold"]
    fout = open("hexamer-energy.txt","w")
    for kmer in list(product(*["ACGT"]*k)):
        kmer = "".join(kmer)
        kmer_rcs = get_rc(kmer)
        for kmer_rc in kmer_rcs:
            if kmer > kmer_rc:
                continue
            if (kmer, kmer_rc) in energies:
                continue
            sequence = f"N{kmer}N&N{kmer_rc}N"
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stdin=subprocess.PIPE)
            lines = proc.communicate(sequence.encode())[0].decode()
            line = lines.strip().split("\n")[-1]
            energy = line[line.rfind("(")+1:line.rfind(")")]
            print(kmer,kmer_rc,energy,sep="\t",file=fout)
            fout.flush()
            energies[(kmer,kmer_rc)] = float(energy)
            proc.wait()
    fout.close()    
             

if __name__ == "__main__":
    main()
