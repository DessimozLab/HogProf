#1900000000

from bloom_filter2 import BloomFilter
import multiprocessing as mp
import functools
# make a bloom filter with all of the interactions
bloom = BloomFilter(max_elements=10**10, error_rate=0.001, filename=('./bloom.bin', -1))

def pairgen(infile, startline= 0, reportlines = 100000000, verbose =True):
    with open(infile,'r') as inlines:
        for i,l in enumerate(inlines):
            if i > startline:
                words = l.split()
                if i % reportlines == 0 and verbose:
                    print(i)
                yield words[0]+words[1]
        print('done')

gen = pairgen('/scratch/dmoi/datasets/STRING/protein.physical.links.detailed.v11.5.txt')
[bloom.add(p) for p in gen]
print('DONE')
