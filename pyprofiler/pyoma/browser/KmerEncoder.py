#!/usr/bin/env python
'''
    DNA / AA to integer packer. i.e., base 5/21 <-> base 10.
    Heavily adapted, but based on https://gist.github.com/bitsandbooks/2649444

    NOTE: any sequences that are being decoded should be sanitised first.

    -- Alex Warwick Vesztrocy, May-June 2017
'''
import numpy as np


# "digits"
DIGITS_AA = np.fromstring('ACDEFGHIKLMNPQRSTVWXY', dtype='S1')
DIGITS_DNA = np.fromstring('ACGTX', dtype='S1')


class KmerEncoder(object):
    def __init__(self, k, is_protein=True):
        '''
            Initialise the kmer converter. k is the kmer length.
            If is_dna=True then DNA else AA.
        '''
        self.digits = DIGITS_AA if is_protein else DIGITS_DNA
        self.k = int(k)  # Cast incase np-type for n below.
        self.max = (len(self.digits) ** self.k) - 1
        self.n = self.decode(self.digits[-1] * self.k) + 1
        self._prot = np.zeros((self.k,), dtype='S1')

    def __len__(self):
        '''
            Return the maximum integer-representation of the kmer length
            in this converter.
        '''
        return self.n

    def encode(self, seq):
        '''
            Encode integer kmer in protein chars.
        '''
        if seq <= self.max:
            self._prot[:] = self.digits[0]
            i = -1
            while seq > 0:
                self._prot[i] = self.digits[seq % self.digits.size]
                seq //= self.digits.size
                i -= 1
            return self._prot.tostring()
        else:
            raise ValueError('{} Larger than largest kmer of size {}'
                             .format(seq, self.k))

    def decode(self, seq):
        '''
            Decode a protein kmer -> integer. NOTE: sanitisation to a byte
            string required first.
        '''
        x = 0
        for digit in seq[:self.k].decode('ascii'):
            x = (x * self.digits.size) + np.searchsorted(self.digits, digit)
        return x

    def decompose(self, seq):
        '''
            Decompose a sequence into counts of its constituent (decoded) kmers.
        '''
        for i in range(len(seq) - self.k + 1):
            yield self.decode(seq[i:(i + self.k)])
