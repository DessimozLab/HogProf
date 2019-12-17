import numpy
import collections
from lark import Lark, Transformer, ParseError
from tables import dtype_from_descr
from .tablefmt import LocusTable
import logging


logger = logging.getLogger(__name__)

"""This package is intended to parse the darwin locus structure
and create a numpy recarray out of it. """


Exon = collections.namedtuple('Exon', ['start', 'end', 'strand'])

grammar = '''?locus : join | complement | complement_join | location
             join  : "join" "(" (complement | location ) ("," (complement | location ))+ ")"
             complement : "complement" "(" location ")"
             complement_join : "complement" "(" "join" "(" location ("," location)+ ")" ")"
             location : pos [".." pos ] | "FromElsewhere" "('" _SEQID "'," pos [".." pos] ")" 
             ?pos : num | "Before" "(" num ")" | "After" "(" num ")"
             ?num : NUMBER      -> number
             _SEQID: /[A-Za-z0-9._-]+/
             
             %import common.NUMBER
             %import common.WS
             %ignore WS'''


class LocusTransformer(Transformer):
    def number(self, vals):
        return int(vals[0])

    def location(self, value):
        return Exon(value[0], value[1] if len(value) > 1 else value[0], 1)

    def complement(self, value):
        rev = [e._replace(strand=-1*e.strand) for e in value]
        if len(rev) == 1:
            return rev[0]
        else:
            return rev

    def complement_join(self, value):
        return self.complement(value)

    def join(self, values):
        return values


class LocusParser(object):
    def __init__(self):
        self.parser = Lark(grammar, start='locus')
        self.locus_transformer = LocusTransformer()
        self.dtype = dtype_from_descr(LocusTable)

    def parse(self, locus_string, entry_nr=0):
        try:
            tree = self.parser.parse(locus_string)
        except ParseError as e:
            raise ValueError("cannot parse '{}' locus string".format(locus_string))
        data = self.locus_transformer.transform(tree)
        nr_exons = 1 if isinstance(data, Exon) else len(data)
        locus_data = numpy.empty(nr_exons, dtype=self.dtype)
        locus_data[['Start', 'End', 'Strand']] = data
        locus_data['EntryNr'] = entry_nr
        return locus_data
