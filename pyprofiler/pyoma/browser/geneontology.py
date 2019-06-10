from builtins import int, bytes, str
import collections
import csv
import logging
import math
import re
import numpy

"""
IMPORTANT NOTE:
---------------
This module has been copied from the dessimoz zoo library
directly. If you want to add functionality to this module,
make sure it is also integrated into the zoo. At the moment
we don't want to depend with pyoma on the zoo library as it
has many dependencies that are difficult to maintain.


Gene ontology module defining classes and methods to parse
and navigate the gene ontology DAG as well as to parse GO
annotations.


:author: Adrian Altenhoff
:institute: ETH Zurich
"""

NUM_ONT = 3
NUM_GO_ID_DIGITS = 7
UP_RELS = frozenset(['is_a', 'part_of'])
_REV_RELS = {'is_a': 'can_be', 'part_of': 'has_part'}


def reverse_name_of_rels(rels):
    down = frozenset([_REV_RELS[z] for z in rels])
    return down


def validate_go_id(term):
    if isinstance(term, (int, numpy.integer)):
        return int(term)

    term = term.strip()
    if term.startswith('GO:'):
        digits = term[3:]
    else:
        digits = term
    if not digits.isdigit() or len(digits) > NUM_GO_ID_DIGITS:
        raise ValueError("GO ID {} is not a valid go term".format(term))
    return int(digits)


class GOAspect(object):
    aspects = dict(molecular_function=0, biological_process=1, cellular_component=2)
    aspect2char = {0: 'F', 1: 'P', 2: 'C'}

    @classmethod
    def from_string(cls, aspect):
        return cls.aspects[aspect]

    @classmethod
    def to_string(cls, aspectnr):
        for o, i in cls.aspects.items():
            if i == aspectnr:
                return o
        raise KeyError('aspect number not found: ' + str(aspectnr))

    @classmethod
    def to_char(cls, aspectnr):
        # Converts an encoded aspect to the character required for GOA files
        return cls.aspect2char[aspectnr]


class GOterm(object):
    """A class representing a single Gene Ontology term.

    This class can serve as a factory for the OntologyParser. For that,
    pass it as a factory on 'Term'. """

    def __init__(self, stanza):
        self.id = validate_go_id(stanza['id'][0])
        self.name = ' '.join(stanza['name'])
        self.definition = ' '.join(stanza['def'])
        self.aspect = GOAspect.from_string(' '.join(stanza['namespace']))
        self.is_a = [validate_go_id(parent) for parent in stanza['is_a']]
        for rel in stanza['relationship']:
            reltype, partner = rel.strip().split()
            if not reltype in self.__dict__.keys():
                self.__dict__[reltype] = list()
            self.__dict__[reltype].append(validate_go_id(partner))

    def replace_parentnames_by_refs(self, ont):
        for rel in [('is_a', 'can_be'), ('part_of', 'has_part')]:
            if rel[0] in self.__dict__.keys():
                for i, parent_id in enumerate(self.__dict__[rel[0]]):
                    parent_obj = ont[parent_id]
                    self.__dict__[rel[0]][i] = parent_obj
                    parent_obj._add_relation(self, rel[1])

    def _add_relation(self, term, rel):
        if rel not in self.__dict__.keys():
            self.__dict__[rel] = list()
        self.__dict__[rel].append(term)

    def get_parents(self, rels=None):
        """iterate over the direct parent GO terms.

        by default "is_a" and "part_of" relations are followed. This can be overwritten
        with the `rels`.

        :param rels: a set of relations to follow."""
        if rels is None:
            rels = UP_RELS
        for rel in rels:
            try:
                for term in getattr(self, rel):
                    yield term
            except AttributeError:
                pass

    def __str__(self):
        fmt = "GO:{{0:0{}d}}".format(NUM_GO_ID_DIGITS)
        return fmt.format(self.id)


class AbstractParser(object):
    def __init__(self, fp):
        self.close_fp = False
        if isinstance(fp, str):
            if fp.endswith('.gz'):
                from gzip import GzipFile
                fp = GzipFile(fp, 'rb')
            else:
                fp = open(fp, 'r')
            self.close_fp = True
        self.fp = fp
        self._read_headers()

    def _read_headers(self):
        pass

    def close_if_opened(self):
        if self.close_fp:
            self.fp.close()


class OntologyParser(AbstractParser):
    """A general purpose Ontolgoy parser

    Any ontology in the OBO format can be parsed with this object. The
    stanzas are converted to objects using the factories passed in the
    initializer."""
    def __init__(self, fp, factories=None):
        """creates an ontology parser

        :param fp: a filehandle or path to file (either plaintext or
            gzipped) containing the ontology.
        :param factories: a dictionary containing per stanza class
            (e.g. [Term]) a factory that returns an object from the
            data. The data is passed as dict to the factory"""
        if not factories:
            factories = dict(Term=GOterm)
        super(OntologyParser, self).__init__(fp)
        self.factories = factories
        self.tag_value_pair_re = re.compile(r"\s*(?P<tag>[^:]+):\s*(?P<value>[^!]*)")
        self.stanza_name_re = re.compile(r"\[(?P<name>[^]]*)\]")

    def stanzas(self):
        """iterates over the stanzas in the ontology yielding 
        objects according to the factory parameter provided 
        in the constructor."""
        curStanza = None
        for line in self.fp:
            line = line.strip()
            if not line or line[0] == '!':
                continue

            # check whether a new stanza starts
            match = self.stanza_name_re.match(line)
            if match is not None:
                obj = self._create_obj_from_stanza(curStanza)
                if obj is not None:
                    yield obj
                curStanza = collections.defaultdict(list)
                curStanza['_name'] = match.group("name")
            elif curStanza is not None:
                # it has to be value-pair line. add it to the stanza
                match = self.tag_value_pair_re.match(line)
                curStanza[match.group('tag')].append(match.group('value'))
        obj = self._create_obj_from_stanza(curStanza)
        if obj is not None:
            yield obj
        self.close_if_opened()

    def _create_obj_from_stanza(self, stanza):
        """method which creates the appropriate object from a given
        stanza or returns None, e.g. for stanza types without a provided
        factory or obsolete terms or ..."""
        res = None
        if stanza is not None:
            try:
                factory = self.factories[stanza['_name']]
            except KeyError:
                # no factory for this stanza type. ignore
                pass
            else:
                if ((not "is_obsolete" in stanza) or
                        (not 'true' in stanza['is_obsolete'])):
                    res = factory(stanza)
        return res

    def __iter__(self):
        return self.stanzas()


class GeneOntology(object):
    """The GeneOntology object contains the whole ontology in an internal
    format and allows to traverse it efficiently."""
    def __init__(self, parser, rels=None):
        if rels is None:
            rels = UP_RELS
        if not isinstance(parser, OntologyParser):
            raise Exception('requires an OntologyParser instance')
        self.parser = parser
        self.up_rels = UP_RELS.intersection(rels)
        self.down_rels = reverse_name_of_rels(self.up_rels)

    def parse(self):
        """parse the ontology data.

        This method should be called after instantiation and before any traversal"""
        self.terms = dict()
        for cur_term in self.parser:
            self.terms[cur_term.id] = cur_term

        # replace parents nrs by the references to the GO terms objects
        # this can be only done once all the terms have been created
        for term in self.terms.values():
            term.replace_parentnames_by_refs(self.terms)

    def ensure_term(self, term):
        """returns the term object associated with term. if term is already
        a GOterm object, it is simply return. term_ids can be either numeric
        ids of existing GO-terms or propper GO-terms ids, i.e. GO:000002

        :param term: a term_id or a GOterm object"""
        if isinstance(term, GOterm):
            return term
        else:
            return self.term_by_id(term)

    def term_by_id(self, term_id):
        """Returns the term object associated with term_id.

        :param term_id: a GO-term number or a GO-term id (GO:008150)."""
        try:
            term = self.terms[validate_go_id(term_id)]
            return term
        except KeyError:
            raise ValueError(str(term_id) + ' is an invalid GO term.')

    def get_superterms_incl_queryterm(self, term, max_steps=-1):
        """returns a set with all the superterms of a query term.

        :param max_steps: The search can be limited to contain only
            terms that are at most 'max_steps' upwards. If set to -1, no
            limit is applied and the search goes up to the root."""
        term = self.ensure_term(term)
        return self._traverseGraph(term, max_steps, self.up_rels)

    def get_subterms(self, term, max_steps=-1):
        term = self.ensure_term(term)
        return self._traverseGraph(term, max_steps, self.down_rels)

    def _traverseGraph(self, node, max_steps, rels):
        """_traverseGraph traverses the graph in a breath first manner
        and reports all the nodes reachable within max_steps."""
        remain = set([node])
        found = set()
        while len(remain) > 0 and max_steps != 0:
            novel = set()
            for t in remain:
                for rel in rels:
                    try:
                        novel.update(t.__dict__[rel])
                    except KeyError:
                        pass
            found.update(remain)
            remain = novel.difference(found)
            max_steps -= 1
        return found


class FreqAwareGeneOntology(GeneOntology):
    """GO hierarchy represents the Gene Ontology vocabulary. 
    
    It gets loaded from the xml file and, in conjunction with 
    an annotation file (GOA) the relative frequencies per term get
    estimated. this estimation respects the hierarchy of the 
    vocabulary.
    Further, this class provides methods to traverse the hierarchy
    in an easy way."""

    def __init__(self, fp, rels=UP_RELS):
        super(FreqAwareGeneOntology, self).__init__(fp, rels=rels)
        self.reset_freqs()

    def reset_freqs(self):
        self.cnts = dict()
        self.tot_cnts = [0] * NUM_ONT

    def estimate_freqs(self, annotations):
        for anno in annotations:
            try:
                self._update_counts(self.term_by_id(anno['TermNr']))
            except ValueError:
                logging.info("invalid annotation term_id in freq estim:" +
                             str(anno.term_id))

    def _update_counts(self, term):
        for cur_term in self.get_superterms_incl_queryterm(term):
            self.cnts[cur_term.id] = self.cnts.get(cur_term.id, 0) + 1
        self.tot_cnts[cur_term.aspect] += 1

    def get_term_frequency(self, term):
        term = self.ensure_term(term)
        try:
            freq = self.cnts.get(term.id, 0) / self.tot_cnts[term.aspect]
            return freq
        except ZeroDivisionError:
            return 0

    def last_common_ancestor(self, *terms):
        cand = self.get_superterms_incl_queryterm(terms[0])
        for t in terms[1:]:
            cand.intersection_update(self.get_superterms_incl_queryterm(t))
        lca = min(cand, key=self.get_term_frequency)
        return lca

    def lin_similarity(self, term1, term2):
        term1 = self.ensure_term(term1)
        term2 = self.ensure_term(term2)
        if term1.aspect != term2.aspect:
            # early abort, since the two terms will be by 
            # definition not similar
            sim = 0
        else:
            lca = self.last_common_ancestor(term1, term2)
            sim = (2 * math.log(self.get_term_frequency(lca)) /
                   (math.log(self.get_term_frequency(term1)) +
                    math.log(self.get_term_frequency(term2))))
        return sim


class AnnotationFilter(object):
    EXP_CODES = frozenset(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'])
    TRUST_IEA_REFS = frozenset([
        'GO_REF:0000002', 'GOA:interpro', 'GOA:interpro|GO_REF:0000002',  # InterPro
        'GO_REF:0000003', 'GOA:spec', 'GOA:spec|GO_REF:0000003'  # EC number
        'GO_REF:0000004', 'GOA:spkw', 'GOA:spkw|GO_REF:0000004',
        'GO_REF:0000037', 'GO_REF:0000038'  # SwissProt Keywords
        'GO_REF:0000023', 'GOA:spsl', 'GOA:spsl|GO_REF:0000023',
        'GO_REF:0000039', 'GO_REF:0000040',  # UniProtKB Subcellular Location
    ])

    @staticmethod
    def is_negated(a):
        return a.qualifier.find('NOT') >= 0

    @classmethod
    def is_exp_annotation(cls, a):
        return a.evidence in cls.EXP_CODES and not cls.is_negated(a)

    @classmethod
    def is_trusted_electronic(cls, a):
        return a.evidence == 'IEA' and a.db_ref in cls.TRUST_IEA_REFS and not cls.is_negated(a)

    @classmethod
    def is_exp_or_trusted_electronic(cls, a):
        return cls.is_exp_annotation(a) or cls.is_trusted_electronic(a)


# definition of the GOA Annotations. for performance reasons, we 
# keep this as a namedtuple collection. 
GOA_Annotation = collections.namedtuple('GOA_Annotation',
         ['db', 'db_obj_id', 'db_obj_sym', 'qualifier', 'term_id',
          'db_ref', 'evidence', 'with_from', 'aspect', 'db_obj_name',
          'db_obj_syn', 'db_obj_typ', 'taxon', 'date', 'assigned_by',
          'ext', 'gene_product_from_id'])


class AnnotationParser(object):
    def __init__(self, fp, factory=GOA_Annotation._make):
        self._needs_close = False
        if isinstance(fp, str):
            if fp.endswith('.gz'):
                from gzip import GzipFile
                fp = GzipFile(fp, 'rb')
                self._needs_close = True
            else:
                fp = open(fp, 'rb')
        self.fp = fp
        self.factory = factory

        self._read_headers()

    def _read_headers(self):
        pass

    def annotations(self):
        """Iterates over the annotations in the file yielding objects
        constructed by the factory argument passed to the constructor
        of this class for each annotation."""
        csv_reader = csv.reader((l for l in self.fp if not l.startswith('!')),
                                delimiter='\t')
        for row in csv_reader:
            yield self.factory(row)
        if self._needs_close:
            self.fp.close()

    def __iter__(self):
        return self.annotations()


