from __future__ import division

import collections

import time


def format_sciname(sci, short=False):
    p = set([sci.find(x) for x in ['(', 'serogroup', 'serotype', 'serovar',
                                   'biotype', 'subsp', 'pv.', 'bv.']])
    if sci.startswith('Escherichia coli'):
        p.add(sci.find('O'))
    p.discard(-1)
    p = min(p) if len(p) > 0 else len(sci)
    return {'species': sci[0:p], 'strain': sci[p:]}


class LazyProperty(object):
    """Decorator to evaluate a property only on access.

    Compute the attribute value and caches it in the instance.
    Python Cookbook (Denis Otkidach) http://stackoverflow.com/users/168352/denis-otkidach
    This decorator allows you to create a property which can be computed once and
    accessed many times."""

    def __init__(self, method, name=None):
        # record the unbound-method and the name
        self.method = method
        self.name = name or method.__name__
        self.__doc__ = method.__doc__

    def __get__(self, inst, cls):
        if inst is None:
            return self
        # compute, cache and return the instance's attribute value
        result = self.method(inst)
        # setattr redefines the instance's attribute so this doesn't get called again
        setattr(inst, self.name, result)
        return result


class KeyWrapper(object):
    '''
        Enables the use of functions, e.g. bisect, with a key function.
    '''
    def __init__(self, it, key):
        self.it = it
        self.key = key

    def __getitem__(self, i):
        return self.key(self.it[i])

    def __len__(self):
        return len(self.it)


class Singleton(type):
    """A meta-class to enforce a Singleton, e.g. a class that can be
    instantiated only exactly once.

    Modified from Python Cookbook, 3rd Edition, p 357ff.

    :Example:

        class Foo(metaclass=Singleton):
            def __init__(self):
                pass  #This part is executed only once
    """
    def __init__(self, *args, **kwargs):
        self.__instance = None
        super(Singleton, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        if self.__instance is None:
            self.__instance = super(Singleton, self).__call__(*args, **kwargs)
        return self.__instance


class ProteinEntry(object):
    """Model for a protein object

    This class provides an easy to use interface for a given protein
    form the database.

    If instantiated with an entry_nr only, no data is loaded until a
    property or method is accessed. Properties that need to access
    additional data or loaded lazily and are cached in the object
    (but not kept after deletion of object)."""
    def __init__(self, db, e):
        self._stored_entry = e
        self._db = db

    @LazyProperty
    def _entry(self):
        return (self._db.entry_by_entry_nr(self._stored_entry)
                if isinstance(self._stored_entry, int)
                else self._stored_entry)

    @classmethod
    def from_entry_nr(cls, db, eNr):
        # e = db.entry_by_entry_nr(eNr)
        return cls(db, int(eNr))

    @property
    def entry_nr(self):
        return int(self._entry['EntryNr'])

    @property
    def locus_start(self):
        return int(self._entry['LocusStart'])

    @property
    def locus_end(self):
        return int(self._entry['LocusEnd'])

    @property
    def strand(self):
        return int(self._entry['LocusStrand'])

    @LazyProperty
    def exons(self):
        return ExonStructure.from_entry_nr(self._db, self.entry_nr)

    @property
    def oma_group(self):
        return int(self._entry['OmaGroup'])

    @property
    def oma_hog(self):
        return self._entry['OmaHOG'].decode()

    @property
    def chromosome(self):
        return self._entry['Chromosome'].decode()

    @property
    def canonicalid(self):
        return self._entry['CanonicalId'].decode()

    @property
    def sequence_md5(self):
        return self._entry['MD5ProteinHash'].decode()

    @LazyProperty
    def genome(self):
        g = self._db.id_mapper['OMA'].genome_of_entry_nr(self._entry['EntryNr'])
        return Genome(self._db, g)

    @LazyProperty
    def omaid(self):
        return self._db.id_mapper['OMA'].map_entry_nr(self._entry['EntryNr'])

    @LazyProperty
    def cdna(self):
        return self._db.get_cdna(self._entry).decode()

    @property
    def gc_content(self):
        cdna = self.cdna
        cnts = list(map(cdna.count, 'GCAT'))
        try:
            return sum(cnts[0:2])/sum(cnts)
        except ZeroDivisionError:
            return 0

    @LazyProperty
    def sequence(self):
        return self._db.get_sequence(self._entry).decode()

    @property
    def sequence_length(self):
        return int(self._entry['SeqBufferLength']) - 1

    @LazyProperty
    def description(self):
        return self._db.get_description(self._entry).decode()

    @property
    def subgenome(self):
        return self._entry['SubGenome'].decode()

    @LazyProperty
    def hog_family_nr(self):
        from .db import Singleton as HOGSingleton
        try:
            fam = self._db.hog_family(self._entry)
        except HOGSingleton:
            fam = 0
        return fam

    @property
    def is_main_isoform(self):
        return (self._entry['AltSpliceVariant'] == 0 or
                self._entry['AltSpliceVariant'] == self._entry['EntryNr'])

    @LazyProperty
    def alternative_isoforms(self):
        return [ProteinEntry(self._db, e)
                for e in self._db.get_splicing_variants(self._entry)
                if e['EntryNr'] != self.entry_nr]

    def __repr__(self):
        return "<{}({}, {})>".format(self.__class__.__name__, self.entry_nr, self.omaid)

    def __len__(self):
        return self.sequence_length


class Genome(object):
    def __init__(self, db, g):
        self._genome = g
        self._db = db

    @property
    def ncbi_taxon_id(self):
        return int(self._genome['NCBITaxonId'])

    @property
    def uniprot_species_code(self):
        return self._genome['UniProtSpeciesCode'].decode()

    @property
    def sciname(self):
        return self._genome['SciName'].decode()

    @property
    def common_name(self):
        try:
            return self._genome['CommonName'].decode()
        except ValueError:
            return ""

    @property
    def synonym_name(self):
        return self._genome['SynName'].decode()

    @LazyProperty
    def species_and_strain_as_dict(self):
        return format_sciname(self.sciname)

    def species(self):
        return self.species_and_strain_as_dict['species']

    def strain(self):
        return self.species_and_strain_as_dict['strain']

    @property
    def url(self):
        return self._genome['Url'].decode()

    @property
    def source(self):
        return self._genome['Source'].decode()

    @property
    def release(self):
        return self._genome['Release'].decode()

    @property
    def last_modfied_timestamp(self):
        return self._genome['Date']

    @property
    def last_modified(self):
        return self.modification_date("%Y-%b-%d")

    def modification_date(self, fmt):
        if self._db.db_schema_version >= (3, 2):
            return time.strftime(fmt, time.localtime(self.last_modfied_timestamp))
        else:
            return 'n/a'

    @property
    def nr_entries(self):
        return int(self._genome['TotEntries'])

    @property
    def entry_nr_offset(self):
        return int(self._genome['EntryOff'])

    @LazyProperty
    def kingdom(self):
        # TODO: store directly in db
        return self._db.tax.get_parent_taxa(self._genome['NCBITaxonId'])[-1]['Name'].decode()

    @property
    def is_polyploid(self):
        return self._genome['IsPolyploid']

    @LazyProperty
    def lineage(self):
        return [lev['Name'].decode() for lev in self._db.tax.get_parent_taxa(
            self._genome['NCBITaxonId'])]

    @LazyProperty
    def chromosomes(self):
        chrs = collections.defaultdict(list)
        entry_tab = self._db.get_hdf5_handle().get_node('/Protein/Entries')
        for row in entry_tab.where('(EntryNr > {}) & (EntryNr <= {})'
                .format(self.entry_nr_offset, self.entry_nr_offset+self.nr_entries)):
            chrs[row['Chromosome'].decode()].append(row['EntryNr'])
        return chrs

    def __repr__(self):
        return "<{}({}, {})>".format(self.__class__.__name__, self.uniprot_species_code,
                                     self.ncbi_taxon_id)

    def __len__(self):
        return self.nr_entries


class PairwiseRelation(object):
    def __init__(self, db, relation):
        self._relation = relation
        self._db = db

    @property
    def distance(self):
        return float(self._relation['Distance'])

    @property
    def score(self):
        return float(self._relation['Score'])

    @property
    def alignment_overlap(self):
        return float(self._relation['AlignmentOverlap'])

    @property
    def synteny_conservation_local(self):
        return float(self._relation['SyntenyConservationLocal'])

    @property
    def confidence(self):
        return float(self._relation['Confidence'])

    @LazyProperty
    def rel_type(self):
        if not isinstance(self._relation['RelType'], str):
            type_map = self._db._get_pw_tab(self._relation['EntryNr1'], 'VPairs').get_enum("RelType")
            return type_map(self._relation['RelType'])
        else:
            return self._relation['RelType']

    @LazyProperty
    def entry_1(self):
        return ProteinEntry(self._db, self._db.entry_by_entry_nr(self._relation['EntryNr1']))

    @LazyProperty
    def entry_2(self):
        return ProteinEntry(self._db, self._db.entry_by_entry_nr(self._relation['EntryNr2']))


class GeneOntologyAnnotation(object):
    def __init__(self, db, anno):
        self.db = db
        self.anno = anno

    @LazyProperty
    def term(self):
        return self.db.gene_ontology.term_by_id(self.anno['TermNr'])

    @property
    def evidence(self):
        return self.anno['Evidence'].decode()

    @property
    def reference(self):
        return self.anno['Reference'].decode()

    @property
    def entry_nr(self):
        return int(self.anno['EntryNr'])

    @LazyProperty
    def aspect(self):
        from .geneontology import GOAspect
        return GOAspect.to_string(self.term.aspect)


class ExonStructure(object):
    def __init__(self, db, exons):
        self._stored = exons
        self._db = db

    @LazyProperty
    def _exons(self):
        return (self._db.get_exons(self._stored)
                if isinstance(self._stored, int)
                else self._stored)

    @classmethod
    def from_entry_nr(cls, db, eNr):
        return cls(db, int(eNr))

    def _iter_exons(self):
        if self._exons['Strand'][0] < 0:
            self._exons[::-1].sort(order='Start')
        else:
            self._exons.sort(order='Start')
        for exon in self._exons:
            yield Exon(exon)

    def __len__(self):
        return len(self._exons)

    def __repr__(self):
        return "<{}(entry_nr={}, nr_exons={})>"\
            .format(self.__class__.__name__,
                    self._exons[0]['EntryNr'], len(self))

    def __str__(self):
        exs = list(str(e) for e in self._iter_exons())
        if len(exs) > 1:
            return "join({})".format(", ".join(exs))
        else:
            return exs[0]


class Exon(object):
    def __init__(self, exon):
        self.exon = exon

    def __str__(self):
        if self.exon['Strand'] < 0:
            template = "complement({}..{})"
        else:
            template = "{}..{}"
        return template.format(self.exon['Start'], self.exon['End'])
