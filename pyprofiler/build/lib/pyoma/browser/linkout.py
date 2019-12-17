import collections
import operator
import os
import logging
import math
import re
import ftplib

from tqdm import tqdm
from lxml import etree
from .db import Database

logger = logging.getLogger(__name__)

"""Module to generate external data and crosslinks for NCBI.

"""


class NCBILinkOutXML(object):
    root_node = "not_set"
    provider_id = "9822"

    def __init__(self):
        root = etree.Element(self.root_node)
        for key, value in self.root_children().items():
            root.append(self.text_elemement(key, value))
        self.tree = etree.ElementTree(root)
        self._add_doctype()

    def root_children(self):
        return {}

    def _add_doctype(self):
        self.tree.docinfo.public_id = '-//NLM//DTD LinkOut 1.0//EN'
        self.tree.docinfo.system_url = 'https://www.ncbi.nlm.nih.gov/projects/linkout/doc/LinkOut.dtd'

    def text_elemement(self, tag, text):
        el = etree.Element(tag)
        el.text = text
        return el

    def write(self, fh):
        fh.write(etree.tostring(self.tree, pretty_print=True, xml_declaration=True, encoding='utf-8'))


class Provider(NCBILinkOutXML):
    root_node = "Provider"

    def root_children(self):
        elements = collections.OrderedDict(
            [("ProviderId", self.provider_id),
             ("Name", "OMA Browser: Orthologous MAtrix"),
             ("NameAbbr", "OMA"),
             #("SubjectType", "taxonomy/phylogenetic"),
             ("Url", "http://omabrowser.org"),
             ("Brief", "OMA is a method and database for the inference of orthologs among complete genomes. "
                       "We provide browsable orthology predictions, APIs, flat file downloads among thousands "
                       "of genomes.")])
        return elements


class Resource(NCBILinkOutXML):
    root_node = "LinkSet"
    link_id = 1

    def _add_objs(self, accs):
        objsel = etree.Element("ObjectSelector")
        objsel.append(self.text_elemement("Database", self.database()))
        objlst = etree.Element("ObjectList")
        objsel.append(objlst)
        for acc in accs:
            objlst.append(self.object_node(acc))
        return objsel

    def object_node(self, acc):
        return self.text_elemement("ObjId", acc)

    def _add_url_section(self, acc):
        el = etree.Element('ObjectUrl')
        el.append(self.text_elemement('Base', self.base_url()))
        nxt = rule = etree.Element("Rule")
        for k, rule_part in enumerate(self.rule_url(acc)):
            if isinstance(rule_part, str):
                if k == 0:
                    nxt.text = rule_part
                else:
                    nxt.tail = rule_part
            elif rule_part.tag == etree.Entity:
                nxt.append(rule_part)
                nxt = rule_part
        el.append(rule)
        el.append(self.text_elemement('SubjectType', "taxonomy/phylogenetic"))
        return el

    def add_link(self, accs):
        lnk = etree.Element("Link")
        lnk.append(self.text_elemement('LinkId', str(self.link_id)))
        lnk.append(self.text_elemement('ProviderId', self.provider_id))
        lnk.append(self._add_objs(accs))
        lnk.append(self._add_url_section(accs))
        self.tree.getroot().append(lnk)
        self._bump_link_id()

    @classmethod
    def _bump_link_id(cls):
        cls.link_id += 1

    def database(self):
        return "not set"

    def base_url(self):
        return "https://omabrowser.org/oma/hogs/"

    def rule_url(self, acc):
        return "",


class GenesResource(Resource):
    DISKSIZE_HEADER = 200
    DISKSIZE_PER_LINK = 435
    base_name = 'resource_genes'

    def base_url(self):
        return "https://omabrowser.org/cgi-bin/gateway.pl/"

    def rule_url(self, acc):
        return "?f=DisplayEntry&p1=" + next(iter(acc.values())),

    def database(self):
        return "Gene"


class ProteinResource(Resource):
    DISKSIZE_HEADER = 500
    DISKSIZE_PER_LINK = 45
    base_name = 'resource_protein'

    def base_url(self):
        return "https://omabrowser.org/oma/hogs/"

    def object_node(self, acc):
        return self.text_elemement("Query", "{}[accn]".format(acc))

    def rule_url(self, acc):
        return etree.Entity("lo.pacc"), "/vis/"

    def database(self):
        return "Protein"


class TaxonomyResource(Resource):
    DISKSIZE_HEADER = 200
    DISKSIZE_PER_LINK = 435
    base_name = 'resource_taxonomy'

    def database(self):
        return "taxonomy"

    def base_url(self):
        return "https://omabrowser.org/cgi-bin/gateway.pl/"

    def rule_url(self, acc):
        return "?f=DisplayOS&p1=" + next(iter(acc.values())),


class LinkoutBuffer(object):
    def __init__(self, resource, outdir='/tmp', bulk_add=True, max_file_size=20*2**20):
        self.max_records = math.floor((max_file_size - resource.DISKSIZE_HEADER) /
                                      resource.DISKSIZE_PER_LINK)
        self.cur_nr = 0
        self.buf = []
        self.bulk_add = bulk_add
        self.resource_type = resource
        self.outdir = outdir
        logger.info('Setup Linkout buffer for {} with max {} records ({}bytes) per file, bulk_add={}'
                    .format(resource.__name__, self.max_records, max_file_size, bulk_add))

    def add(self, obj):
        self.buf.append(obj)
        if len(self.buf) >= self.max_records:
            self.flush()

    def flush(self):
        res = self.resource_type()
        if self.bulk_add:
            res.add_link(self.buf)
        else:
            for obj in self.buf:
                res.add_link(obj)
        fn = os.path.join(self.outdir,
                          '{}_{:02d}.xml'.format(res.base_name, self.cur_nr))
        with open(fn, 'wb') as fh:
            res.write(fh)
        self.cur_nr += 1
        self.buf = []


class GenesPriorizationHandler(object):
    """Adapter to LinkoutBuffer to select only a limited number of crossrefs.
    NCBI linkout caps at 10%"""

    def __init__(self, max_linkouts=None, db=None, **kwargs):
        self.max_links = int(max_linkouts) if max_linkouts else 20357436//10  # obtained in Jan2018
        logger.info('Limiting Genes to {} links max'.format(self.max_links))
        self.genes_buffer = LinkoutBuffer(GenesResource, **kwargs)
        self.genes = []
        self.db = db

    def add(self, key, value):
        self.genes.append((key, value))

    def _genome_size_map(self):
        gs = self.db.get_hdf5_handle().get_node('/Genome').read()
        return {row['UniProtSpeciesCode'].decode(): row['TotEntries'] for row in gs}

    def flush(self):
        priority_prefixes = ['HUMAN', 'MOUSE', 'RATNO', 'PIGXX', 'DRO', 'SCH', 'YEAST', 'ARA',
                             'WHEAT', 'PLAF', 'ECO', 'BAC', 'PANTR', 'ORY', 'GOSHI', 'BRA',
                             'DANRE', 'CAE', 'MYC', 'STR', 'MAIZE', 'GORGO', 'PANTR', 'PONAB',
                             'MACMU', 'YARLI', 'PEDHC', 'TRICA', 'XENTR', 'YERPE', 'POPTR']
        pat = re.compile(r"^({})".format('|'.join(priority_prefixes)))
        if len(self.genes) > self.max_links:
            # final sort order will be 'priority genome', genome size and proteins within genome
            self.genes.sort(key=operator.itemgetter(1))
            if self.db is not None:
                genome_size = self._genome_size_map()
                self.genes.sort(key=lambda x: genome_size[x[1][0:5]], reverse=True)
            self.genes.sort(key=lambda x: pat.match(x[1]) is None)
        for link_acc, link_target in self.genes[0:self.max_links]:
            self.genes_buffer.add({link_acc: link_target})
        self.genes_buffer.flush()

        c = collections.defaultdict(int)
        for acc, target in self.genes[self.max_links:]:
            c[target[0:5]] += 1
        logger.info('Skipping genes link in the following species: {}'.format(c))


def prepare_linkout_files(outdir='/tmp', infile='../pyomabrowser/OmaServer.h5'):
    prov = Provider()
    with open(os.path.join(outdir, 'provider.xml'), 'wb') as fh:
        prov.write(fh)

    db = Database(infile)
    xrefs = db.get_hdf5_handle().get_node('/XRef')
    xref_source_enum = xrefs.get_enum('XRefSource')

    protein_buffer = LinkoutBuffer(ProteinResource, outdir=outdir, bulk_add=True)
    genes_buffer = GenesPriorizationHandler(db=db, outdir=outdir, bulk_add=False)
    for xref in tqdm(xrefs):
        if xref['XRefSource'] == xref_source_enum['RefSeq']:
            protein_buffer.add(xref['XRefId'].decode())
        elif xref['XRefSource'] == xref_source_enum['EntrezGene']:
            genes_buffer.add(xref['XRefId'].decode(),
                             db.id_mapper['OMA'].map_entry_nr(xref['EntryNr']))
    protein_buffer.flush()
    genes_buffer.flush()

    with open(os.path.join(outdir, 'resource_taxonomy.xml'), 'wb') as fh:
        taxs = TaxonomyResource()
        for row in db.id_mapper['OMA'].genome_table:
            taxs.add_link({str(row['NCBITaxonId']): row['UniProtSpeciesCode'].decode()})
        taxs.write(fh)


def copy_to_ncbi(dir, password, host='ftp-private.ncbi.nlm.nih.gov', user='omabrow'):
    with ftplib.FTP(host, user, password) as session:
        session.cwd('/holdings')

        for fname in os.listdir(dir):
            if fname.endswith('.xml'):
                with open(os.path.join(dir, fname), 'rb') as fh:
                    cmd = "STOR {}".format(fname)
                    session.storbinary(cmd, fp=fh)
                logger.info('finished transfering '+fname)