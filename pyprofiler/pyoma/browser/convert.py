from __future__ import division, print_function
from builtins import str, chr, range, object, super, bytes

import pandas
from future.standard_library import hooks
from PySAIS import sais
from tempfile import NamedTemporaryFile
from tqdm import tqdm
import csv
import resource
import tables
import numpy
import numpy.lib.recfunctions
import os
import subprocess
import errno
import json
import time
import familyanalyzer
import re
import multiprocessing as mp
import lxml.html
import collections
import gzip
import hashlib
import itertools
import operator
import fileinput

from .. import common
from . import locus_parser
from . import tablefmt
from .KmerEncoder import KmerEncoder
from .OrthoXMLSplitter import OrthoXMLSplitter
from .geneontology import GeneOntology, OntologyParser
from .synteny import SyntenyScorer
from .homoeologs import HomeologsConfidenceCalculator

with hooks():
    import urllib.request


class DarwinException(Exception):
    pass


def callDarwinExport(func, drwfile=None):
    """Function starts a darwin session, loads convert.drw file
    and calls the darwin function passed as argument. The output
    is expected to be written by darwin in json format into the
    file specified by 'outfn'.
    This function returns the parsed json datastructure"""

    with NamedTemporaryFile(suffix='.dat') as tmpfile:
        if drwfile is None:
            drwfile = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")
        # with open(os.devnull, 'w') as DEVNULL:
        stacksize = resource.getrlimit(resource.RLIMIT_STACK)
        common.package_logger.info('current stacklimit: {}'.format(stacksize))
        common.package_logger.info('setting stacklimit: {}'.format((max(stacksize)-1, stacksize[1])))
        resource.setrlimit(resource.RLIMIT_STACK, (min(stacksize), stacksize[1]))
        p = subprocess.Popen(['darwin', '-q', '-E', '-B'], stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        drw_cmd = "outfn := '{}': ReadProgram('{}'): {}; done;".format(
            tmpfile.name,
            drwfile,
            func).encode('utf-8')
        common.package_logger.debug('calling darwin function: {}'.format(func))
        (stdout, stderr) = p.communicate(input=drw_cmd)
        if p.returncode > 0:
            raise DarwinException(p.stderr.read())

        trans_tab = "".join(str(chr(x)) for x in range(128)) + " " * 128
        with open(tmpfile.name, 'r') as jsonData:
            rawdata = jsonData.read()
            return json.loads(rawdata.translate(trans_tab))


def uniq(seq):
    """return uniq elements of a list, preserving order

    :param seq: an iterable to be analyzed
    """
    seen = set()
    return [x for x in seq if not (x in seen or seen.add(x))]


def silentremove(filename):
    """Function to remove a given file. No exception is raised if the
    file does not exist. Other errors are passed to the user.
    :param filename: the path of the file to be removed"""
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


def gz_is_empty(fname):
    """Test if gzip file fname is empty

    Return True if the uncompressed data in fname has zero length
    or if fname itself has zero length
    Raises OSError if fname has non-zero length and is not a gzip file
    """
    with gzip.open(fname, 'rb') as f:
        data = f.read(1)
    return len(data) == 0


def load_tsv_to_numpy(args):
    fn, off1, off2, swap = args
    rel_desc = tablefmt.PairwiseRelationTable
    # we need to get the enum as a dict to be able to extend it
    # with the reversed labels, i.e. n:1
    relEnum = rel_desc.columns['RelType'].enum._names
    relEnum['n:1'] = relEnum['m:1']
    relEnum['1:m'] = relEnum['1:n']
    relEnum['n:m'] = relEnum['m:n']
    read_dir = -1 if swap else 1
    tsv_dtype = [('EntryNr1', 'i4'), ('EntryNr2', 'i4'), ('Score', 'f4'), ('RelType', 'i1'),
                 ('AlignmentOverlap', 'f2'), ('Distance', 'f4')]
    for curNr, curFn in enumerate([fn, fn.replace('.ext.', '.')]):
        try:
            if gz_is_empty(curFn):
                return numpy.empty(0, dtype=tables.dtype_from_descr(rel_desc))
            with gzip.GzipFile(curFn) as fh:
                data = numpy.genfromtxt(fh, dtype=tsv_dtype,
                                        names=[_[0] for _ in tsv_dtype],
                                        delimiter='\t',
                                        usecols=(0, 1, 2, 3, 4, 5),
                                        converters={'EntryNr1': lambda nr: int(nr) + off1,
                                                    'EntryNr2': lambda nr: int(nr) + off2,
                                                    'RelType': lambda rel: (relEnum[rel[::read_dir].decode()]
                                                                            if len(rel) <= 3
                                                                            else relEnum[rel.decode()]),
                                                    'Score': lambda score: float(score) / 100})
                break
        except OSError as e:
            if curNr < 1:
                common.package_logger.info('tried to load {}'.format(curFn))
                pass
            else:
                raise e

    if swap:
        reversed_cols = tuple(data.dtype.names[z] for z in (1, 0, 2, 3, 4, 5))
        data.dtype.names = reversed_cols
    full_table = numpy.empty(data.size, dtype=tables.dtype_from_descr(rel_desc))
    common_cols = list(data.dtype.names)
    full_table[common_cols] = data[common_cols]
    for col_not_in_tsv in set(full_table.dtype.names) - set(data.dtype.names):
        full_table[col_not_in_tsv] = rel_desc.columns[col_not_in_tsv].dflt
    return full_table


def read_vps_from_tsv(gs, ref_genome):
    ref_genome_idx = gs.get_where_list('(UniProtSpeciesCode=={!r})'.
                                       format(ref_genome))[0]
    job_args = []
    for g in range(len(gs)):
        if g == ref_genome_idx:
            continue
        g1, g2 = sorted((g, ref_genome_idx,))
        off1, off2 = gs.read_coordinates(numpy.array((g1, g2)), 'EntryOff')
        fn = os.path.join(os.environ['DARWIN_OMADATA_PATH'], 'Phase4',
                          gs.cols.UniProtSpeciesCode[g1].decode(),
                          gs.cols.UniProtSpeciesCode[g2].decode() + ".orth.txt.gz")
        tup = (fn, off1, off2, g1 != ref_genome_idx)
        common.package_logger.info('adding job: {}'.format(tup))
        job_args.append(tup)

    pool = mp.Pool(processes=min(os.cpu_count(), 10))
    all_pairs = pool.map(load_tsv_to_numpy, job_args)
    pool.close()
    return numpy.lib.recfunctions.stack_arrays(all_pairs, usemask=False)


class DataImportError(Exception):
    pass


def _load_taxonomy_without_ref_to_itselfs(data):
    dtype = tables.dtype_from_descr(tablefmt.TaxonomyTable)
    arr = numpy.array([tuple(x) for x in data], dtype=dtype)
    clean = arr[numpy.where(arr['NCBITaxonId'] != arr['ParentTaxonId'])]
    return clean


def compute_ortholog_types(data, genome_offs):
    """this function computes the type of orthologs from the data and sets in
    the RelType column.

    :param data: a numpy recarray corresponding to the `numpy.dtype` of
           `tablefmt.PairwiseRelationTable`
    :param genome_offs: a numpy array with the genome offsets, i.e. the entry
           numbers where the next genome starts

    :returns: a modified version of data
    """
    typEnum = tablefmt.PairwiseRelationTable.columns.get('RelType').enum
    query_type = {val: 'm' if cnt > 1 else '1'
                  for val, cnt in zip(*numpy.unique(data['EntryNr2'],
                                                    return_counts=True))}

    def genome_idx(enr):
        return numpy.searchsorted(genome_offs, enr - 1, side='right')

    g0 = genome_idx(data[0]['EntryNr2'])
    it = numpy.nditer(data, flags=['c_index'], op_flags=['readwrite'])
    while not it.finished:
        row0 = it[0]
        i1 = it.index + 1
        # we move i1 forward to the row where the next genome starts, i.e. the
        # current query changes the species or the query itself changes
        while i1 < len(data):
            row1 = data[i1]
            g1 = genome_idx(row1['EntryNr2'])
            if g1 != g0 or row0['EntryNr1'] != row1['EntryNr1']:
                break
            i1 += 1
        subj_type = 'n' if i1 - it.index > 1 else '1'
        while not it.finished and it.index < i1:
            typ = '{}:{}'.format(query_type[int(it[0]['EntryNr2'])], subj_type)
            it[0]['RelType'] = typEnum[typ]
            it.iternext()
        g0 = g1


def get_or_create_tables_node(h5, path, desc=None):
    """return the node of a given path from the h5 file

    If the node does not yet exist, it is created (including potential
    inexistant internal nodes).

    :param h5: Handle to the hdf5 object
    :param str path: Path of the node to return
    :param str desc: Description to be added to the node"""
    try:
        grp = h5.get_node(path)
    except tables.NoSuchNodeError:
        base, name = os.path.split(path)
        grp = h5.create_group(base, name, title=desc, createparents=True)
    return grp


class DarwinExporter(object):
    DB_SCHEMA_VERSION = '3.2'
    DRW_CONVERT_FILE = os.path.abspath(os.path.splitext(__file__)[0] + '.drw')

    def __init__(self, path, logger=None, mode=None):
        self.logger = logger if logger is not None else common.package_logger
        fn = os.path.normpath(os.path.join(
            os.getenv('DARWIN_BROWSERDATA_PATH', ''),
            path))
        if mode is None:
            mode = 'append' if os.path.exists(fn) else 'write'
        self._compr = tables.Filters(complevel=6, complib='zlib', fletcher32=True)
        self.h5 = tables.open_file(fn, mode=mode[0], filters=self._compr)
        self.logger.info("opened {} in {} mode, options {}".format(
            fn, mode, str(self._compr)))
        if mode == 'write':
            self.h5.root._f_setattr('convertion_start', time.strftime("%c"))

    def call_darwin_export(self, func):
        return callDarwinExport(func, self.DRW_CONVERT_FILE)

    def _get_or_create_node(self, path, desc=None):
        return get_or_create_tables_node(self.h5, path, desc)

    def create_table_if_needed(self, parent, name, drop_data=False, **kwargs):
        """create a table if needed.

        The function only checks whether a table exists with that name,
        but not if it is compatible with the passed arguments.
        if you pass data with the `obj` argument, this data is appended
        to the table. If you set `drop_data` to True, data that was
        previously in the existing table is dropped prior to adding new
        data."""
        try:
            tab = self.h5.get_node(parent, name=name)
            if drop_data:
                tab.remove_rows(0, tab.nrows)
            if 'obj' in kwargs:
                tab.append(kwargs['obj'])
        except tables.NoSuchNodeError:
            tab = self.h5.create_table(parent, name, **kwargs)
        return tab

    def get_version(self):
        """return version of the dataset.

        Default implementation searches for 'mname' in Matrix or matrix_stats.drw files.
        """
        for fname in ('Matrix', 'matrix_stats.drw'):
            with open(os.path.join(os.environ['DARWIN_BROWSERDATA_PATH'], fname), 'r') as fh:
                for i, line in enumerate(fh):
                    if line.startswith('mname :='):
                        match = re.match(r'mname := \'(?P<version>[^\']*)\'', line)
                        return match.group('version')
                    if i > 1000:
                        break
        raise DataImportError('No version information found')

    def add_version(self):
        version = self.get_version()
        self.h5.set_node_attr('/', 'oma_version', version)
        self.h5.set_node_attr('/', 'pytables', tables.get_pytables_version())
        self.h5.set_node_attr('/', 'hdf5_version', tables.get_hdf5_version())
        self.h5.set_node_attr('/', 'db_schema_version', self.DB_SCHEMA_VERSION)

    def add_species_data(self):
        cache_file = os.path.join(
            os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
            'pyoma', 'gs.json')
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as fd:
                data = json.load(fd)
        else:
            data = self.call_darwin_export('GetGenomeData();')
        gstab = self.h5.create_table('/', 'Genome', tablefmt.GenomeTable,
                                     expectedrows=len(data['GS']))
        gs_data = self._parse_date_columns(data['GS'], gstab)
        self._write_to_table(gstab, gs_data)
        gstab.cols.NCBITaxonId.create_csindex(filters=self._compr)
        gstab.cols.UniProtSpeciesCode.create_csindex(filters=self._compr)
        gstab.cols.EntryOff.create_csindex(filters=self._compr)

        taxtab = self.h5.create_table('/', 'Taxonomy', tablefmt.TaxonomyTable,
                                      expectedrows=len(data['Tax']))
        self._write_to_table(taxtab, _load_taxonomy_without_ref_to_itselfs(data['Tax']))
        taxtab.cols.NCBITaxonId.create_csindex(filters=self._compr)

    def _parse_date_columns(self, data, tab):
        """convert str values in a date column to epoch timestamps"""
        time_cols = [i for i, col in enumerate(tab.colnames) if tab.coldescrs[col].kind == 'time']
        dflts = [tab.coldflts[col] for col in tab.colnames]

        def map_data(col, data):
            try:
                val = data[col]
                if col in time_cols and isinstance(val, str):
                    for fmt in ('%b %d, %Y', '%B %d, %Y', '%d.%m.%Y', '%Y%m%d'):
                        try:
                            date = time.strptime(val, fmt)
                            return time.mktime(date)
                        except ValueError:
                            pass
                    raise ValueError("Cannot parse date of '{}'".format(val))
                return val
            except IndexError:
                return dflts[col]

        arr = numpy.empty(len(data), dtype=tab.dtype)
        for i, row in enumerate(data):
            as_tup = tuple(map_data(c, row) for c in range(len(dflts)))
            arr[i] = as_tup
        return arr

    def _convert_to_numpyarray(self, data, tab):
        """convert a list of list dataset into a numpy rec array that
        corresponds to the table definition of `tab`.

        :param data: the data to be converted.
        :param tab: a pytables table node."""

        enum_cols = {i: tab.get_enum(col) for (i, col) in enumerate(tab.colnames)
                     if tab.coltypes[col] == 'enum'}
        dflts = [tab.coldflts[col] for col in tab.colnames]

        def map_data(col, data):
            try:
                val = data[col]
                return enum_cols[col][val]
            except IndexError:
                return dflts[col]
            except KeyError:
                return val

        arr = numpy.empty(len(data), dtype=tab.dtype)
        for i, row in enumerate(data):
            as_tup = tuple(map_data(c, row) for c in range(len(dflts)))
            arr[i] = as_tup
        return arr

    def add_orthologs(self):
        genome_offs = self.h5.root.Genome.col('EntryOff')
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'VPairs' not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                    'pyoma', 'vps', '{}.json'.format(genome))
                if os.path.exists(cache_file):
                    with open(cache_file, 'r') as fd:
                        data = json.load(fd)
                elif ((not os.getenv('DARWIN_OMADATA_PATH') is None) and
                      os.path.exists(os.path.join(
                           os.environ['DARWIN_OMADATA_PATH'], 'Phase4'))):
                    # try to read from Phase4 in parallel.
                    data = read_vps_from_tsv(self.h5.root.Genome,
                                             genome.encode('utf-8'))
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export('GetVPsForGenome({})'.format(genome))

                vp_tab = self.h5.create_table(rel_node_for_genome, 'VPairs', tablefmt.PairwiseRelationTable,
                                              expectedrows=len(data))
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                if numpy.any(data['RelType'] >= tablefmt.PairwiseRelationTable.columns.get('RelType').enum['n/a']):
                    compute_ortholog_types(data, genome_offs)
                self._write_to_table(vp_tab, data)
                vp_tab.cols.EntryNr1.create_csindex()

    def add_same_species_relations(self):
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'within' not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                    'pyoma', 'cps', '{}.json'.format(genome))
                if os.path.exists(cache_file):
                    with open(cache_file, 'r') as fd:
                        data = json.load(fd)
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export('GetSameSpeciesRelations({})'.format(genome))

                ss_tab = self.h5.create_table(rel_node_for_genome, 'within', tablefmt.PairwiseRelationTable,
                                              expectedrows=len(data))
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, ss_tab)
                self._write_to_table(ss_tab, data)
                ss_tab.cols.EntryNr1.create_csindex()

    def add_synteny_scores(self):
        """add synteny scores of pairwise relations to database.

        Current implementation only computes synteny scores for
        homoeologs, but easy to extend. Question is rather if we
        need synteny scores for all genome pairs, and if not, how
        to select.

        The computations of the scores are done using :mod:`synteny`
        module of this package."""
        # TODO: compute for non-homoeologs relation as well.
        self.logger.info("Adding synteny scores for polyploid genomes")
        polyploid_genomes = self.h5.root.Genome.where('IsPolyploid==True')
        for genome in polyploid_genomes:
            genome_code = genome['UniProtSpeciesCode'].decode()
            self.logger.info('compute synteny score for {}'.format(genome_code))
            synteny_scorer = SyntenyScorer(self.h5, genome_code)
            rels = synteny_scorer.compute_scores()
            self._callback_store_rel_data(
                genome_code, rels, [('SyntenyConservationLocal', 'mean_synteny_score')])

    def add_homoeology_confidence(self):
        """adds the homoeology confidence scores to the database.

        This method should be called only after the synteny scores have
        been computed and added to the database.

        The computations are done using :mod:`homoeologs` module."""
        self.logger.info("Adding homoeolog confidence scores")
        polyploid_genomes = self.h5.root.Genome.where('IsPolyploid==True')
        for genome in polyploid_genomes:
            genome_code = genome['UniProtSpeciesCode'].decode()
            self.logger.info("compute homoeolog confidence for {}".format(genome_code))
            homoeolg_scorer = HomeologsConfidenceCalculator(self.h5, genome_code)
            rels = homoeolg_scorer.calculate_scores()
            self._callback_store_rel_data(
                genome_code, rels, [("Confidence", "fuzzy_confidence_scaled")])

    def _callback_store_rel_data(self, genome, rels_df, assignments):
        tab = self.h5.get_node('/PairwiseRelation/{}/within'.format(genome))
        df_all = pandas.DataFrame(tab.read())
        if 'entry_nr1' in list(rels_df):
            enr_col_names = ['entry_nr1', 'entry_nr2']
        else:
            enr_col_names = ['EntryNr1', 'EntryNr2']
        merged = pandas.merge(df_all, rels_df, how="left", left_on=['EntryNr1', 'EntryNr2'],
                              right_on=enr_col_names, validate='one_to_one')

        for target, source in assignments:
            # replace NaN in column from rels_df by the default value of the target column
            merged.loc[merged[source].isnull(), source] = tab.coldescrs[target].dflt
            # update the data in the target hdf5 column by the source column data
            tab.modify_column(column=merged[source].as_matrix(), colname=target)
        tab.flush()

    def _add_sequence(self, sequence, row, sequence_array, off, typ="Seq"):
        # add ' ' after each sequence (Ascii is smaller than
        # any AA, allows to build PAT array with split between
        # sequences.
        seqLen = len(sequence) + 1
        row[typ + 'BufferOffset'] = off
        row[typ + 'BufferLength'] = seqLen
        seqNumpyObj = numpy.ndarray((seqLen,),
                                    buffer=(sequence + " ").encode('utf-8'),
                                    dtype=tables.StringAtom(1))
        sequence_array.append(seqNumpyObj)
        if typ == "Seq":
            row['MD5ProteinHash'] = hashlib.md5(sequence.encode('utf-8')).hexdigest()
        return seqLen

    def add_proteins(self):
        gsNode = self.h5.get_node('/Genome')
        nrProt = sum(gsNode.cols.TotEntries)
        nrAA = sum(gsNode.cols.TotAA)
        protGrp = self._get_or_create_node('/Protein', "Root node for protein (oma entries) information")
        protTab = self.h5.create_table(protGrp, 'Entries', tablefmt.ProteinTable,
                                       expectedrows=nrProt)
        seqArr = self.h5.create_earray(protGrp, 'SequenceBuffer',
                                       tables.StringAtom(1), (0,), 'concatenated protein sequences',
                                       expectedrows=nrAA + nrProt)
        cdnaArr = self.h5.create_earray(protGrp, 'CDNABuffer',
                                        tables.StringAtom(1), (0,), 'concatenated cDNA sequences',
                                        expectedrows=3 * nrAA + nrProt)
        seqOff = cdnaOff = 0
        loc_parser = locus_parser.LocusParser()
        for gs in gsNode.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            cache_file = os.path.join(
                os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                'pyoma', 'prots', '{}.json'.format(genome))
            if os.path.exists(cache_file):
                with open(cache_file, 'r') as fd:
                    data = json.load(fd)
            else:
                data = self.call_darwin_export('GetProteinsForGenome({})'.format(genome))

            if len(data['seqs']) != gs['TotEntries']:
                raise DataImportError('number of entries ({:d}) does '
                                      'not match number of seqs ({:d}) for {}'
                                      .format(len(data['seqs']), gs['TotEntries'], genome))

            locTab = self.h5.create_table('/Protein/Locus',
                                          genome, tablefmt.LocusTable, createparents=True,
                                          expectedrows=gs['TotEntries'] * 4)

            for nr in range(gs['TotEntries']):
                eNr = data['off'] + nr + 1
                protTab.row['EntryNr'] = eNr
                protTab.row['OmaGroup'] = data['ogs'][nr]

                seqOff += self._add_sequence(data['seqs'][nr], protTab.row, seqArr, seqOff)
                cdnaOff += self._add_sequence(data['cdna'][nr], protTab.row, cdnaArr, cdnaOff, 'CDNA')

                protTab.row['Chromosome'] = data['chrs'][nr]
                protTab.row['AltSpliceVariant'] = data['alts'][nr]
                protTab.row['OmaHOG'] = b" "  # will be assigned later
                protTab.row['CanonicalId'] = b" "  # will be assigned later

                locus_str = data['locs'][nr]
                try:
                    locus_tab = loc_parser.parse(locus_str, eNr)
                    locTab.append(locus_tab)
                    len_cds = sum(z['End'] - z['Start']+1 for z in locus_tab)
                    if len_cds != protTab.row['CDNABufferLength']-1:
                        self.logger.warning("sum of exon lengths differ with cdna sequence for {}: {} vs {}"
                                            .format(eNr, len_cds, protTab.row['CDNABufferLength']-1))

                    protTab.row['LocusStart'] = locus_tab['Start'].min()
                    protTab.row['LocusEnd'] = locus_tab['End'].max()
                    protTab.row['LocusStrand'] = locus_tab[0]['Strand']
                except ValueError as e:
                    self.logger.warning(e)
                protTab.row['SubGenome'] = data['subgenome'][nr].encode('ascii')
                protTab.row.append()
            protTab.flush()
            seqArr.flush()
            for n in (protTab, seqArr, locTab):
                if n.size_in_memory != 0:
                    self.logger.info('worte %s: compression ratio %3f%%' %
                                     (n._v_pathname, 100 * n.size_on_disk / n.size_in_memory))
        protTab.cols.EntryNr.create_csindex(filters=self._compr)
        protTab.cols.MD5ProteinHash.create_csindex(filters=self._compr)

    def _write_to_table(self, tab, data):
        if len(data)>0:
            tab.append(data)
        self.logger.info('wrote %s : compression ratio %.3f%%' %
                         (tab._v_pathname, 100 * tab.size_on_disk / tab.size_in_memory))

    def add_hogs(self):
        hog_path = os.path.normpath(os.path.join(
            os.environ['DARWIN_NETWORK_SCRATCH_PATH'],
            'pyoma', 'split_hogs'))
        entryTab = self.h5.get_node('/Protein/Entries')
        tree_filename = os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            'speciestree.nwk')
        if not os.path.exists(hog_path):
            hog_file = os.path.join(os.environ['DARWIN_BROWSERDATA_PATH'],
                                    '..', 'downloads', 'oma-hogs.orthoXML.gz')
            splitter = OrthoXMLSplitter(hog_file, cache_dir=hog_path)
            splitter()
        hog_converter = HogConverter(entryTab)
        hog_converter.attach_newick_taxonomy(tree_filename)
        hogTab = self.h5.create_table('/', 'HogLevel', tablefmt.HOGsTable,
                                      'nesting structure for each HOG', expectedrows=1e8)
        self.orthoxml_buffer = self.h5.create_earray('/OrthoXML', 'Buffer',
                                                     tables.StringAtom(1), (0,), 'concatenated orthoxml files',
                                                     expectedrows=1e9, createparents=True)
        self.orthoxml_index = self.h5.create_table('/OrthoXML', 'Index', tablefmt.OrthoXmlHogTable,
                                                   'Range index per HOG into OrthoXML Buffer', expectedrows=5e6)
        for root, dirs, filenames in os.walk(hog_path):
            for fn in filenames:
                try:
                    levels = hog_converter.convert_file(os.path.join(root, fn))
                    hogTab.append(levels)
                    fam_nrs = set([z[0] for z in levels])
                    self.add_orthoxml(os.path.join(root, fn), fam_nrs)
                except Exception as e:
                    self.logger.error('an error occured while processing ' + fn + ':')
                    self.logger.exception(e)

        hog_converter.write_hogs()

    def add_orthoxml(self, orthoxml_path, fam_nrs):
        """append orthoxml file content to orthoxml_buffer array and add index for the HOG family"""
        if len(fam_nrs) > 1:
            self.logger.warning('expected only one family per HOG file, but found {}: {}'
                                .format(len(fam_nrs), fam_nrs))
            self.logger.warning(' --> the orthoxml files per family will be not correct, '
                                'i.e. they will contain all families of this file.')
        with open(orthoxml_path, 'r') as fh:
            orthoxml = fh.read().encode('utf-8')
            offset = len(self.orthoxml_buffer)
            length = len(orthoxml)
            self.orthoxml_buffer.append(numpy.ndarray((length,),
                                                      buffer=orthoxml, dtype=tables.StringAtom(1)))
            for fam in fam_nrs:
                row = self.orthoxml_index.row
                row['Fam'] = fam
                row['HogBufferOffset'] = offset
                row['HogBufferLength'] = length
                offset += length
                row.append()

    def xref_databases(self):
        return os.path.join(os.environ['DARWIN_BROWSERDATA_PATH'], 'ServerIndexed.db')

    def add_xrefs(self):
        self.logger.info('start extracting XRefs, EC and GO annotations')
        db_parser = DarwinDbEntryParser()
        xref_tab = self.h5.create_table('/', 'XRef', tablefmt.XRefTable,
                                        'Cross-references of proteins to external ids / descriptions',
                                        expectedrows=1e8)

        ec_tab = self.h5.create_table('/Annotations', 'EC', tablefmt.ECTable, 'Enzyme Commission annotations',
                                      expectedrows=1e7, createparents=True)
        gs = self.h5.get_node('/Genome').read()
        with DescriptionManager(self.h5, '/Protein/Entries', '/Protein/DescriptionBuffer') as de_man, \
             GeneOntologyManager(self.h5, '/Annotations/GeneOntology', '/Ontologies/GO') as go_man:
            xref_importer = XRefImporter(db_parser, gs, xref_tab, ec_tab, go_man, de_man)
            files = self.xref_databases()
            dbs_iter = fileinput.input(files=files)
            db_parser.parse_entrytags(dbs_iter)
            xref_importer.flush_buffers()
            xref_importer.build_suffix_index()

    def add_group_metadata(self):
        m = OmaGroupMetadataLoader(self.h5)
        m.add_data()

    def close(self):
        self.h5.root._f_setattr('conversion_end', time.strftime("%c"))
        self.h5.close()
        self.logger.info('closed {}'.format(self.h5.filename))

    def create_indexes(self):
        self.logger.info('creating indexes for HogLevel table')
        hogTab = self.h5.get_node('/HogLevel')
        for col in ('Fam', 'ID', 'Level'):
            if not hogTab.colindexed[col]:
                hogTab.colinstances[col].create_csindex()
        orthoxmlTab = self.h5.get_node('/OrthoXML/Index')
        orthoxmlTab.cols.Fam.create_csindex()

        self.logger.info('creating missing indexes for Entries table')
        entryTab = self.h5.get_node('/Protein/Entries')
        for col in ('EntryNr', 'OmaHOG', 'OmaGroup', 'MD5ProteinHash'):
            if not entryTab.colindexed[col]:
                entryTab.colinstances[col].create_csindex()

        self.logger.info('creating index for xrefs (EntryNr and XRefId)')
        xrefTab = self.h5.get_node('/XRef')
        xrefTab.cols.EntryNr.create_csindex()
        xrefTab.cols.XRefId.create_csindex()

        self.logger.info('creating index for go (EntryNr and TermNr)')
        goTab = self.h5.get_node('/Annotations/GeneOntology')
        goTab.cols.EntryNr.create_csindex()
        goTab.cols.TermNr.create_index()

        self.logger.info('creating index for EC (EntryNr)')
        ec_tab = self.h5.get_node('/Annotations/EC')
        ec_tab.cols.EntryNr.create_csindex()

        self.logger.info('creating index for domains (EntryNr)')
        domtab = self.h5.get_node('/Annotations/Domains')
        domtab.cols.EntryNr.create_csindex()

        self.logger.info('creating indexes for HOG to prevalent domains '
                         '(Fam and DomainId)')
        dom2hog_tab = self.h5.get_node('/HOGAnnotations/Domains')
        dom2hog_tab.cols.DomainId.create_csindex()
        domprev_tab = self.h5.get_node('/HOGAnnotations/DomainArchPrevalence')
        domprev_tab.cols.Fam.create_csindex()

    def _iter_canonical_xref(self):
        """extract one canonical xref id for each protein.

        We take the first valid xref per gene with the ordering of xrefsources
        as given in the xrefsource_order."""
        xrefsource_order = ('UniProtKB/SwissProt', 'UniProtKB/TrEMBL',
                            'Ensembl Gene', 'Ensembl Protein', 'FlyBase',
                            'WormBase', 'EnsemblGenomes', 'RefSeq', 'SourceID')

        xrefs = self.h5.get_node('/XRef')
        source_enum = xrefs.get_enum('XRefSource')
        canonical_sources = [source_enum[z] for z in xrefsource_order]
        current_protein = None
        past_proteins = set([])
        for xref in xrefs:
            if xref['EntryNr'] != current_protein:
                if current_protein:
                    past_proteins.add(current_protein)
                    yield (current_protein, current_xref[1])
                current_protein = xref['EntryNr']
                current_xref = (1000, b'')  # init with a sentinel
                if current_protein in past_proteins:
                    raise DataImportError('Data in /XRef is not grouped w.r.t. EntryNr')
            try:
                rank = canonical_sources.index(xref['XRefSource'])
                if rank < current_xref[0]:
                    current_xref = (rank, xref['XRefId'])
            except ValueError:
                pass
        if current_protein:
            yield (current_protein, current_xref[1])

    def add_canonical_id(self):
        """add one canonical xref id to the /Protein/Entries table."""
        self.logger.info('adding canonical ids for each protein...')
        prot_tab = self.h5.get_node('/Protein/Entries')
        canonical_ids = numpy.chararray(shape=(len(prot_tab),), itemsize=prot_tab.cols.CanonicalId.dtype.itemsize)
        for eNr, canonical_id in self._iter_canonical_xref():
            row_nr = eNr - 1
            row = prot_tab[row_nr]
            if row['EntryNr'] != eNr:
                self.logger.warn('Entries table not properly sorted: {}, expected {}'.format(row['EntryNr'], eNr))
                raise DataImportError('Entries table not properly sorted')
            canonical_ids[row_nr] = canonical_id
        prot_tab.modify_column(0, len(prot_tab), 1, column=canonical_ids, colname='CanonicalId')
        prot_tab.flush()

    def add_domain_info(self, domains):
        self.logger.info('adding domain information...')
        domtab = self.h5.create_table('/Annotations', 'Domains', tablefmt.DomainTable, createparents=True,
                                      expectedrows=1e7)
        entrytab = self.h5.get_node('/Protein/Entries')
        md5_to_enr = collections.defaultdict(list)
        for e in entrytab:
            md5_to_enr[e['MD5ProteinHash']].append(e['EntryNr'])

        buffer = []
        for i, domain in enumerate(domains):
            for entry_nr in md5_to_enr[domain.md5.encode('utf-8')]:
                buffer.append((entry_nr, domain.id, domain.coords))
                if len(buffer) > 5000:
                    domtab.append(buffer)
                    buffer = []
            if i % 50000 == 0:
                self.logger.info('processed {:d} domain annotations so far'.format(i))
        if len(buffer) > 0:
            domtab.append(buffer)
        domtab.flush()

    def add_domainname_info(self, domainname_infos):
        self.logger.info('adding domain name information...')
        dom_name_tab = self.h5.create_table('/Annotations', 'DomainDescription', tablefmt.DomainDescriptionTable,
                                            createparents=True, expectedrows=2e5)
        buffer = []
        for i, dom_info in enumerate(domainname_infos):
            buffer.append(dom_info)
            if len(buffer) > 5000:
                self._write_to_table(dom_name_tab, buffer)
                buffer = []
            if i % 50000 == 0:
                self.logger.info('processed {:d} domain name descriptions so far'.format(i))
        if len(buffer) > 0:
            self._write_to_table(dom_name_tab, buffer)
        dom_name_tab.flush()

    def update_summary_stats(self):
        """update the summary statistics of xrefs & go.

        The function analyses the well-known xref sources as well as
        GO annotations and computes aggregated counts for
        all / in OMA Group / in HOGs for all of them.
        """
        for tab_name, sum_fun in [('/Annotations/GeneOntology', self.count_xref_summary),
                                  ('/XRef', self.count_xref_summary)]:
            summary = sum_fun()
            tab = self.h5.get_node(tab_name)
            for attr, val in summary.items():
                tab.set_attr(attr, val)

        group_sizes = self.collect_group_sizes()
        summary = self._get_or_create_node('/Summary', 'Various Summary Statistics')
        for group_type in group_sizes.keys():
            grp_size_tab = self.create_table_if_needed(
                summary, '{}_size_hist'.format(group_type),
                description=tablefmt.GroupsizeHistogram,
                drop_data=True)
            data = sorted(group_sizes[group_type].items())
            grp_size_tab.append(data)

        cov_fracs = self.add_domain_covered_sites_counts()
        cov_hist, bins = numpy.histogram(cov_fracs[cov_fracs > 0], bins=numpy.linspace(0, 1, 51))
        cov_hist_data = numpy.zeros(50, dtype=[('BinEndValue', 'f4'), ('Counts', 'i4')])
        cov_hist_data['BinEndValue'] = bins[1:]
        cov_hist_data['Counts'] = cov_hist
        dom_cov_hist_tab = self.create_table_if_needed(summary, 'Domain_coverage_hist',
                 drop_data=True, obj=cov_hist_data)
        dom_cov_hist_tab.set_attr('frac_genes_w_domain', len(cov_fracs[cov_fracs > 0]) / len(cov_fracs))
        dom_cov_hist_tab.set_attr('mean_coverage_overall', numpy.mean(cov_fracs))
        dom_cov_hist_tab.set_attr('mean_coverage_w_domain', numpy.mean(cov_fracs[cov_fracs > 0]))

    def count_gene_ontology_summary(self):
        self.logger.info('Bulding gene ontology annotations summary info')
        go_tab = self.h5.get_node('/Annotations/GeneOntology')
        prot_tab = self.h5.get_node('/Protein/Entries')
        exp_codes = frozenset([b'EXP', b'IDA', b'IPI', b'IMP', b'IGI' b'IEP'])
        cnts = collections.Counter()
        cur_enr = None
        for (enr, term), row_iter in itertools.groupby(go_tab, operator.itemgetter('EntryNr','TermNr')):
            evidences = {row['Evidence'] for row in row_iter}
            is_iea = b'IEA' in evidences
            evidences.discard(b'IEA')
            is_exp = not exp_codes.isdisjoint(evidences)
            is_cur = len(evidences.difference(exp_codes)) > 0
            cnts['annotations_any'] += 1
            if is_exp:
                cnts['annotations_exp'] += 1
            if is_cur:
                cnts['annotations_currated'] += 1
            if is_iea:
                cnts['annotations_iea'] += 1
            if cur_enr != enr:
                e = next(prot_tab.where('EntryNr == {}'.format(enr))).fetch_all_fields()
                cnts['proteins_any'] += 1
                if e['OmaGroup'] != 0:
                    cnts['protein_OmaGroup'] += 1
                if len(e['OmaHOG']) > 0:
                    cnts['protein_HOG'] += 1
                cur_enr = enr
        return cnts

    def count_xref_summary(self):
        self.logger.info('Building cross-ref summary info')
        xref_tab = self.h5.get_node('/XRef')
        prot_tab_iter = iter(self.h5.get_node('/Protein/Entries'))
        source = xref_tab.get_enum('XRefSource')
        trusted = frozenset(['UniProtKB/SwissProt', 'UniProtKB/TrEMBL', 'RefSeq', 'EntrezGene', 'Ensembl Gene', 'Ensembl Protein'])
        if len(trusted.difference(source._names.keys())) > 0:
            raise ValueError('set of trusted xrefs is invalid')
        cnts = collections.Counter()

        entry = next(prot_tab_iter)
        for enr, xref_it in itertools.groupby(xref_tab, operator.itemgetter('EntryNr')):
            while entry['EntryNr'] < enr:
                entry = next(prot_tab_iter)
            sources_all = [source._values[x['XRefSource']] for x in xref_it]
            cnts += collections.Counter(sources_all)
            has_trusted_xref = len(trusted.intersection(sources_all)) > 0
            if has_trusted_xref:
                cnts['trusted_all'] += 1
                if entry['OmaGroup'] != 0:
                    cnts['trusted_OmaGroup'] += 1
                if len(entry['OmaHOG']) > 0:
                    cnts['trusted_HOG'] += 1
        return cnts

    def collect_group_sizes(self):
        self.logger.info("Building grouping size histograms")
        groupings = ('OmaHOG', 'OmaGroup')
        memb_cnts = {grp: collections.defaultdict(int) for grp in groupings}
        fam_re = re.compile(br'([A-Z]+:)?(?P<fam>[0-9]+).*')
        prot_tab = self.h5.get_node('/Protein/Entries')
        for row in prot_tab:
            for grp in groupings:
                if grp == 'OmaHOG':
                    m = fam_re.match(row[grp])
                    if m is None:
                        continue
                    grp_id = int(m.group('fam'))
                else:
                    grp_id = int(row[grp])
                    if grp_id == 0:
                        continue
                memb_cnts[grp][grp_id] += 1
        sizes = {grp: collections.defaultdict(int) for grp in groupings}
        for grp in groupings:
            for grp_size in memb_cnts[grp].values():
                sizes[grp][grp_size] += 1
        return sizes

    def compute_domaincovered_sites(self):
        dom_tab = self.h5.get_node('/Annotations/Domains')
        domains = pandas.DataFrame.from_records(dom_tab[:])

        def dlen(coords):
            doms = [int(pos) for pos in coords.split(b':')]
            return sum((doms[i + 1] - doms[i] + 1 for i in range(0, len(doms), 2)))

        # sum all parts of each domain region and store total length in DLen column
        domains = domains.assign(DLen=domains['Coords'].apply(dlen))
        # sum over all domains per protein
        cov_sites = domains.groupby('EntryNr').agg({'DLen': sum})
        return cov_sites

    def add_domain_covered_sites_counts(self):
        """Stores the number of AA covered by a DomainAnnotation.

        This method adds to the hdf5 file a /Protein/DomainCoverage array that
        contains the number of AA sites covered by a domain. The position
        corresponds to the protein entry numbers in /Protein/Entries.

        :Note: The method assumes that the domains are all non-overlapping.
            If they are not, the reported coverage will be too high!

        :return: covered fractions by domains for each protein
        :rtype: numpy.array"""
        self.logger.info("Counting covered sites by domains")
        cov_sites_df = self.compute_domaincovered_sites()

        prot_tab = self.h5.get_node('/Protein/Entries')
        enr_col = prot_tab.col('EntryNr')
        assert numpy.all(numpy.equal(enr_col, numpy.arange(1, len(prot_tab)+1)))

        cov_sites = numpy.zeros(len(prot_tab), dtype=numpy.uint32)
        for eNr, coverage in zip(cov_sites_df.index, cov_sites_df.DLen.values):
            cov_sites[eNr-1] = coverage
        create_node = False
        try:
            dom_cov_tab = self.h5.get_node('/Protein/CoveredSitesByDomains')
            if len(dom_cov_tab) != len(cov_sites):
                self.h5.remove_node('/Protein/CoveredSitesByDomains')
                create_node = True
        except tables.NoSuchNodeError:
            create_node = True
        if create_node:
            dom_cov_tab = self.h5.create_carray('/Protein', 'CoveredSitesByDomains',
                                                tables.UInt32Atom(), (len(cov_sites),))
        dom_cov_tab[0:len(cov_sites)] = cov_sites
        return cov_sites / (prot_tab.col('SeqBufferLength') - 1)

    def add_sequence_suffix_array(self, k=6, fn=None, sa=None):
        '''
            Adds the sequence suffix array to the database. NOTE: this
            (obviously) requires A LOT of memory for large DBs.
        '''
        # Ensure we're run in correct order...
        assert ('Protein' in self.h5.root), 'Add proteins before calc. SA!'
        idx_compr = tables.Filters(complevel=6, complib='blosc', fletcher32=True)

        # Add to separate file if fn is set.
        if fn is None:
            db = self.h5
        else:
            fn = os.path.normpath(os.path.join(
                    os.getenv('DARWIN_BROWSERDATA_PATH', ''),
                    fn))
            db = tables.open_file(fn, 'w', filters=idx_compr)
            db.create_group('/', 'Protein')
            db.root._f_setattr('conversion_start', time.strftime("%c"))
            self.logger.info('opened {}'.format(db.filename))

        # Load sequence buffer to memory - this is required to calculate the SA.
        # Do it here (instead of in PySAIS) so that we can use it for computing
        # the split points later.
        seqs = self.h5.get_node('/Protein/SequenceBuffer')[:].tobytes()
        n = len(self.h5.get_node('/Protein/Entries'))

        # Compute & save the suffix array to DB. TODO: work out what compression
        # works best!
        if sa is None:
            sa = sais(seqs)
            sa[:n].sort()  # Sort delimiters by position.
        db.create_carray('/Protein',
                         name='SequenceIndex',
                         title='concatenated protein sequences suffix array',
                         obj=sa,
                         filters=idx_compr)

        # Create lookup table for fa2go
        dtype = (numpy.uint32 if (n < numpy.iinfo(numpy.uint32).max) else
                 numpy.uint64)
        idx = numpy.zeros(sa.shape, dtype=dtype)
        mask = numpy.zeros(sa.shape, dtype=numpy.bool)

        # Compute mask and entry index for sequence buff
        for i in range(n):
            s = (sa[i - 1] if i > 0 else -1) + 1
            e = (sa[i] + 1)
            idx[s:e] = i + 1
            mask[(e - k):e] = True  # (k-1) invalid and delim.

        # Mask off those we don't want...
        sa = sa[~mask[sa]]

        # Reorder the necessary elements of entry index
        idx = idx[sa]

        # Initialise lookup array
        atom = (tables.UInt32Atom if dtype is numpy.uint32 else tables.UInt64Atom)
        kmers = KmerEncoder(k, is_protein=True)
        kmer_lookup_arr = db.create_vlarray('/Protein',
                              name='KmerLookup',
                              atom=atom(shape=()),
                              title='kmer entry lookup table',
                              filters=idx_compr,
                              expectedrows=len(kmers))
        kmer_lookup_arr._f_setattr('k', k)

        # Now find the split points and construct lookup ragged array.
        ii = 0
        for kk in tqdm(range(len(kmers)), desc='Constructing kmer lookup'):
            kmer = kmers.encode(kk)
            if (ii < len(sa)) and (seqs[sa[ii]:(sa[ii] + k)] == kmer):
                jj = ii + 1
                while (jj < len(sa)) and (seqs[sa[jj]:(sa[jj] + k)] == kmer):
                    jj += 1
                kmer_lookup_arr.append(idx[ii:jj])
                # New start
                ii = jj
            else:
                # End or not found
                kmer_lookup_arr.append([])

        if db.filename != self.h5.filename:
            self.logger.info('storing external links to SequenceIndex and KmerLookup')
            self.h5.create_external_link('/Protein', 'KmerLookup',
                                         self._relative_path_to_external_node(kmer_lookup_arr))
            self.h5.create_external_link('/Protein', 'SequenceIndex',
                                         self._relative_path_to_external_node(db.root.Protein.SequenceIndex))
            db.root._f_setattr('conversion_end', time.strftime("%c"))
            db.close()
            self.logger.info('closed {}'.format(db.filename))

    def _relative_path_to_external_node(self, node):
        rel_path = os.path.relpath(node._v_file.filename, os.path.dirname(self.h5.filename))
        return str(rel_path + ":" + node._v_pathname)

    def add_hog_domain_prevalence(self):
        # Check that protein entries / domains are added already to the DB
        assert True  # TODO

        # Used later
        hl_tab = self.h5.get_node('/HogLevel')
        if not hl_tab.colindexed['Fam']:
            hl_tab.colinstances['Fam'].create_csindex()

        # Load the HOG -> Entry table to memory
        prot_tab = self.h5.root.Protein.Entries
        # TODO: work out how to do this in a neater way
        df = pandas.DataFrame.from_records(((z['EntryNr'], z['OmaHOG'], z['SeqBufferLength'])
                                        for z in prot_tab.iterrows()),
                                       columns=['EntryNr', 'OmaHOG', 'SeqBufferLength'])
        # Strip singletons
        df = df[~(df['OmaHOG'] == b'')]

        # Reformat HOG ID to plain-integer for top-level grouping only
        df['OmaHOG'] = df['OmaHOG'].apply(lambda i: int(i[4:].split(b'.')[0]))

        # Load domains
        domains = pandas.DataFrame.from_records(self.h5.root.Annotations.Domains[:])

        # Ensure sorted by coordinate - TODO: move this to DA import function
        domains['start'] = domains['Coords'].apply(lambda c:
                                                   int(c.split(b':')[0]))
        domains.sort_values(['EntryNr', 'start'], inplace=True)
        domains = domains[['EntryNr', 'DomainId']]

        # Merge domains / entry-hog tables. Keep entries with no domains
        # so that we can count the size of the HOGs.
        df = pandas.merge(df, domains, on='EntryNr', how='left')

        # Gather entry-domain for each HOG.
        hog2dom = []
        hog2info = []
        for (hog_id, hdf) in tqdm(df.groupby('OmaHOG')):
            size = len(set(hdf['EntryNr']))

            hdf = hdf[~hdf['DomainId'].isnull()]
            cov = len(set(hdf['EntryNr']))  # Coverage with any DA

            if (size > 2) and (cov > 1):
                # There are some annotations
                da = collections.defaultdict(list)
                for (enum, edf) in hdf.groupby('EntryNr'):
                    d = edf['DomainId']
                    d = tuple(d) if (type(d) != bytes) else (d,)
                    da[d].append(enum)

                da = sorted(da.items(), key=lambda i: len(i[1]), reverse=True)
                c = len(da[0][1])  # Count of prev. DA
                if c > 1:
                    # DA exists in more than one member.
                    cons_da = da[0][0]
                    repr_entry = da[0][1][0]
                    tl = hl_tab.read_where('Fam == {}'.format(hog_id))[0]['Level'].decode('ascii')
                    rep_len = hdf[hdf['EntryNr'] == repr_entry]['SeqBufferLength']
                    rep_len = int(rep_len if len(rep_len) == 1 else list(rep_len)[0])

                    # Save the consensus DA
                    off = len(hog2info)  # Offset in the information table.
                    hog2dom += [(off, d) for d in cons_da]

                    # Save required information about this group for the web
                    # view.
                    hog2info.append((hog_id,      # HOG ID
                                     repr_entry,  # Repr. entry
                                     rep_len,     # Repr. entry length
                                     tl,          # Top level of HOG
                                     size,        # HOG size
                                     c))          # Prevalence

        # Create tables in file -- done this way as these end up being pretty
        # small tables (<25MB)
        tab = self.h5.create_table('/HOGAnnotations',
                                   'DomainArchPrevalence',
                                   tablefmt.HOGDomainArchPrevalenceTable,
                                   createparents=True,
                                   expectedrows=len(hog2info))
        self._write_to_table(tab, hog2info)
        tab.flush()  # Required?

        # HOG <-> Domain table
        tab = self.h5.create_table('/HOGAnnotations',
                                   'Domains',
                                   tablefmt.HOGDomainPresenceTable,
                                   createparents=True,
                                   expectedrows=len(hog2dom))
        self._write_to_table(tab, hog2dom)
        tab.flush()  # Required?


def download_url_if_not_present(url, force_copy=False):
    if url.startswith('file://') and not force_copy:
        fname = url[len('file://'):]
        if os.path.exists(fname):
            common.package_logger.info('using file "{}" directly from source without copying.'.format(url))
            return fname
    tmpfolder = os.path.join(os.getenv('DARWIN_NETWORK_SCRATCH_PATH', '/tmp'), "Browser", "xref")
    basename = url.split('/')[-1]
    fname = os.path.join(tmpfolder, basename)
    if not os.path.exists(tmpfolder):
        os.makedirs(tmpfolder)
    if not os.path.exists(fname):
        common.package_logger.info("downloading {} into {}".format(url, fname))
        try:
            urllib.request.urlretrieve(url, fname)
        except urllib.request.URLError:
            common.package_logger.warn('cannot download {}'.format(url))
    return fname


def iter_domains(url):
    DomainTuple = collections.namedtuple('DomainTuple', ('md5', 'id', 'coords'))

    fname = download_url_if_not_present(url)
    with gzip.open(fname, 'rt') as uncompressed:
        dialect = csv.Sniffer().sniff(uncompressed.read(4096))
        uncompressed.seek(0)
        csv_reader = csv.reader(uncompressed, dialect)
        col_md5, col_id, col_coord = (None,) * 3
        coord_fromat_trans = str.maketrans('-,', '::')

        for lineNr, row in enumerate(csv_reader):
            if col_md5 is None:
                # identify which tuples to use.
                if len(row) >= 9:
                    # representative_proteins format. use columns 5-7
                    col_md5, col_id, col_coord = 4, 5, 6
                elif len(row) == 3:
                    # additionally created ones, minimal format
                    col_md5, col_id, col_coord = 0, 1, 2
                else:
                    raise DataImportError("Unknown Domain Annotation format in {}".format(uncompressed.filename))
            try:
                dom = DomainTuple(row[col_md5], row[col_id], row[col_coord].translate(coord_fromat_trans))
                if lineNr < 10:
                    # do some sanity checks on the first few lines
                    if re.match(r'[0-9a-f]{32}$', dom.md5) is None:
                        raise DataImportError("md5 hash of line {:d} has unexpected values: {}"
                                              .format(lineNr, dom.md5))
                    if re.match(r'([1-4]\.\d+\.\d+\.\d+|PF\d+)$', dom.id) is None:
                        raise DataImportError("Domain-ID of line {:d} has unexpected value: {}"
                                              .format(lineNr, dom.id))
                    if re.match(r'\d+:\d+', dom.coords) is None:
                        raise DataImportError("Domain coordinates in line {:d} has unexpected value: {}"
                                              .format(lineNr, dom.coords))
                yield dom
            except Exception:
                common.package_logger.exception('cannot create tuple from line {}'.format(lineNr))


def only_pfam_or_cath_domains(iterable):
    cath_re = re.compile(r'[1-4]\.')
    for dom in iterable:
        if dom.id.startswith('PF') or cath_re.match(dom.id) is not None:
            yield dom


def filter_duplicated_domains(iterable):
    """filter duplicated domain annotations that come from different proteins
    with the exact same sequence."""
    seen = set([])
    ignored = 0
    for dom in iterable:
        if not dom in seen:
            seen.add(dom)
            yield dom
        else:
            ignored += 1
    common.package_logger.info("skipped {} duplicated domains. {} distinct domains yielded"
                               .format(ignored, len(seen)))


class OmaGroupMetadataLoader(object):
    """OMA Group Meta data extractor.

    This class provides the means to import the Keywords and Fingerprints
    of the OMA Groups into the hdf5 database. The data is stored under
    in the node defined by :attr:`meta_data_path`, which defaults to
    /OmaGroups/MetaData.
    """
    keyword_name = "Keywords.drw"
    finger_name = "Fingerprints"

    meta_data_path = '/OmaGroups/MetaData'

    def __init__(self, db):
        self.db = db

    def add_data(self):
        common.package_logger.info('adding OmaGroup Metadata')
        nr_groups = self._get_nr_of_groups()
        has_meta_data = self._check_textfiles_avail()
        if has_meta_data:
            data = self._load_data()
            fingerprints = data['Fingerprints']
            keywords = data['Keywords']
        else:
            common.package_logger.warning('No fingerprint nor keyword information available')
            fingerprints = [b'n/a'] * nr_groups
            keywords = [b''] * nr_groups
        if nr_groups != len(fingerprints) or nr_groups != len(keywords):
            raise DataImportError('nr of oma groups does not match the number of fingerprints and keywords')

        grptab, keybuf = self._create_db_objects(nr_groups)
        self._fill_data_into_db(fingerprints, keywords, grptab, keybuf)
        grptab.modify_column(column=self._get_group_member_counts(), colname='NrMembers')
        self._create_indexes(grptab)

    def _create_db_objects(self, nrows):
        key_path = os.path.join(os.path.dirname(self.meta_data_path), 'KeywordBuffer')
        try:
            self.db.get_node(self.meta_data_path)
            self.db.remove_node(self.meta_data_path)
            self.db.remove_node(key_path)
        except tables.NoSuchNodeError:
            pass
        root, name = self.meta_data_path.rsplit('/', 1)
        grptab = self.db.create_table(root, name, tablefmt.OmaGroupTable,
                                      expectedrows=nrows, createparents=True)
        buffer = self.db.create_earray(root, "KeywordBuffer", tables.StringAtom(1), (0,),
                                       'concatenated group keywords  descriptions',
                                              expectedrows=500 * nrows)
        return grptab, buffer

    def _fill_data_into_db(self, stable_ids, keywords, grp_tab, key_buf):
        row = grp_tab.row
        buf_pos = 0
        for i in range(len(stable_ids)):
            row['GroupNr'] = i+1
            row['Fingerprint'] = stable_ids[i]
            row['KeywordOffset'] = buf_pos
            row['KeywordLength'] = len(keywords[i])
            row.append()
            key = numpy.ndarray((len(keywords[i]),), buffer=keywords[i],
                                dtype=tables.StringAtom(1))
            key_buf.append(key)
            buf_pos += len(keywords[i])
        grp_tab.flush()
        key_buf.flush()

    def _create_indexes(self, grp_tab):
        grp_tab.cols.Fingerprint.create_csindex()
        grp_tab.cols.GroupNr.create_csindex()

    def _parse_darwin_string_list_file(self, fh):
        data = fh.read()
        start, end = data.find(b'['), data.rfind(b', NULL]')
        if end == -1:
            end = data.rfind(b']:')
        part = data[start:end] + b']'
        as_json = part.replace(b"''", b"__apos__").replace(b"'", b'"')\
                      .replace(b'__apos__', b"'")
        as_list = json.loads(as_json.decode())
        return [el.encode('utf8') for el in as_list]

    def _load_data(self):
        return callDarwinExport('GetGroupData()')

    def _get_nr_of_groups(self):
        etab = self.db.get_node('/Protein/Entries')
        try:
            return etab[etab.colindexes['OmaGroup'][-1]]['OmaGroup']
        except KeyError:
            return max(etab.col('OmaGroup'))

    def _get_group_member_counts(self):
        grp_nr, cnts = numpy.unique(self.db.get_node('/Protein/Entries').col('OmaGroup'), return_counts=True)
        if grp_nr[0] == 0:
            cnts = cnts[1:]
        assert(len(cnts) == self._get_nr_of_groups())
        return cnts

    def _check_textfiles_avail(self):
        rootdir = os.getenv('DARWIN_BROWSERDATA_PATH','')
        fn1 = os.path.join(rootdir, self.keyword_name)
        fn2 = os.path.join(rootdir, self.finger_name)
        return os.path.exists(fn1) and os.path.exists(fn2)


class DescriptionManager(object):
    def __init__(self, db, entry_path, buffer_path):
        self.db = db
        self.entry_path = entry_path
        self.buffer_path = buffer_path

    def __enter__(self):
        self.entry_tab = self.db.get_node(self.entry_path)
        if not numpy.all(numpy.equal(self.entry_tab.col('EntryNr'),
                                     numpy.arange(1, len(self.entry_tab) + 1))):
            raise RuntimeError('entry table is not sorted')

        root, name = os.path.split(self.buffer_path)
        self.desc_buf = self.db.create_earray(root, name,
                                              tables.StringAtom(1), (0,), 'concatenated protein descriptions',
                                              expectedrows=len(self.entry_tab) * 100)
        self.cur_eNr = None
        self.cur_desc = []
        bufindex_dtype = numpy.dtype([(col, self.entry_tab.coldtypes[col])
                                      for col in ('DescriptionOffset', 'DescriptionLength')])
        # columns to be stored in entry table with buffer index data
        self.buf_index = numpy.zeros(len(self.entry_tab), dtype=bufindex_dtype)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.cur_eNr:
            self._store_description()
        self.desc_buf.flush()
        self.entry_tab.modify_columns(columns=self.buf_index,
                                      names=self.buf_index.dtype.names)
        self.entry_tab.flush()

    def add_description(self, eNr, desc):
        """stages a description for addition. Note that the descriptions
        must be ordered according to the entryNr, i.e. all descriptions
        related to eNr X must be staged before changeing to another eNr."""
        if self.cur_eNr and self.cur_eNr != eNr:
            self._store_description()
            self.cur_desc = []
        self.cur_eNr = eNr
        self.cur_desc.append(desc)

    def _store_description(self):
        buf = "; ".join(self.cur_desc).encode('utf-8')
        buf = buf[0:2 ** 16 - 1]  # limit to max value of buffer length field
        len_buf = len(buf)
        idx = self.cur_eNr - 1
        self.buf_index[idx]['DescriptionOffset'] = len(self.desc_buf)
        self.buf_index[idx]['DescriptionLength'] = len_buf
        self.desc_buf.append(numpy.ndarray((len_buf,), buffer=buf, dtype=tables.StringAtom(1)))


class GeneOntologyManager(object):
    ontology_url = "http://purl.obolibrary.org/obo/go/go-basic.obo"

    def __init__(self, db, annotation_path, ontology_path):
        self.db = db
        self.annotation_path = annotation_path
        self.ontology_path = ontology_path
        self._go_buf = []
        self.quote_re = re.compile(r'([[,])([\w_:]+)([,\]])')

    def __enter__(self):
        go_obo_file = download_url_if_not_present(self.ontology_url)
        # check that ontology file is not broken. if we can build it, it should be ok
        self.go = GeneOntology(OntologyParser(go_obo_file))
        self.go.parse()

        with open(go_obo_file, 'rb') as fh:
            go_obo = fh.read()
        root, name = os.path.split(self.ontology_path)
        obo = self.db.create_carray(root, name, title='Gene ontology hierarchy definition', createparents=True,
                              obj=numpy.ndarray(len(go_obo), buffer=go_obo, dtype=tables.StringAtom(1)))
        obo._f_setattr('ontology_release', self._get_obo_version(obo))

        root, name = os.path.split(self.annotation_path)
        self.go_tab = self.db.create_table(root, name, tablefmt.GeneOntologyTable,
                                      'Gene Ontology annotations', expectedrows=1e8, createparents=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._flush_buffers()
        self.go_tab.flush()

    def _get_obo_version(self, obo_arr):
        header = obo_arr[0:1000].tobytes()
        rel_info = re.search(b'data-version:\s*(?P<version>[\w/_ -]+)', header)
        if rel_info is not None:
            rel_info = rel_info.group('version').decode()
        return rel_info

    def _flush_buffers(self):
        common.package_logger.info('flushing go annotations buffers')
        if len(self._go_buf) > 0:
            self.go_tab.append(self._go_buf)
        self._go_buf = []

    def add_annotations(self, enr, gos):
        """parse go annotations and add them to the go buffer"""
        if not (isinstance(enr, int) and isinstance(gos, str)):
            raise ValueError('input data invalid')
        for t in gos.split('; '):
            t = t.strip()
            try:
                term, rem = t.split('@')
            except ValueError as e:
                common.package_logger.warning('cannot parse GO annotation: ' + t)
                continue

            try:
                term_nr = self.go.term_by_id(term).id
            except ValueError:
                common.package_logger.warning('invalid GO term for entry {:d}: {:s} (likely obsolete)'
                                              .format(enr, term))
                continue
            rem = rem.replace('{', '[')
            rem = rem.replace('}', ']')
            rem = self.quote_re.sub('\g<1>"\g<2>"\g<3>', rem)
            for evi, refs in eval(rem):
                for ref in refs:
                    self._go_buf.append((enr, term_nr, evi, ref.encode('utf-8')))
            if len(self._go_buf) > 2e6:
                self._flush_buffers()


class GroupAnnotatorInclGeneRefs(familyanalyzer.GroupAnnotator):
    def _annotateGroupR(self, node, og, idx=0):
        if familyanalyzer.OrthoXMLQuery.is_geneRef_node(node):
            node.set('og', og)
        else:
            super()._annotateGroupR(node, og, idx)


class HogConverter(object):
    def __init__(self, entry_tab):
        self.fam_re = re.compile(r'HOG:(?P<fam_nr>\d+)')
        self.hogs = numpy.zeros(shape=(len(entry_tab) + 1,), dtype=entry_tab.cols.OmaHOG.dtype)
        self.entry_tab = entry_tab

    def attach_newick_taxonomy(self, tree):
        self.taxonomy = familyanalyzer.NewickTaxonomy(tree)

    def _assert_hogid_has_correct_prefix(self, fa_parser):
        for grp in fa_parser.getToplevelGroups():
            if not grp.get('id').startswith('HOG:'):
                grp.set('id', 'HOG:{:07d}'.format(int(grp.get('id'))))

    def convert_file(self, fn):
        p = familyanalyzer.OrthoXMLParser(fn)
        self._assert_hogid_has_correct_prefix(p)
        if hasattr(self, 'taxonomy'):
            p.augmentTaxonomyInfo(self.taxonomy)
        else:
            p.augmentTaxonomyInfo(familyanalyzer.TaxonomyFactory.newTaxonomy(p))
        GroupAnnotatorInclGeneRefs(p).annotateDoc()

        levs = []
        for fam in p.getToplevelGroups():
            m = self.fam_re.match(fam.get('og'))
            fam_nr = int(m.group('fam_nr'))
            levs.extend([(fam_nr, n.getparent().get('og'), n.get('value'),) + self.get_hog_scores(n.getparent())
                         for n in p._findSubNodes('property', root=fam)
                         if n.get('name') == "TaxRange"])

        geneNodes = p.root.findall('.//{{{ns0}}}geneRef'.
                                   format(**familyanalyzer.OrthoXMLParser.ns))
        for x in geneNodes:
            self.hogs[int(x.get('id'))] = x.get('og')

        return levs

    def write_hogs(self):
        """update the Entry Table with the newly collected OmaHOG values for all
        the proteins at once.

        .. note: This method will overwrite any previous value of the OmaHOG column"""
        self.entry_tab.modify_column(0, len(self.entry_tab), 1, self.hogs[1:], 'OmaHOG')
        self.entry_tab.flush()

    def get_hog_scores(self, og_node):
        """extract the scores associated with an orthologGroup node

        only scores that are defined in HOGsTable are extract. The method
        returns a tuple with the scores in the order of the score fields."""
        scores = collections.OrderedDict([(score, tablefmt.HOGsTable.columns[score].dflt)
                                          for score in ('CompletenessScore', 'ImpliedLosses')])
        for score in og_node.iterfind('{*}score'):
            score_id = score.get("id")
            if score_id == "CompletenessScore":
                scores['CompletenessScore'] = float(score.get('value'))
            elif score_id == "ImpliedLosses":
                scores['ImpliedLosses'] = int(score.get('value'))
        return tuple(scores.values())


class XRefImporter(object):
    """Object to import various types of crossreferences into hdf5.

    The XRefImporter registers at a db_parser object various handlers
    to import the various types of xrefs, namely ids, go-terms,
    EC annotations and descriptions."""
    def __init__(self, db_parser, genomes_tab, xref_tab, ec_tab, go_manager, desc_manager):
        self.xrefs = []
        self.ec = []
        self.xref_tab = xref_tab
        self.ec_tab = ec_tab
        self.go_manager = go_manager
        self.desc_manager = desc_manager

        self.verif_enum = tablefmt.XRefTable.columns.get('Verification').enum
        xrefEnum = tablefmt.XRefTable.columns.get('XRefSource').enum
        tag_to_enums = {
            'GI': (xrefEnum['GI'], 'exact'),
            'EntrezGene': (xrefEnum['EntrezGene'], 'exact'),
            'WikiGene': (xrefEnum['WikiGene'], 'unchecked'),
            'IPI': (xrefEnum['IPI'], 'unchecked'),
            'Refseq_ID': (xrefEnum['RefSeq'], 'exact'),
            'SwissProt': (xrefEnum['UniProtKB/SwissProt'], 'exact'),
            'GeneName': (xrefEnum['Gene Name'], 'unchecked'),
            'ORFNames': (xrefEnum['ORF Name'], 'unchecked'),
            'OrderedLocusNames': (xrefEnum['Ordered Locus Name'], 'unchecked'),
            'ProtName': (xrefEnum['Protein Name'], 'unchecked'),
            'Synonyms': (xrefEnum['Synonym'], 'unchecked'),
            'HGNC_Id': (xrefEnum['HGNC'], 'unchecked'),
            'PMP': (xrefEnum['PMP'], 'exact'),
            'PDB': (xrefEnum['PDB'], 'unchecked'),
            'EMBL': (xrefEnum['EMBL'], 'unchecked'),
            'ID': (xrefEnum['SourceID'], 'exact'),
            'AC': (xrefEnum['SourceAC'], 'exact'),
        }
        for tag, enumval in tag_to_enums.items():
            db_parser.add_tag_handler(
                tag,
                lambda key, enr, typ=enumval: self.multi_key_handler(key, enr, typ[0], typ[1]))
        db_parser.add_tag_handler('DE',
                                  lambda key, enr: self.description_handler(key, enr))
        db_parser.add_tag_handler('GO', self.go_handler)
        db_parser.add_tag_handler('ID', self.assign_source_handler)
        db_parser.add_tag_handler('AC', self.assign_source_handler)
        db_parser.add_tag_handler('EC', self.ec_handler)

        for tag in ['SwissProt_AC', 'UniProt']:  # UniProt/TrEMBL tag is cut to UniProt!
            db_parser.add_tag_handler(tag,
                                      lambda key, enr, typ=xrefEnum['UniProtKB/TrEMBL']:
                                      self.remove_uniprot_code_handler(key, enr, typ))

        # register the potential_flush as end_of_entry_notifier
        db_parser.add_end_of_entry_notifier(self.potential_flush)

        self.db_parser = db_parser
        self.xrefEnum = xrefEnum
        self.ENS_RE = re.compile(r'ENS(?P<species>[A-Z]{0,3})(?P<typ>[GTP])(?P<num>\d{11})')
        self.FB_RE = re.compile(r'FB(?P<typ>[gnptr]{2})(?P<num>\d{7})')
        self.NCBI_RE = re.compile(r'[A-Z]{3}\d{5}\.\d$')
        self.WB_RE = re.compile(r'WBGene\d{8}$')
        self.EC_RE = re.compile(r'\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)')
        self.ENSGENOME_RE = re.compile(b'Ensembl (Metazoa|Plant|Fungi|Protist|Bacteria)', re.IGNORECASE)

        self.FLUSH_SIZE = 5e6

        # info about current genome
        self.genomes_tab = genomes_tab
        self._cur_genome = None

    def _get_genome_info(self, entry_nr):
        if not (self._cur_genome is not None and self._cur_genome['EntryOff'] < entry_nr <=
                self._cur_genome['EntryOff'] + self._cur_genome['TotEntries']):
            self._cur_genome = self.genomes_tab[self.genomes_tab['EntryOff'].searchsorted(entry_nr+1)-1]
        return self._cur_genome

    def from_EnsemblGenome(self, entry_nr):
        genome_info = self._get_genome_info(entry_nr)
        return self.ENSGENOME_RE.search(genome_info['Release']) is not None

    def flush_buffers(self):
        common.package_logger.info('flushing xrefs and ec buffers')
        if len(self.xrefs) > 0:
            self.xref_tab.append(sorted(uniq(self.xrefs)))
            self.xrefs = []
        if len(self.ec) > 0:
            self.ec_tab.append(sorted(uniq(self.ec)))
            self.ec = []

    def potential_flush(self):
        if len(self.xrefs) > self.FLUSH_SIZE:
            self.flush_buffers()

    def _add_to_xrefs(self, eNr, enum_nr, key, verif='unchecked'):
        if not isinstance(eNr, int):
            raise ValueError('eNr is of wrong type:' + str(eNr))
        self.xrefs.append((eNr, enum_nr, key.encode('utf-8'), self.verif_enum[verif], ))

    def key_value_handler(self, key, eNr, enum_nr, verif='unchecked'):
        """basic handler that simply adds a key (the xref) under a given enum_nr"""
        self._add_to_xrefs(eNr, enum_nr, key, verif)

    def multi_key_handler(self, multikey, eNr, enum_nr, verif='unchecked'):
        """try to split the myltikey field using '; ' as a delimiter and add each
        part individually under the passed enum_nr id type."""
        for key in multikey.split('; '):
            if key.startswith('Rep'):
                continue
            pos = key.find('.Rep')
            if pos > 0:
                key = key[0:pos]
            self._add_to_xrefs(eNr, enum_nr, key, verif)

    def assign_source_handler(self, multikey, eNr):
        """handler that splits the multikey field at '; ' locations and
        tries to guess for each part the id_type. If a type could be
        identified, it is added under with this id type, otherwise left out."""
        for key in multikey.split('; '):
            ens_match = self.ENS_RE.match(key)
            if ens_match is not None:
                typ = ens_match.group('typ')
                if typ == 'P':
                    enum_nr = self.xrefEnum['Ensembl Protein']
                elif typ == 'G':
                    enum_nr = self.xrefEnum['Ensembl Gene']
                elif typ == 'T':
                    enum_nr = self.xrefEnum['Ensembl Transcript']
                common.package_logger.debug(
                    'ensembl: ({}, {}, {})'.format(key, typ, enum_nr))
                self._add_to_xrefs(eNr, enum_nr, key, 'exact')

            for enum, regex in {'FlyBase': self.FB_RE, 'NCBI': self.NCBI_RE, 'WormBase': self.WB_RE}.items():
                match = regex.match(key)
                if match is not None:
                    enum_nr = self.xrefEnum[enum]
                    self._add_to_xrefs(eNr, enum_nr, key, 'unchecked')
            if self.from_EnsemblGenome(eNr):
                self._add_to_xrefs(eNr, self.xrefEnum.EnsemblGenomes, key, 'exact')

    def go_handler(self, gos, enr):
        self.go_manager.add_annotations(enr, gos)

    def ec_handler(self, ecs, enr):
        for t in ecs.split('; '):
            t = t.strip()
            acc_match = self.EC_RE.match(t)
            if acc_match is not None:
                self.ec.append((enr, acc_match.group(0)))

    def description_handler(self, de, eNr):
        self.desc_manager.add_description(eNr, de)

    def remove_uniprot_code_handler(self, multikey, eNr, enum_nr):
        """remove the species part (sep by '_') of a uniprot long accession to the short acc"""
        common.package_logger.debug(
            'remove_uniprot_code_handler called ({}, {},{})'.format(multikey, eNr, enum_nr))
        for key in multikey.split('; '):
            pos = key.find('_')
            if pos > 0:
                self._add_to_xrefs(eNr, enum_nr, key[0:pos], 'exact')
            else:
                self._add_to_xrefs(eNr, enum_nr, key, 'exact')

    def build_suffix_index(self, force=False):
        parent, name = os.path.split(self.xref_tab._v_pathname)
        file_ = self.xref_tab._v_file
        idx_node = get_or_create_tables_node(file_, os.path.join(parent, "{}_Index".format(name)))
        for arr_name, typ in (('buffer', tables.StringAtom(1)), ('offset', tables.UInt32Atom())):
            try:
                n = idx_node._f_get_child(arr_name)
                if not force:
                    raise tables.NodeError("Suffix index for xrefs does already exist. Use 'force' to overwrite")
                n.remove()
            except tables.NoSuchNodeError:
                pass
            file_.create_earray(idx_node, arr_name, typ, (0,), expectedrows=100e6)
        buf, off = (idx_node._f_get_child(node) for node in ('buffer', 'offset'))
        self._build_lowercase_xref_buffer(buf, off)
        sa = sais(buf)
        try:
            idx_node._f_get_child('suffix').remove()
        except tables.NoSuchNodeError:
            pass
        file_.create_carray(idx_node, 'suffix', obj=sa)

    def _build_lowercase_xref_buffer(self, buf, off):
        cur_pos = 0
        for xref_row in tqdm(self.xref_tab):
            lc_ref = xref_row['XRefId'].lower()
            ref = numpy.ndarray((len(lc_ref),), buffer=lc_ref, dtype=tables.StringAtom(1))
            buf.append(ref)
            off.append([cur_pos])
            cur_pos += len(lc_ref)


class DarwinDbEntryParser:
    def __init__(self):
        """Initializes a Parser for SGML formatted darwin database file
        """
        self.tag_handlers = collections.defaultdict(list)
        self.end_of_entry_notifier = []

    def add_tag_handler(self, tag, handler):
        """add a callback handler for a certain tag"""
        self.tag_handlers[tag].append(handler)
        common.package_logger.debug('# handlers for {}: {}'.format(tag, len(self.tag_handlers[tag])))

    def add_end_of_entry_notifier(self, handler):
        self.end_of_entry_notifier.append(handler)

    def parse_entrytags(self, fh):
        """ AC, CHR, DE, E, EMBL, EntrezGene, GI, GO, HGNC_Name, HGNC_Sym,
        ID, InterPro, LOC, NR , OG, OS, PMP, Refseq_AC, Refseq_ID, SEQ,
        SwissProt, SwissProt_AC, UniProt/TrEMBL, WikiGene, flybase_transcript_id

        :param fh: an already opened file handle to the darwin database
                   file to be parsed."""
        eNr = 0
        for line in fh:
            line = line.strip()
            if not line.startswith('<E>'):
                common.package_logger.debug('skipping line:' + line)
                continue

            eNr += 1
            common.package_logger.debug('entry {}: {}'.format(eNr, line.encode('utf-8')))
            entry = lxml.html.fragment_fromstring(line)
            for tag, handlers in self.tag_handlers.items():
                common.package_logger.debug('tag {} ({} handlers)'.format(tag, len(handlers)))
                tag_text = [t.text for t in entry.findall('./' + tag.lower())]
                for value in tag_text:
                    # common.package_logger.debug('value of tag: {}'.format(value.encode('utf-8')))
                    if value is None:
                        continue
                    for handler in handlers:
                        handler(value, eNr)
                        # common.package_logger.debug('called handler {} with ({},{})'.format(
                        #    handler, value.encode('utf-8'), eNr))
            for notifier in self.end_of_entry_notifier:
                notifier()


DomainDescription = collections.namedtuple('DomainDescription',
                                           tables.dtype_from_descr(tablefmt.DomainDescriptionTable).names)


class CathDomainNameParser(object):
    re_pattern = re.compile(r'(?P<id>[0-9.]*)\s{3,}\w{7}\s{3,}:\s*(?P<desc>.*)')
    source = b'CATH/Gene3D'

    def __init__(self, url):
        self.fname = download_url_if_not_present(url)

    def parse(self):
        open_lib = gzip.open if self.fname.endswith('.gz') else open
        with open_lib(self.fname, 'rt') as fh:
            for line in fh:
                match = self.re_pattern.match(line)
                if match is not None:
                    yield DomainDescription(DomainId=match.group('id').encode('utf-8'),
                                            Source=self.source,
                                            Description=match.group('desc').encode('utf-8'))


class PfamDomainNameParser(CathDomainNameParser):
    re_pattern = re.compile(r'(?P<id>\w*)\t\w*\t\w*\t\w*\t(?P<desc>.*)')
    source = b'Pfam'


def augment_genomes_json_download_file(fpath, h5, backup='.bak'):
    """Augment the genomes.json file in the download section with additional info

    This function stores the ncbi taxonomy identifiers of internal nodes and adds
    the number of ancestral genes to the internal nodes.

    :param fpath: path to genomes.json file
    :param h5: hdf5 database handle."""
    common.package_logger.info("Augmenting genomes.json file with Nr of HOGs per level")
    # load nr of ancestral genomes at each level
    ancestral_hogs = collections.Counter()
    step = 2**15
    hog_tab = h5.get_node('/HogLevel')
    for start in range(0, len(hog_tab), step):
        ancestral_hogs.update((l.decode() for l in hog_tab.read(start, stop=start+step, field='Level')))
    # load taxonomy and sorter by Name
    tax = h5.get_node('/Taxonomy').read()
    sorter = numpy.argsort(tax['Name'])
    with open(fpath, 'rt') as fh:
        genomes = json.load(fh)
    os.rename(fpath, fpath + '.bak')

    def traverse(node):
        if 'children' not in node:
            return
        for child in node['children']:
            traverse(child)
        try:
            node['nr_hogs'] = ancestral_hogs[node['name']]
        except KeyError as e:
            common.package_logger.warning('no ancestral hog counts for '+node['name'])
            node['nr_hogs'] = 0

        try:
            n = node['name'].encode('utf-8')
            idx = numpy.searchsorted(tax['Name'], n, sorter=sorter)
            if tax['Name'][sorter[idx]] == n:
                node['taxid'] = int(tax['NCBITaxonId'][sorter[idx]])
            else:
                raise ValueError('not in taxonomy: {}'.format(n))
        except Exception:
            common.package_logger.exception('Cannot identify taxonomy id')

    traverse(genomes)
    with open(fpath, 'wt') as fh:
        json.dump(genomes, fh)


def getLogger(level='DEBUG'):
    import logging

    log = logging.getLogger('pyoma')
    if isinstance(level, str):
        level = logging.getLevelName(level.upper())
        if not isinstance(level, int):
            level = logging.DEBUG
    log.setLevel(level)
    logHandler = logging.StreamHandler()
    logHandler.setLevel(level)
    logHandler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    log.addHandler(logHandler)
    return log


def main(name="OmaServer.h5", k=6, idx_name=None, domains=None, log_level='INFO'):
    idx_name = (name + '.idx') if idx_name is None else idx_name

    log = getLogger(log_level)
    x = DarwinExporter(name, logger=log)
    x.add_version()
    x.add_species_data()
    x.add_orthologs()
    x.add_same_species_relations()
    x.add_proteins()
    x.add_hogs()
    x.add_xrefs()
    x.add_synteny_scores()
    x.add_homoeology_confidence()
    if domains is None:
        domains = ["file://dev/null"]
    x.add_domain_info(filter_duplicated_domains(only_pfam_or_cath_domains(itertools.chain(
        iter_domains('ftp://orengoftp.biochem.ucl.ac.uk/gene3d/CURRENT_RELEASE/' +
                     'representative_uniprot_genome_assignments.csv.gz'),
        iter_domains('file://{}/additional_domains.mdas.csv.gz'.format(os.getenv('DARWIN_BROWSERDATA_PATH', '')))
    ))))
    x.add_domainname_info(itertools.chain(
        CathDomainNameParser('http://download.cathdb.info/cath/releases/latest-release/'
                             'cath-classification-data/cath-names.txt').parse(),
        PfamDomainNameParser('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz').parse()))
    x.add_canonical_id()
    x.add_group_metadata()
    x.add_hog_domain_prevalence()
    x.close()

    x = DarwinExporter(name, logger=log)
    x.create_indexes()
    x.add_sequence_suffix_array(k=k, fn=idx_name)
    x.update_summary_stats()

    genomes_json_fname = os.path.normpath(os.path.join(
        os.path.dirname(name), '..', 'downloads', 'genomes.json'))
    augment_genomes_json_download_file(genomes_json_fname, x.h5)
    x.close()
