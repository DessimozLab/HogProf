from .convert import *
from pyoma.browser import OrthoXMLSplitter
import os


class StandaloneExporter(DarwinExporter):
    DRW_CONVERT_FILE = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")

    def __init__(self, root, name, **kwargs):
        os.environ['DARWIN_BROWSERDATA_PATH'] = os.path.abspath(root)
        super(StandaloneExporter, self).__init__(name, **kwargs)
        self.transformed = False
        self.cache_dir = os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), 'pyoma')

    def add_homologs(self):
        self.assert_cached_results()
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'homologs' not in rel_node_for_genome:
                pass

    def get_version(self):
        # TODO: obtain real version
        return "OmaStandalone; 1.0.x"

    def assert_cached_results(self):
        if not self.transformed:
            res = self.call_darwin_export("TransformDataToCache('{}');".format(
                    self.cache_dir))
            if res != 'success':
                raise DarwinException('could not transform data from darwin')
            self.transformed = True
            os.environ['DARWIN_NETWORK_SCRATCH_PATH'] = os.getenv('DARWIN_BROWSERDATA_PATH')

    def add_orthologs(self):
        self.assert_cached_results()
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'VPairs' not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                    'pyoma', 'vps', '{}.txt.gz'.format(genome))
                if os.path.exists(cache_file):
                    data = load_tsv_to_numpy((cache_file, 0, 0, False,))
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export('GetVPsForGenome({})'.format(genome))

                vp_tab = self.h5.create_table(rel_node_for_genome, 'VPairs', tablefmt.PairwiseRelationTable,
                                              expectedrows=len(data))
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                self._write_to_table(vp_tab, data)
                vp_tab.cols.EntryNr1.create_csindex()

    def add_hogs(self):
        hog_path = os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'], 'Output')

        entryTab = self.h5.get_node('/Protein/Entries')

        tree_filename = os.path.join(
            os.environ['DARWIN_BROWSERDATA_PATH'],
            'EstimatedSpeciesTree.nwk')

        hog_converter = HogConverter(entryTab)

        if os.path.exists(tree_filename):
            hog_converter.attach_newick_taxonomy(tree_filename)

        fn = 'HierarchicalGroups.orthoxml'

        # Split the OrthoXML up (puts in cache_dir/split_hog).
        hog_cache_dir = os.path.join(self.cache_dir, 'split_hogs')
        ortho_splitter = OrthoXMLSplitter.OrthoXMLSplitter(os.path.join(hog_path, fn), cache_dir=hog_cache_dir)
        ortho_splitter()

        hogTab = self.h5.create_table('/', 'HogLevel', tablefmt.HOGsTable,
                                      'nesting structure for each HOG', expectedrows=1e8)
        self.orthoxml_buffer = self.h5.create_earray('/OrthoXML', 'Buffer',
                                                     tables.StringAtom(1), (0,), 'concatenated orthoxml files',
                                                     expectedrows=1e9, createparents=True)
        self.orthoxml_index = self.h5.create_table('/OrthoXML', 'Index', tablefmt.OrthoXmlHogTable,
                                                   'Range index per HOG into OrthoXML Buffer', expectedrows=5e6)

        try:
            levels = hog_converter.convert_file(os.path.join(hog_path, fn))
            hogTab.append(levels)
            fam_nrs = set([z[0] for z in levels])
            for fam_nr in fam_nrs:
                hog_fn = "HOG{:06d}.orthoxml".format(fam_nr)
                self.add_orthoxml(os.path.join(hog_cache_dir, hog_fn), [fam_nr])
        except Exception as e:
            self.logger.error('an error occured while processing ' + fn + ':')
            self.logger.exception(e)

        hog_converter.write_hogs()

    def _get_genome_database_paths(self):
        return self.call_darwin_export('GetGenomeFileNames();')

    def xref_databases(self):
        return self._get_genome_database_paths()


def import_oma_run(path, outfile, add_domains=True, log_level='INFO'):
    log = getLogger(log_level)
    x = StandaloneExporter(path, outfile, logger=log, mode='write')
    x.add_version()
    x.add_species_data()
    x.add_orthologs()
    x.add_proteins()
    x.add_hogs()
    x.add_xrefs()
    domain_url = ('ftp://orengoftp.biochem.ucl.ac.uk/gene3d/CURRENT_RELEASE/'+
                  'representative_uniprot_genome_assignments.csv.gz')
    if not add_domains:
        domain_url = 'file:///dev/null'
    x.add_domain_info(only_pfam_or_cath_domains(iter_domains(domain_url)))
    x.add_domainname_info(itertools.chain(
        CathDomainNameParser('http://download.cathdb.info/cath/releases/latest-release/'
                             'cath-classification-data/cath-names.txt').parse(),
        PfamDomainNameParser('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz').parse()))
    x.add_canonical_id()
    x.add_group_metadata()
    x.add_hog_domain_prevalence()
    x.close()

    x = StandaloneExporter(path, outfile, logger=log)
    x.create_indexes()
    x.add_sequence_suffix_array()
    x.update_summary_stats()
    x.close()


if __name__ == "__main__":
    import_oma_run('~/Repositories/OmaStandalone', 'oma.h5')
