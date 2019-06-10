import pandas
import tables
import logging
try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x, **kwargs: x
logger = logging.getLogger(__name__)


class SyntenyScorer(object):
    def __init__(self, h5_handle, genome, windowsize=10):
        self.h5_handle = h5_handle
        self.genome = genome
        self.windowsize = windowsize
        if isinstance(h5_handle, tables.File):
            self.h5_handle = h5_handle
        elif isinstance(h5_handle, (str, bytes)):
            self.h5_handle = tables.open_file(h5_handle, 'r')
        else:
            raise TypeError("expected h5_handle to be either h5-file handle or a path to file")

        genome_row = next(self.h5_handle.root.Genome.where('UniProtSpeciesCode == genome'))
        self.genome_range = (int(genome_row['EntryOff']) + 1,
                             int(genome_row['EntryOff'] + genome_row['TotEntries']))
        genome_df = pandas.DataFrame(self.h5_handle.root.Protein.Entries.read_where(
            '(EntryNr >= {}) & (EntryNr <= {})'.format(*self.genome_range)))
        self.genome_df = genome_df[(genome_df['AltSpliceVariant'] == 0) | (genome_df['AltSpliceVariant'] == genome_df['EntryNr'])]
        self.genome_df.reset_index(inplace=True)
        self.relations_df = self._load_pairwise_relations()

    def _load_pairwise_relations(self):
        df = pandas.DataFrame(self.h5_handle.get_node('/PairwiseRelation/{}/within'.format(self.genome)).read_where('RelType == 5'))
        return df[['EntryNr1', 'EntryNr2', 'SyntenyConservationLocal']]

    def get_neighbor_genes(self, query):
        q = self.genome_df[self.genome_df['EntryNr'] == query]
        if len(q) == 0:
            logger.error("querying neighbor genes for non-primary variant (EntryNr: {})".format(query))
            return []
        query_chr = q['Chromosome']
        neighbor = self.genome_df[max(0, q.index.item() - self.windowsize // 2): q.index.item() + self.windowsize//2 + 1]
        return neighbor[neighbor['Chromosome'] == query_chr.item()]

    def score_of_pair(self, entry1, entry2):
        neigh1 = self.get_neighbor_genes(entry1)
        neigh2 = self.get_neighbor_genes(entry2)
        if len(neigh1) <= 1 or len(neigh2) <= 1:
            raise TooSmallChromosome("too few genes on chromosome: {}, {}".format(len(neigh1), len(neigh2)))

        rels_among_windows = self.relations_df[
            (self.relations_df['EntryNr1'] >= neigh1.iloc[0]['EntryNr']) &
            (self.relations_df['EntryNr1'] <= neigh1.iloc[-1]['EntryNr']) &
            (self.relations_df['EntryNr2'] >= neigh2.iloc[0]['EntryNr']) &
            (self.relations_df['EntryNr2'] <= neigh2.iloc[-1]['EntryNr'])]

        score1 = (len(set(rels_among_windows['EntryNr1'])) - 1) / (len(neigh1) - 1)
        score2 = (len(set(rels_among_windows['EntryNr2'])) - 1) / (len(neigh2) - 1)
        res = {'entry_nr1': int(entry1), 'entry_nr2': int(entry2),
               'chr1': neigh1.iloc[0]['Chromosome'].decode(),
               'chr2': neigh2.iloc[0]['Chromosome'].decode(),
               'nr_genes_window1': len(neigh1), 'nr_genes_window2': len(neigh2),
               'synteny_score_1': score1, 'synteny_score_2': score2,
               'mean_synteny_score': (score1 + score2) / 2}
        return res

    def compute_scores(self):
        res = []
        for idx, rel in tqdm(self.relations_df.iterrows(), total=len(self.relations_df)):
            try:
                res.append(self.score_of_pair(rel['EntryNr1'], rel['EntryNr2']))
            except TooSmallChromosome as e:
                logging.info("Skipping {}/{}: {}".format(int(rel['EntryNr1']), int(rel['EntryNr2']), e))
                pass
        return pandas.DataFrame(res)


class TooSmallChromosome(Exception):
    pass


#### MAIN ####
if __name__ == "__main__":
    import argparse
    # Get arguments from command line
    parser = argparse.ArgumentParser(
        description='Returns windows and their proportion of homoeolog VPs for a given chromosome')
    parser.add_argument('--h5', help='name of h5 file, full path', required=True)
    parser.add_argument('--window_genes', help='window size in genes', default=10)
    parser.add_argument('--genome', help='5 letter code of polyploid genome to analyze')
    parser.add_argument('--outfile', help='name where results will be stored (file name created to include parameters)', \
                        default="synteny_results.tsv")

    args = parser.parse_args()
    h5file_path = args.h5
    logging.basicConfig(level=logging.INFO)

    scorer = SyntenyScorer(tables.open_file(h5file_path), args.genome)
    data = scorer.compute_scores()
    columns = ['entry_nr1', 'chr1', 'nr_genes_window1', 'entry_nr2', 'chr2', 'nr_genes_window2', 'synteny_score_1',
               'synteny_score_2', 'mean_synteny_score']
    data[columns].to_csv(args.outfile, sep='\t', header=True, index=True)
