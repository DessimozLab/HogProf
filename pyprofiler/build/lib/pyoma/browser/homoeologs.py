import pandas
import logging
import collections
import numpy as np
import matplotlib
matplotlib.use('agg')
from skfuzzy import control as ctrl
from skfuzzy import gaussmf
import sklearn
import sklearn.preprocessing
import tables

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x, **kwargs: x
logger = logging.getLogger(__name__)


def define_universe(df):
    # New Antecedent/Consequent objects hold universe variables and membership functions
    distance = ctrl.Antecedent(np.arange(0, df['Distance'].max() + .01, .01), 'distance')
    synteny = ctrl.Antecedent(np.arange(0, 1.01, .01), 'synteny_score')
    total_nb_hom = ctrl.Antecedent(np.arange(2, df['TotalCopyNr'].max() + 1, 1), 'total_nb_homoeologs')
    conf = ctrl.Consequent(np.arange(0, 101, 1), 'conf')
    return distance, synteny, total_nb_hom, conf


def create_fuzzy_rules(distance, synteny, total_nb_hom, conf):
    """Takes the antecedent and consequent objects as input"""

    # very low confidence
    rule1 = ctrl.Rule(synteny['low'] & distance['high'] & total_nb_hom['high'], conf['very_low'])

    # low confidence
    rule2 = ctrl.Rule((synteny['low'] & distance['high'] & total_nb_hom['low']) |
                      (synteny['low'] & distance['high'] & total_nb_hom['med']) |
                      (synteny['low'] & distance['med'] & total_nb_hom['high']) |
                      (synteny['med'] & distance['high'] & total_nb_hom['high']),
                      conf['low'])

    # medium confidence
    rule3 = ctrl.Rule((synteny['high'] & distance['high'] & total_nb_hom['high']) |

                      (synteny['low'] & distance['med'] & total_nb_hom['low']) |
                      (synteny['low'] & distance['med'] & total_nb_hom['med']) |
                      (synteny['med'] & distance['high'] & total_nb_hom['low']) |
                      (synteny['med'] & distance['high'] & total_nb_hom['med']) |
                      (synteny['med'] & distance['med'] & total_nb_hom['high']) |
                      (synteny['med'] & distance['med'] & total_nb_hom['low']) |
                      (synteny['med'] & distance['med'] & total_nb_hom['med']) |
                      (synteny['low'] & distance['low'] & total_nb_hom['low']) |
                      (synteny['low'] & distance['low'] & total_nb_hom['med']) |
                      (synteny['low'] & distance['low'] & total_nb_hom['high']),
                      conf['med'])

    # high confidence
    rule4 = ctrl.Rule((synteny['high'] & distance['high'] & total_nb_hom['low']) |
                      (synteny['high'] & distance['high'] & total_nb_hom['med']) |
                      (synteny['high'] & distance['low'] & total_nb_hom['high']) |
                      (synteny['high'] & distance['low'] & total_nb_hom['med']) |
                      (synteny['high'] & distance['med'] & total_nb_hom['high']) |
                      (synteny['high'] & distance['med'] & total_nb_hom['low']) |
                      (synteny['high'] & distance['med'] & total_nb_hom['med']) |
                      (synteny['med'] & distance['low'] & total_nb_hom['high']) |
                      (synteny['med'] & distance['low'] & total_nb_hom['med']) |
                      (synteny['med'] & distance['low'] & total_nb_hom['low']),
                      conf['high'])

    # very high confidence
    rule5 = ctrl.Rule(synteny['high'] & distance['low'] & total_nb_hom['low'],
                      conf['very_high'])
    return [rule1, rule2, rule3, rule4, rule5]


def get_distance_mf(df, distance):
    # here, the first numnber is the universe, second number the central point, third the standard deviation
    distance['low'] = gaussmf(distance.universe, 0, (df['Distance'].max() / 10))

    distance['med'] = gaussmf(distance.universe,
                                   (df['Distance'].max() / 4),
                                   (df['Distance'].max() / 10))

    distance['high'] = gaussmf(distance.universe,
                                    df['Distance'].max(),
                                    (df['Distance'].max() / 2.5))
    return distance


def get_synteny_mf(df, synteny, view=False):
    # synteny (gaussian)
    synteny['low'] = gaussmf(synteny.universe, 0, .15)
    synteny['med'] = gaussmf(synteny.universe, .3, .15)
    synteny['high'] = gaussmf(synteny.universe, .7, .25)
    return synteny


def get_total_nb_hom_mf(df, total_nb_hom):
    copy_nr_median = df['TotalCopyNr'].median()
    total_nb_hom['low'] = gaussmf(total_nb_hom.universe,
                                       copy_nr_median, copy_nr_median)

    total_nb_hom['med'] = gaussmf(total_nb_hom.universe,
                                       4 * copy_nr_median,
                                       1.5 * copy_nr_median)

    total_nb_hom['high'] = gaussmf(total_nb_hom.universe,
                                        df['TotalCopyNr'].max(),
                                        df['TotalCopyNr'].max() / 2.5)
    return total_nb_hom


def get_conf_mf(df, conf):
    # confidence (gaussian)
    conf['very_low'] = gaussmf(conf.universe, 0, 20)
    conf['low'] = gaussmf(conf.universe, 50, 10)
    conf['med'] = gaussmf(conf.universe, 70, 10)
    conf['high'] = gaussmf(conf.universe, 90, 10)
    conf['very_high'] = gaussmf(conf.universe, 100, 10)
    return conf


def get_conf_score(simulation, input_dic):
    """This function takes the simulation and outputs confidence score
    'input_dic' is a dictionary of the inputs for a homoeolog pair"""

    simulation.inputs(input_dic)
    simulation.compute()
    return simulation.output['conf']


class HomeologsConfidenceCalculator(object):
    def __init__(self, h5_handle, genome):
        self.h5_handle = h5_handle
        self.genome = genome
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
        self.genome_df = genome_df[
            (genome_df['AltSpliceVariant'] == 0) | (genome_df['AltSpliceVariant'] == genome_df['EntryNr'])]
        self.genome_df.reset_index(inplace=True)
        self.relations_df = self._load_pairwise_relations()

    def _load_pairwise_relations(self):
        """load the homoeologous relations of the cannonical splice variants only
        The method returns a pandas dataframe with the relations."""
        df = pandas.DataFrame(
            self.h5_handle.get_node('/PairwiseRelation/{}/within'.format(self.genome)).read_where('RelType == 5'))
        df = df[df['EntryNr1'].isin(self.genome_df['EntryNr']) & df['EntryNr2'].isin(self.genome_df['EntryNr'])]
        return df[['EntryNr1', 'EntryNr2', 'SyntenyConservationLocal', 'Distance']]

    def _count_homeologs_per_entry(self, df):
        return collections.Counter(df['EntryNr1'])

    def _augment_dataframe_with_all_features(self, df):
        counts = self._count_homeologs_per_entry(df)
        df['TotalCopyNr'] = df.apply(lambda x: counts[x['EntryNr1']] + counts[x['EntryNr2']], axis=1)
        df.loc[df.SyntenyConservationLocal < 0, 'SyntenyConservationLocal'] = 0
        return df

    def calculate_scores(self):
        # load dataframe
        df = self.relations_df
        df = self._augment_dataframe_with_all_features(df)

        distanceObj, syntenyObj, total_nb_homObj, confObj = define_universe(df)
        distance = get_distance_mf(df, distanceObj)
        synteny = get_synteny_mf(df, syntenyObj)
        total_nb_hom = get_total_nb_hom_mf(df, total_nb_homObj)
        conf = get_conf_mf(df, confObj)

        # create simulation
        rules = create_fuzzy_rules(distance, synteny, total_nb_hom, conf)
        control_system = ctrl.ControlSystem(rules)
        simulation = ctrl.ControlSystemSimulation(control_system)

        def defuzzify(row):
            return get_conf_score(simulation,
                                  {'distance': row['Distance'],
                                   'synteny_score': row['SyntenyConservationLocal'],
                                   'total_nb_homoeologs': row['TotalCopyNr']})

        df['fuzzy_confidence'] = df.apply(defuzzify, axis=1)

        # scale the confidence between minimum value and 100
        min_max_scaler = sklearn.preprocessing.MinMaxScaler(feature_range=(df['fuzzy_confidence'].min(), 100))
        df['fuzzy_confidence_scaled'] = min_max_scaler.fit_transform(df['fuzzy_confidence'].values.reshape(-1, 1))
        return df


class HomeologsConfidenceCalculatorFromTSV(HomeologsConfidenceCalculator):
    def __init__(self, infile):
        self.relations_df = pandas.read_csv(infile, sep='\t')
        expected_columns = ['EntryNr1', 'EntryNr2', 'SyntenyConservationLocal', 'Distance']
        if len(set(expected_columns) - set(self.relations_df.columns.values)) > 0:
            raise KeyError("provided inputfile does not have all expected columns. "
                           "Expected columns are {}".format(expected_columns))


if __name__ == "__main__":
    import argparse
    # Get arguments from command line
    parser = argparse.ArgumentParser(
        description='Computes homoeology confidence score using fuzzy logic')
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('--h5', help="name of hdf5 file, full path")
    grp.add_argument('--csv', help="tab-separated file with input data as alternative to hdf5 file")
    parser.add_argument('--genome',
                        help="5 letter code of polyploid genome to analyze. "
                             "Must be specified if used with --h5 option.")
    parser.add_argument('--outfile',
                        help="name where results will be stored (file name created to include parameters)",
                        default="homoeolog_confidence.tsv")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    if args.h5 is not None and args.genome is None:
        import sys
        sys.stderr.write("genomes argument required if using with an hdf5 file as input")
        sys.exit(1)

    if args.h5:
        import tables
        scorer = HomeologsConfidenceCalculator(tables.open_file(args.h5), args.genome)
    else:
        scorer = HomeologsConfidenceCalculatorFromTSV(args.csv)
    data = scorer.calculate_scores()
    data.to_csv(args.outfile, sep='\t', header=True, index=True)
