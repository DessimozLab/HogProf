#import profiler
print("\nNote: this script requires the hogprof environment\n")
import argparse
import pandas as pd
import os
import sys
import pickle
from time import time
from datetime import datetime
import traceback
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gumbel_r, kstest,  beta, probplot, gaussian_kde
import statsmodels.api as sm
import glob
from statsmodels.distributions.empirical_distribution import ECDF
from ete3 import NCBITaxa
ncbi = NCBITaxa()

### get current time
start_time_secs = time() 
start_time = datetime.fromtimestamp(start_time_secs).strftime('%Y-%m-%d %H:%M:%S')
### print start time human readable
print(f"Script started at: {start_time}\n")

def add_profiler_path(profiler_path):
    """Add the directory containing profiler.py to Python's sys.path."""
    if not os.path.isfile(profiler_path):
        raise FileNotFoundError(
            f"Could not find profiler.py at: {profiler_path}"
        )
    profiler_dir = os.path.dirname(os.path.abspath(profiler_path))
    if profiler_dir not in sys.path:
        sys.path.append(profiler_dir)

def create_directory(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
    return dirpath

'''here expect index of df to have subhog ids and a column taxid'''
def create_bins(hogmetadata_df, outputdir, taxidmapper={}, samplesize=21):
    ### reverse dict
    idmapper = {v: k for k, v in taxidmapper.items()}
    ### add taxid column to df from index prefix
    hogmetadata_df['taxid'] = hogmetadata_df.index.to_series().apply(lambda x: str(x).split('_')[1], None)
    print(hogmetadata_df.head())
    ### sort by taxid - is this problematic for random sampling????????????????
    hogmetadata_df.sort_values(by='taxid', inplace=True) 
    # set seed
    np.random.seed(42)
    # Create a new DataFrame to hold the samples
    pairs_df = pd.DataFrame(columns=['subhog'])
    # Create a new DataFrame to hold the info about the chosen subhogs
    subhogs_samples_info_df = pd.DataFrame(columns=['taxid', 'clade','subhog']) ###used to have proteins_num
    #print(species_bins)
    subhogs_samples_info_list = []
    samples_list = []
    for taxid_bin in hogmetadata_df['taxid'].unique():
        # Get the subhogs in the current species bin
        subhogs_in_bin = hogmetadata_df[hogmetadata_df['taxid'] == taxid_bin]
        # If there are at least 2 subhogs, sample pairs
        if len(subhogs_in_bin) >= 2:
            try:
                sample = subhogs_in_bin.sample(n=samplesize, replace=False) ### sample 21 for allvsall, 200 for pairs
            except ValueError as e:
                print(f"Error sampling for taxid bin {taxid_bin}: {e}")
                sample = subhogs_in_bin
                print(f"  Sampled {len(sample)} subhogs for taxid_bin {taxid_bin}")
            for index, row in sample.iterrows():
                subhogs_samples_info_list.append({
                    'taxid': taxid_bin,
                    'subhog': index
                })
            # Create rows for later
            for i in range(0, len(sample)):
                samples_list.append({
                    'subhog': sample.index[i]
                    })
    # Convert the lists to DataFrames
    subhogs_samples_info_df = pd.DataFrame(subhogs_samples_info_list)
    pairs_df = pd.DataFrame(samples_list)
    print(pairs_df.head())
    ### save to csv files
    subhogs_samples_info_df.to_csv(os.path.join(outputdir, "subhogs_samples_info.csv"), index=True)
    pairs_df.to_csv(os.path.join(outputdir, "subhogs_sampled_pairs.csv"), index=False)
    print(f"Saved subhogs_samples_info_df to {os.path.join(outputdir, 'subhogs_samples_info.csv')}")
    print(f"Saved pairs_df to {os.path.join(outputdir, 'subhogs_sampled_pairs.csv')}")
    return idmapper
    

def get_hog_to_hog_jaccards(subhogs_sampled_pairs_df, lshforestfile, hashes_h5, treefile, 
                            fam2orthoxmlfile, outputdir, profiler_path, allvsall=False, suffix=""):
    ### import profiler
    add_profiler_path(profiler_path)
    import profiler

    ### check for NaN, there should not be any
    if not allvsall:
        if subhogs_sampled_pairs_df['subhog1'].isnull().any() or subhogs_sampled_pairs_df['subhog2'].isnull().any():
            print("Error: Some subhog IDs could not be matched.")
            print(f"Subhogs sampled pairs DataFrame has {subhogs_sampled_pairs_df['subhog1'].isnull().sum()} null subhog1 and {subhogs_sampled_pairs_df['subhog2'].isnull().sum()} null subhog2.")
            print(subhogs_sampled_pairs_df[subhogs_sampled_pairs_df['subhog1'].isnull() | subhogs_sampled_pairs_df['subhog2'].isnull()])
            ### remove those rows
            subhogs_sampled_pairs_df = subhogs_sampled_pairs_df.dropna(subset=['subhog1', 'subhog2'])
            #exit()
    else:
        if subhogs_sampled_pairs_df['subhog'].isnull().any():
            print("Error: Some subhog IDs could not be matched.")
            print(f"Subhogs sampled pairs DataFrame has {subhogs_sampled_pairs_df['subhog'].isnull().sum()} null subhogs")
            print(subhogs_sampled_pairs_df[subhogs_sampled_pairs_df['subhog1'].isnull()])
            print(subhogs_sampled_pairs_df[subhogs_sampled_pairs_df['subhog2'].isnull()])
            ### remove those rows
            subhogs_sampled_pairs_df = subhogs_sampled_pairs_df.dropna(subset=['subhog'])
            #exit()

    ### if already exists, just read
    if os.path.exists(os.path.join(outputdir, f"allvsall_hashmat_{suffix}.npy")):
        print(f"allvsall_hashmat_{suffix}.npy already exists, loading from file.")
        hashmat = np.load(os.path.join(outputdir, f"allvsall_hashmat_{suffix}.npy"))
        return hashmat

    ### read lshforest
    p = profiler.Profiler(lshforestpath = lshforestfile, 
                            hashes_h5= hashes_h5, 
                            #mat_path= fam2orthoxmlfile ,
                            oma = False , 
                            nsamples = 256 ,
                            mastertree = treefile,
                            slicesubhogs = True
                            )

    if not allvsall:
        ### get jaccard similarity for each pair using hog_v_hog(self, hogs)
        print("Calculating Jaccard similarity for each pair of subhogs...")
        jaccard_results = []
        for index, row in subhogs_sampled_pairs_df.iterrows():
            subhog1 = row['subhog1']
            subhog2 = row['subhog2']
            try:
                jaccard_value = p.hog_v_hog([subhog1, subhog2])
                jaccard_results.append({
                    'subhog1': subhog1,
                    'subhog2': subhog2,
                    'jaccard': jaccard_value
                })
            except Exception as e:
                print(f"Error calculating Jaccard for {subhog1} and {subhog2}: {e}")
                jaccard_results.append({
                    'subhog1': subhog1,
                    'subhog2': subhog2,
                    'jaccard': None
                })
        # Convert results to DataFrame
        jaccard_df = pd.DataFrame(jaccard_results)
        # Save to CSV
        jaccard_df.to_csv(os.path.join(outputdir, f"sampled_jaccard_results_{suffix}.csv"), index=False)
        print(f"Jaccard results saved to {os.path.join(outputdir, f'sampled_jaccard_results_{suffix}.csv')}\n")
    else:
        ### use pull_hashes and hog_v_hog method
        print("Calculating Jaccard similarity for all pairs of subhogs...")
        jaccard_results = []
        try:
            subhogs_list = pd.unique(
                pd.concat([subhogs_sampled_pairs_df['subhog1'], subhogs_sampled_pairs_df['subhog2']])
            )
        except:
            subhogs_list = pd.unique(subhogs_sampled_pairs_df['subhog'])
        hashes_dict = p.pull_hashes(subhogs_list)
        hashmat = p.allvall_hashes(hashes_dict)
        # Mask out the diagonal (self-self hits)
        allvsall_flat = hashmat[~np.eye(hashmat.shape[0], dtype=bool)].flatten()
        ### save np table to file
        np.save(os.path.join(outputdir, f"allvsall_hashmat_{suffix}.npy"), allvsall_flat)
        return allvsall_flat

def get_jaccard_threshold(jaccard_values, alpha=0.05):
    """
    Returns the Jaccard similarity value above which p < alpha.
    I.e., the (1 - alpha) quantile of the values.

    Parameters:
        jaccard_values (array-like): Array of Jaccard similarity values.
        alpha (float): Significance level (default 0.05).
    
    Returns:
        float: Jaccard threshold value, or None if input is empty.
    """
    if len(jaccard_values) == 0:
        return None
    return np.quantile(jaccard_values, 1 - alpha)

def get_bins_jaccards_pairwise(samples_info_df, bin_name, bin_dict, subhogs_sampled_pairs_df, lshforestfile, hashes_h5, treefile, outputdir,
                               profiler_path):
    ### import profiler
    add_profiler_path(profiler_path)
    import profiler
    bin_counter = 1
    ### read lshforest
    p = profiler.Profiler(lshforestpath = lshforestfile, 
                            hashes_h5= hashes_h5, 
                            #mat_path= fam2orthoxmlfile ,
                            oma = False , 
                            nsamples = 256 ,
                            mastertree = treefile,
                            slicesubhogs = True
                            )

    ### ensure Profiler has fam_dict
    if not hasattr(p, 'fam_dict') or p.fam_dict is None:
        print("Error: fam_dict is None.")
        print(p.__dict__.keys())
        exit(1)

    ### iterate over every bin large enough
    for species_bin, group_df in samples_info_df.groupby(bin_name):
        binname = species_bin
        bin_counter+=1
        bin_dict[binname] = species_bin
        if os.path.exists(os.path.join(outputdir, f"allvsall_hashmat_{binname}.npy")):   
            continue
        ### get group_df with subhog ids in this bin from subhogs_sampled_pairs_df
        subhogs_in_bin = set(group_df['subhog'])
        # Filter pairs where both subhogs are in this bin
        bin_pairs_df = subhogs_sampled_pairs_df[
            subhogs_sampled_pairs_df['subhog'].isin(subhogs_in_bin)
            ]
        if len(bin_pairs_df) < 10:
            continue
        print(f"Processing bin: {species_bin}")
        ### replacing get_hog_to_hog_jaccards
        jaccard_results = []
        ### split bin_pairs_df in half to make comparisons
        bin_paired_for_comparison_df = turn_subhogs_to_pairs(bin_pairs_df)

        for index, row in bin_paired_for_comparison_df.iterrows():
            subhog1 = row['subhog1']
            subhog2 = row['subhog2']
            try:
                jaccard_value = p.hog_v_hog([subhog1, subhog2])
                jaccard_results.append({
                    'subhog1': subhog1,
                    'subhog2': subhog2,
                    'jaccard': jaccard_value
                })
            except Exception as e:
                print(f"Error calculating Jaccard for {subhog1} and {subhog2}: {e}")
                jaccard_results.append({
                    'subhog1': subhog1,
                    'subhog2': subhog2,
                    'jaccard': None
                })
        # Convert results to DataFrame
        jaccard_df = pd.DataFrame(jaccard_results)
        # Save to CSV
        jaccard_df.to_csv(os.path.join(outputdir, f"sampled_jaccard_results_{binname}.csv"), index=False)
        print(f"Jaccard results saved to {os.path.join(outputdir, f'sampled_jaccard_results_{binname}.csv')}\n")

def turn_subhogs_to_pairs(subhogs_sampled_pairs_df):
    # Shuffle the DataFrame to ensure randomness
    shuffled_df = subhogs_sampled_pairs_df.sample(frac=1, random_state=42).reset_index(drop=True)
    # Create pairs by taking two consecutive rows
    pairs_list = []
    for i in range(0, len(shuffled_df) - 1, 2):
        pairs_list.append({
            'subhog1': shuffled_df.iloc[i]['subhog'],
            'subhog2': shuffled_df.iloc[i + 1]['subhog']
        })
    # Convert the list of pairs to a DataFrame
    pairs_df = pd.DataFrame(pairs_list)
    return pairs_df

def taxid_to_name(taxid, expected_root="Metazoa"):
    if taxid == 'internal_0':
        print(f"Expected root: {expected_root}")
        return f'root_{expected_root}'
    try:
        taxid = int(taxid)
    ### NOTE: maybe this exception should be handled differently
    # for now just return empty string
    except (ValueError, TypeError):
        print(f"Invalid taxid: {taxid}\nUsing ''.")
        return ''
    name_dict = ncbi.get_taxid_translator([taxid])#.translate_to_names([taxid])
    #print(name_dict)
    if len(name_dict) > 1:
        print(f"Warning: More than one name found for taxid {taxid}. Using the first one.")
        #print(name_dict)
    return name_dict[taxid]

def main(hogprofoutputfolder, outputdir, profiler_path, expected_root):

    ### get all variables:
    fam2orthoxmlfile = os.path.join(hogprofoutputfolder, "fam2orthoxml.csv")
    hashes_h5 = os.path.join(hogprofoutputfolder, "hashes.h5")
    lshforestfile = os.path.join(hogprofoutputfolder, "newlshforest.pkl")
    treefile = os.path.join(hogprofoutputfolder, "reformatted_tree.nwk")
    hogmetadatafile = os.path.join(hogprofoutputfolder, "subhogize_table.csv")
    taxaindexfile = os.path.join(hogprofoutputfolder, "idmapper.pkl")
    
    ### print input parameters
    print(f"fam2orthoxmlfile: {fam2orthoxmlfile}")
    print(f"hogmetadatafile: {hogmetadatafile}")
    print(f"lshforestfile: {lshforestfile}")
    print(f"hashes_h5: {hashes_h5}")
    print(f"treefile: {treefile}")
    print(f"taxaindexfile: {taxaindexfile}")
    print(f"outputdir: {outputdir}\n")    

    #print("\nWARNING: This script only works locally for now. Do not try on curnagl!\n")    
    print("Reminder: if this script does not run on curnagl, you need to reinstall in editable mode (pip install -e .)")

    ### create output directory if it does not exist
    outputdir = create_directory(outputdir)

    ### get taxidmapper:
    with open(taxaindexfile, 'rb') as f:
        taxidmapper = pickle.load(f)

    ### get subhogids from fam2orthoxml file
    subhogids_df = pd.read_csv(fam2orthoxmlfile)
    ### if it is an OmaServer run it will not have column labels, so add them (fam, subhog_id)
    if 'fam' not in subhogids_df.columns:
        subhogids_df.columns = ['fam', 'subhog_id'] + list(subhogids_df.columns[2:])
    ### merge first column with subhog_id column 
    subhogids_df['subhogid_full'] = subhogids_df['fam'].astype(str) + "_" + subhogids_df['subhog_id'].astype(str)
    ### else it is an OmaServer run with no column labels
    ### use as index
    subhogids_df.set_index('subhogid_full', inplace=True)
    ### add taxid column (first part of subhog_id)
    subhogids_df['taxid'] = subhogids_df['subhog_id'].astype(str).apply(lambda x: x.split('_')[0])

    ### create bins and sample subhogs
    if not os.path.exists(os.path.join(outputdir, "subhogs_sampled_pairs.csv")):
        idmapper = create_bins(subhogids_df, outputdir, taxidmapper, samplesize=500) ### used to be hogmetadata_df ###500 HOGs, later 250 pairs to create distribution
    else:
        idmapper = {v: k for k, v in taxidmapper.items()}

    ### read samples as pairs
    subhogs_sampled_pairs_df = pd.read_csv(os.path.join(outputdir, "subhogs_sampled_pairs.csv"), index_col=None)
    print(f"Reading sampled pairs")
    print(subhogs_sampled_pairs_df.head())

    ### get bin specific all vs all 
    samples_info_df = pd.read_csv(os.path.join(outputdir, "subhogs_samples_info.csv"), index_col=0)
    
    bin_dict = {}
    bin_name = 'taxid' #'species_bin'
    ### generate jaccard similarities for each bin
    get_bins_jaccards_pairwise(samples_info_df, bin_name, bin_dict, subhogs_sampled_pairs_df, lshforestfile, hashes_h5, treefile, outputdir, profiler_path)

    bin_dict = {}
    bin_name = 'taxid' #'species_bin'
    print(f"Number of unique bins in '{bin_name}': {samples_info_df[bin_name].nunique()}")

    bin_counter = 1
    ### iterate over every bin large enough
    for species_bin, group_df in samples_info_df.groupby(bin_name):
        binname = species_bin
        bin_counter+=1
        bin_dict[binname] = species_bin
        if os.path.exists(os.path.join(outputdir, f"allvsall_hashmat_{binname}.npy")):   
            continue
        ### get group_df with subhog ids in this bin from subhogs_sampled_pairs_df
        subhogs_in_bin = set(group_df['subhog'])
        # Filter pairs where both subhogs are in this bin
        bin_pairs_df = subhogs_sampled_pairs_df[
            subhogs_sampled_pairs_df['subhog'].isin(subhogs_in_bin)
            ]
        if bin_pairs_df.shape[0] < 10:
            continue
        print(f"Processing bin: {species_bin}")
        allvsall_hashmat = get_hog_to_hog_jaccards(bin_pairs_df, lshforestfile, hashes_h5, treefile, 
                                fam2orthoxmlfile, outputdir, profiler_path, allvsall=True, suffix=binname)
        
        ### plot all vs all jaccard distribution
        plt.figure(figsize=(8, 5))
        plt.hist(allvsall_hashmat.flatten(), bins=50, color='purple', alpha=0.7)
        plt.xlabel('Jaccard Similarity')
        plt.ylabel('Frequency')
        plt.title(f'All vs All Jaccard Similarity Distribution {species_bin}')
        plt.tight_layout()
        plt.savefig(os.path.join(outputdir, f"allvsall_jaccard_distribution_{binname}.png"))

    if not os.path.exists(os.path.join(outputdir, f"bins_jaccard_thresholds.csv")):
        #### separately deal with each bin
        allvsall_files = sorted(glob.glob(os.path.join(outputdir, "allvsall_hashmat_*.npy")))
        print(f"Found {len(allvsall_files)} bin-specific all-vs-all files:")
        bins_dicts = []
        for allvsallfile in allvsall_files:
            if 'bin' in os.path.basename(allvsallfile):
                continue
            if not os.path.exists(os.path.join(outputdir, f"bins_jaccard_thresholds.csv")):
                binname = os.path.basename(allvsallfile).split("allvsall_hashmat_")[1][:-4]
                if binname == '':
                    continue
                binnum = int(binname.split("bin")[-1])
                ### read all vs all hashmat
                allvsall_hashmat = np.load(allvsallfile)
                # Remove self hits (diagonal)
                if allvsall_hashmat.ndim == 2:
                    allvsall_flat = allvsall_hashmat[~np.eye(allvsall_hashmat.shape[0], dtype=bool)].flatten()
                else:
                    # Already 1D, nothing to mask
                    allvsall_flat = allvsall_hashmat
                ### calculate threshold for bin
                empirical_t = get_jaccard_threshold(allvsall_flat, alpha=0.01)
                bins_dicts.append({
                            'bin': binname,
                            'bin_range': bin_dict[binnum],
                            'empirical_t': empirical_t,
                            'sampled_comparisons': len(allvsall_flat)
                        })
        ### turn to df
        # Convert results to DataFrame
        bins_thresholds_df = pd.DataFrame(bins_dicts)
        print(bins_thresholds_df.head())
        ### here add idmapper to have taxnames together with the thresholds
        # Assumes taxids were used in the species tree and taxnames in the orthoxml file.
        # This is not always the case.
        bins_thresholds_df['taxid'] = bins_thresholds_df['bin'].map(idmapper)
        ### here add actual taxon names
        bins_thresholds_df['taxname'] = bins_thresholds_df['taxid'].apply(taxid_to_name, expected_root=expected_root)
        ### check here if the assumption was correct - if the column is full of ''
        if bins_thresholds_df['taxname'].isnull().all() or (bins_thresholds_df['taxname'] == '').all():
            # remove that column
            bins_thresholds_df.drop(columns=['taxname'], inplace=True)
            # rename taxid column to taxname
            bins_thresholds_df.rename(columns={'taxid': 'taxname'}, inplace=True)
            # replace NaN with root if bin is 0
            bins_thresholds_df.loc[
                (bins_thresholds_df['bin'].astype(str) == '0') & 
                (bins_thresholds_df['taxname'].isna()),
                'taxname'
            ] = f'root_{expected_root}'
        # Save to CSV
        bins_thresholds_df.to_csv(os.path.join(outputdir, f"bins_jaccard_thresholds.csv"), index=False)
        print(f"Bins thresholds saved to {os.path.join(outputdir, f'bins_jaccard_thresholds.csv')}\n")

    ### read bins thresholds df
    bins_thresholds_df = pd.read_csv(os.path.join(outputdir, f"bins_jaccard_thresholds.csv"), index_col=None)
    print(bins_thresholds_df.head())

    ### remove bins that cannot be trusted - setting min sample size based on chatGPT suggestion
    bins_thresholds_df = bins_thresholds_df[bins_thresholds_df['sampled_comparisons']> 150]

    # Convert bin_range to string for categorical plotting
    bins_thresholds_df['bin_range_str'] = bins_thresholds_df['bin_range'].astype(str)
    bins_thresholds_df.sort_values(by='bin_range', inplace=True)

    ### plot bin max to thresholds
    plt.figure(figsize=(max(8, 0.5 * len(bins_thresholds_df)), 5))  # dynamic width
    x = np.arange(len(bins_thresholds_df))
    plt.scatter(x, bins_thresholds_df['empirical_t'], marker='o', label='Empirical threshold')
    plt.xlabel('Bin')
    ### rotate x labels
    plt.xticks(x, bins_thresholds_df['taxname'], rotation=45, ha='right')
    plt.ylabel('Jaccard threshold')
    plt.title('Jaccard Thresholds vs. Bin')
    plt.legend()
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave room at bottom
    plt.subplots_adjust(bottom=0.35)         # guarantee label visibility
    plt.savefig(os.path.join(outputdir, "bin_vs_jaccard_thresholds.png"),dpi=300, bbox_inches='tight')
    plt.close()

def parse_args():
    default_profiler = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "profiler.py",
    )
    parser = argparse.ArgumentParser(
        description="Generate jaccard distributions for theshold testing. Designed for HogProf levels."
    )
    parser.add_argument(
        "--hogprof_folder",
        required=True,
        help="Folder where the result of HogProf LSHbuilder is stored."
    )
    parser.add_argument(
            "--expected_root",
            required=False,
            default="Metazoa",
            help="Expected root taxon (to replace 'internal_0')."
        )
    parser.add_argument(
        "--output",
        required=True,
        help="Output folder."
    )
    parser.add_argument(
        "--profiler-path",
        default=default_profiler,
        help=(
            "Path to profiler.py. "
            "Defaults to profiler.py in the same directory as this script."
        ),
    )
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    main(
        args.hogprof_folder,
        args.output,
        args.profiler_path,
        args.expected_root
        )

### print end time
end_time_secs = time()
end_time = datetime.fromtimestamp(end_time_secs).strftime('%Y-%m-%d %H:%M:%S')
print(f"Script finished at: {end_time}\n")
execution_time = end_time_secs - start_time_secs
### print execution time but reformat it to be human readable
print(f"Total execution time: {round(execution_time, 2)} seconds\n")
print("\nDONE!\n")