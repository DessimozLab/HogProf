import profiler
import argparse
import pandas as pd




def main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, allvsall = False, k= 20, fams_list = [], hogs_list = []):
    """
    Main function to extract hits from the profiler.
    """
    p = profiler.Profiler(lshforestpath = lshforestpath, 
                          hashes_h5= hashes_h5, 
                          #mat_path= outputfolder + 'fam2orthoxml.csv' ,
                          oma = False , 
                          nsamples = 256 ,
                          mastertree = treepath,
                          slicesubhogs = True
                          )
    ### list to be filled with the results
    sorted_hogs_dfs_list = []

    ### if we want all vs all results
    if allvsall:
        ### print warning
        if len(fams_list) > 0:
            print("Warning! You have provided a list of families and you are also running all vs all. " \
            "\nThe list will be ignored.")
        if len(hogs_list) > 0:
            print("Warning! You have provided a list of hogs and you are also running all vs all. " \
            "\nThe list will be ignored.")
        print("\nRunning all vs all")
        print("Warning! This may take a while!")
        ### get list of families from the fam2orthoxml
        fam2orthoxml_df = pd.read_csv(mat_path)
        print(fam2orthoxml_df.head())
        fams_list = list(set(fam2orthoxml_df['fam'].tolist()))
        print(f"Number of families: {len(fams_list)}")
    if fams_list != []:
        if len(hogs_list) > 0:
            print("Warning! You have provided a list of hogs and a list of families (or all vs all). " \
            "\nThe list of hogs will be ignored.")
        ### check if all the fams in the list are integers
        try:
            fams_list = [int(i) for i in fams_list]
        except:
            print("Warning! Not all the fams in the list are integers. " \
            "\nCheck the input or fam2orthoxml table.")
            exit()
        ### if all well:
        for i in fams_list:
            #i+=1
            #try:
            hogdict, sortedhogs = p.hog_query_sorted( hog_id= i , k = k )
            sorted_hogs_dfs_list.append(sortedhogs)
            #except:
            #    break
    elif hogs_list != []:
        ### check if all the hogs in the list are strings
        try:
            hogs_list = [str(i) for i in hogs_list]
        except:
            print("Warning! Not all the hogs in the list are strings. " \
            "\nCheck the input.")
            exit()
        ### if all well:
        for i in hogs_list:
            hogdict, sortedhogs = p.hog_query_sorted( hog_id= i , k = k )
            sorted_hogs_dfs_list.append(sortedhogs)

    #print()
    #if sortedhogs['hit_subhogid'].str.contains('E0712183').any(): ## for local was 0_0_0, for curnagl was 0_4_0
    #    print('got hit!\n')
    #else:
    #    print('Warning! Did not find itself!\n')
    #'''

    ### save the df
    totalsortedhogs = pd.concat(sorted_hogs_dfs_list)
    print(totalsortedhogs.head())
    totalsortedhogs.to_csv(outputfile, index=False)
    print(f"Created file: {outputfile}")
    

#'''
def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="Extract hits from the profiler. Queries can be families or"
    " specific HOGs.")
    parser.add_argument(
        "--lshforestpath",
        type=str,
        required=True,
        help="Path to the LSH forest pickle file (newlshforest.pkl).",
    )
    parser.add_argument(
        "--hashes_h5",
        type=str,
        required=True,
        help="Path to the hashes HDF5 file (hashes.h5).",
    )
    parser.add_argument(
        "--treepath",
        type=str,
        required=True,
        help="Path to the species tree file (newick format).",
    )
    parser.add_argument(
        "--mat_path",
        type=str,
        required=True,
        help="Path to the fam2orthoxml CSV file.",
    )
    parser.add_argument(
        "--outputfile",
        type=str,
        required=True,
        help="Output csv file for the extracted hits.",
    )
    parser.add_argument(
        "--allvsall",
        action="store_true",
        help="Run all vs all mode (may take a long time).",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=20,
        help="Number of top hits to extract for each query.",
    )
    parser.add_argument(
        "--fams_list",
        type=str,
        nargs="*",
        default=[],
        help="List of families to query.",
    )
    parser.add_argument(
        "--hogs_list",
        type=str,
        nargs="*",
        default=[],
        help="List of HOGs to query.",
    )
    return parser.parse_args()

    
if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()

    # Run the main function with the parsed arguments
    main(
        lshforestpath=args.lshforestpath,
        hashes_h5=args.hashes_h5,
        treepath=args.treepath,
        mat_path=args.mat_path,
        outputfile=args.outputfile,
        allvsall=args.allvsall,
        k=args.k,
        fams_list=args.fams_list,
        hogs_list=args.hogs_list,
    )
#'''

'''
### debugging
outputfolder = '/home/agavriil/Documents/venom_project/2a_hogprof_testing/sauria_hogprof/'
inputfolder = '/home/agavriil/Documents/venom_project/hogprof_levels_scripts/test_data/'
lshforestpath = outputfolder + "newlshforest.pkl"
hashes_h5 = outputfolder + "hashes.h5"
mat_path = outputfolder + "fam2orthoxml.csv"
treepath = inputfolder + "Sauria_speciestree_edited.newick"
outputfile = outputfolder + "extracted_hits.csv"
#main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, hogs_list=['0_0_HOG:E0712183_1394_1'])
#main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, allvsall=True)
'''