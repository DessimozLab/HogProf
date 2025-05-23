import profiler
import argparse
import pandas as pd
import h5py
import os


def main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, allvsall = False, k= 20, fams_list = [], 
         hogs_list = [], oma_path=""):
    """
    Main function to extract hits from the profiler.
    """
    
    ### check how many things are present in the hashes file
    with h5py.File(hashes_h5, 'r') as f:
        def visitor(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"\n📦 Dataset: {name}")
                print(f"   - Shape: {obj.shape}")
                print(f"   - Dtype: {obj.dtype}")
                try:
                    data_preview = obj[()]
                    # Print a few rows depending on dimensionality
                    if data_preview.ndim == 1:
                        print("   - Preview:", data_preview[:3])
                    else:
                        print("   - Preview:\n", data_preview[:3])
                except Exception as e:
                    print(f"   - Preview error: {e}")
            elif isinstance(obj, h5py.Group):
                print(f"\n📁 Group: {name}")
        print(f"🔍 Scanning: {hashes_h5}")
        f.visititems(visitor)

    ### check that the tree exists:
    if os.path.exists(treepath) == False:
        print(f"Error! The tree file {treepath} does not exist!")
        exit()

    if oma_path == "":
        print("\nBuilding the profiler")
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
    else:
        print("\nBuilding the profiler with OMA file")
        p = profiler.Profiler(lshforestpath = lshforestpath, 
                            hashes_h5= hashes_h5, 
                            mat_path= mat_path, #outputfolder + 'profilersavingpath.csv' ,
                            oma =  "/home/agavriil/Documents/venom_project/2a_hogprof_testing/local_oma_test/OmaServer.h5", 
                            nsamples = 256 ,
                            mastertree = treepath,
                            slicesubhogs = True
                            )
        ### list to be filled with the results
        sorted_hogs_dfs_list = []
        if fams_list != []:
            print(f"Number of families: {len(fams_list)}")
            if len(hogs_list) > 0:
                print("Warning! You have provided a list of hogs and a list of families (or all vs all). " \
                "\nThe list of hogs will be ignored.")
            ### check if all the fams in the list are integers (works for orthoxmls)
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


    ### save the df
    totalsortedhogs = pd.concat(sorted_hogs_dfs_list)
    print(totalsortedhogs.head())
    totalsortedhogs.to_csv(outputfile, index=False)
    print(f"Created file: {outputfile}")
    

'''
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
        required=False,
        default="",
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
        help="List of families to query (format 1,2,3 as in fam2orthoxml).",
    )
    parser.add_argument(
        "--hogs_list",
        type=str,
        nargs="*",
        default=[],
        help="List of HOGs to query (format 0_0_HOG:E0712183_1394_1).",
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
'''

#'''
### debugging
#outputfolder = '/home/agavriil/Documents/venom_project/2a_hogprof_testing/sauria_hogprof/'
#inputfolder = '/home/agavriil/Documents/venom_project/hogprof_levels_scripts/test_data/'
#lshforestpath = outputfolder + "newlshforest.pkl"
#hashes_h5 = outputfolder + "hashes.h5"
#mat_path = outputfolder + "fam2orthoxml.csv"
#treepath = inputfolder + "Sauria_speciestree_edited.newick"
#outputfile = outputfolder + "extracted_hits.csv"
#main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, hogs_list=['0_0_HOG:E0712183_1394_1'])
#main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, allvsall=True)

outputfolder = '/home/agavriil/Documents/venom_project/2a_hogprof_testing/levels_oma_sauria_250523/'
lshforestpath = outputfolder + "newlshforest.pkl"
hashes_h5 = outputfolder + "hashes.h5"
treepath = outputfolder + "speciestree.nwk"
treepath = outputfolder + "reformatted_tree.nwk"
outputfile = outputfolder + "extracted_hits_fromhogs.csv"
mat_path = outputfolder + "fam2orthoxml.csv"
mat_path = outputfolder + "profilersavingpath.csv"
oma_path = "/home/agavriil/Documents/venom_project/2a_hogprof_testing/local_oma_test/OmaServer.h5"
#main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, hogs_list=['712183_3584_HOG:E0712183_0'], oma_path=oma_path)
outputfile = outputfolder + "extracted_hits_fromfams.csv"
#main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, fams_list=[712183], oma_path=oma_path)

### get all venom hogs and extract results
preliminary_venom_omamer_table = "/home/agavriil/Documents/venom_project/1_omamer_results/uniprot_venom_Metazoa_omamer.hogmap.xlsx"
hogs_df = pd.read_excel(preliminary_venom_omamer_table)
outputfile = outputfolder + "extracted_hits_fromfams_venom.csv"
### get hogids only if not nan
hogs_df = hogs_df[hogs_df['hogid'].notna()]
hogs_df['hogid'] = hogs_df['hogid'].apply(lambda x: str(x).split('.')[0])
hogs_df['hogid'] = hogs_df['hogid'].apply(lambda x: str(x).split(':E')[1])
### turn into list
hogs_list = list(hogs_df['hogid'].unique())
hogs_list = [int(i) for i in hogs_list]
### get results - for 182 fams (not all passed the filter in the lshbuilder) it took at least 35min
main(lshforestpath, hashes_h5, treepath, mat_path, outputfile, fams_list=hogs_list, oma_path=oma_path, k=20)
#'''