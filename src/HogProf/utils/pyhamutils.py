import pyham
import xml.etree.cElementTree as ET
import ete3
import pickle
import traceback
import sys

def get_orthoxml_oma(fam, db_obj):
    orthoxml = db_obj.get_orthoxml(fam).decode()    
    return orthoxml

def get_orthoxml_tar(fam, tar):
    f = tar.extractfile(fam)
    if f is not None:
        return f.read()
    else:
        raise Exception( member + ' : not found in tarfile ')
    return orthoxml

def get_species_from_orthoxml(orthoxml):
    NCBI_taxid2name = {}
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            NCBI_taxid2name[child.attrib['NCBITaxId']] = child.attrib['name']
    return NCBI_taxid2name

def switch_name_ncbi_id(orthoxml , mapdict = None  ):
    #swap ncbi taxid for species name to avoid ambiguity
    #mapdict should be a mapping from species name to taxid if the info isnt in the orthoxmls
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            child.attrib['name'] = child.attrib['NCBITaxId']
        elif mapdict:
            child.attrib['name'] = mapdict[child.attrib['name']]
        
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml

def reformat_treenames( tree , mapdict = None  ):   
    #tree is an ete3 tree instance
    #replace ( ) - / . and spaces with underscores
    #iterate over all nodes
    for node in tree.traverse():
        if mapdict:
            node.name = mapdict[node.name]
        else:
            node.name = node.name.replace('(', '').replace(')', '').replace('-', '_').replace('/', '_')
    return tree

def reformat_names_orthoxml(orthoxml , mapdict = None  ):
    #replace ( ) - / . and spaces with underscores
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            child.attrib['name'] = child.attrib['name'].replace('(', '').replace(')', '').replace('-', '_').replace('/', '_')
        elif mapdict:
            child.attrib['name'] = mapdict[child.attrib['name']]
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml

def create_nodemapping(tree):
    #create a mapping from node name to node
    nodemapping = {}
    for i,node in enumerate(tree.traverse()):
        nodemapping[node.name] = str(i)
    
    #assert number of nodes is equal to number of unique names
    assert len([n for n in tree.traverse()]) == len(set(nodemapping.values()))

    return nodemapping

def tree2numerical(tree):
    mapper = create_nodemapping(tree)
    for i,node in enumerate(tree.traverse()):
        node.name = mapper[node.name]
    return tree , mapper

def orthoxml2numerical(orthoxml , mapper):
    #print(orthoxml)
    try:
        root = ET.fromstring(orthoxml)
    except Exception as e:
        with open( orthoxml , 'r') as f:
            orthoxml = f.read()
        root = ET.fromstring(orthoxml)
    for child in root:
        if 'name' in child.attrib:
            child.attrib['name'] = mapper[child.attrib['name']]
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml 

def get_ham_treemap_from_row(row, tree , levels = None , swap_ids = True , orthoXML_as_string = True , use_phyloxml = False , use_internal_name = True ,reformat_names= True, orthomapper = None,
                             limit_species = 10, limit_events = 0, dataset_nodes = None):  
    fam, orthoxml = row
    format = 'newick_string'
    if use_phyloxml:
        format = 'phyloxml'
    if orthoxml:
        if swap_ids == True and orthoXML_as_string == True:
            orthoxml = switch_name_ncbi_id(orthoxml)
            quoted = False
        elif reformat_names == True:
            orthoxml = orthoxml2numerical(orthoxml , orthomapper)
            orthoXML_as_string = True
            quoted = False
        else:
            quoted = True
        try:
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string )            
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0]) 
            if dataset_nodes is not None:
                if tp.treemap.name not in dataset_nodes:
                    return None
            #check for losses / events and n leaves 
            return tp.treemap
        except Exception as e:
            # Capture the exception and format the traceback
            full_error_message = str(e)
            if 'TypeError: species name ' in full_error_message and 'maps to an ancestral name, not a leaf' in full_error_message:
                print('error' , full_error_message)
                #species name from bullshit error
                #TypeError: species name '3515' maps to an ancestral name, not a leaf of the taxono
                species = full_error_message.split('species name ')[1].split(' ')[0].replace('\'','')
                print( 'trim tree'+species)
                tree = ete3.Tree(tree , format = 1)
                #select all nodes with name = species
                nodes = tree.search_nodes(name = species)
                #get the first node
                node = nodes[0]
                for c in node.children:
                    #delete all children
                    c.delete()
                #rerun with trimmed tree
                ham_obj = pyham.Ham(tree.write(format=1), orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string )
            else:
                print('error' , full_error_message)
                return None



def get_subhog_ham_treemaps_from_row(row, tree , levels = None , swap_ids = True , orthoXML_as_string = True , use_phyloxml = False , use_internal_name = True ,reformat_names= True, orthomapper = None,
                                     limit_species = 10, limit_events = 0, dataset_nodes = None, hogid_for_all = None):  
    fam, orthoxml = row
    format = 'newick_string'
    def check_limits(treenode, limit_species, limit_events, subhogname):
                ###removed because already covered in generates_dataframes!!!!!!!
                ### leaves counting failed e.g. in HOG:E0712183.1e, counts more than there is
                #print(dir(treenode))
                leaves_num = sum(1 for node in treenode.traverse() if node.is_leaf())
                if leaves_num < limit_species:
                    #print(treenode.name,subhogname)
                    return False

                total_dupl = 0
                total_loss = 0
                for node in treenode.traverse():
                    try:
                        total_dupl += node.dupl
                    except:
                        total_dupl += 0
                    try:
                        total_loss += node.lost
                    except:
                        total_loss += 0
                if total_dupl > limit_events or total_loss > limit_events:
                    return True
                #print(treenode.name,subhogname, total_dupl, total_loss)
                return False
    if use_phyloxml:
        format = 'phyloxml'
    if orthoxml:
        if swap_ids == True and orthoXML_as_string == True and reformat_names == False:
            orthoxml = switch_name_ncbi_id(orthoxml)
            quoted = False
        elif reformat_names == True:
            orthoxml = orthoxml2numerical(orthoxml , orthomapper)
            orthoXML_as_string = True
            quoted = False
        else:
            quoted = True
        try:
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string ) 
            #print(dir(ham_obj)) 
            ### Create tree profile for the top-level HOG
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0]) 
            ### save root name
            #rootname = tp.treemap.name + '_0'
            
            ### get all subhogs
            subhogs  = tp.hog.get_all_descendant_hogs()
            hogid_for_all = subhogs[0].hog_id
            rootname = subhogs[0].genome.name + '_' + str(hogid_for_all) + '_0'
            ### try to get the HOG id to use as part of the subhog name
            #print(dir(subhogs[0]))
            #try:  
            if hogid_for_all is None:
                hogid_for_all = fam
            #hogs = { subhog.genome.name +'_' + str(subhog.hog_id) + '_' + str(i + 1):  ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }
            hogs = { subhog.genome.name +'_' + str(subhog.hog_id) + '_' + str(i + 1):  ham_obj.create_tree_profile(hog=subhog) for i,subhog in enumerate(subhogs) }
            ### manually add roothog cause apparently we are not including it
            #print(f'Subhogs: {len(hogs)}')
            hogs[rootname] = tp
            print(f'RootHOG: {rootname}')
            print(f'Subhogs total: {len(hogs)}')
            ### filter out the small ones and turn into treemaps
            hogs = {subhogname: hogs[subhogname].treemap for subhogname in hogs if len(hogs[subhogname].hog.get_all_descendant_genes()) > limit_species}
            print(f'Subhogs large enough: {len(hogs)}')
            ### it wont be possible in some cases, so just use the genome name (taxnode )
            #except Exception as e:     
            #    hogs = { subhog.genome.name + '_' + str(i):  ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }
            #print(hogs)

            ### If dataset_nodes are specified, avoid calculating unnecessary subhogs
            if dataset_nodes is not None:
                #print('fam', fam)
                #print('before',len(hogs.keys()))
                print([subhogname for subhogname in hogs])
                hogs = {subhogname: hogs[subhogname] for subhogname in hogs if subhogname.split('_')[0] in dataset_nodes}
                #print('after',len(hogs.keys()))
                if len(hogs) == 0:
                    print('no suitable subhogs')
                    return {}
            #'''
            ### first check rootHOG to see if there will be at least one hog returned
            ### if dataset_nodes is specified, this step cannot be done
            if dataset_nodes is None and not check_limits(hogs[rootname], limit_species, limit_events, 'root'):
                print('no suitable rootHOG')
                return {}

            ### then check subhogs and remove the ones that do not meet the limits
            hogs = {subhogname: hogs[subhogname] for subhogname in hogs if check_limits(hogs[subhogname], limit_species, limit_events, subhogname)}
            if len(hogs) == 0:
                print('no suitable subhogs')
            #print(hogs)
            #print(f'Subhogs after filtering: {len(hogs)}')
            #'''
            return hogs
        
        #'''
        except Exception as e:
            # Capture the exception and format the traceback
            full_error_message = str(e)
            if 'TypeError: species name ' in full_error_message and 'maps to an ancestral name, not a leaf' in full_error_message:
                print('error' , full_error_message, file=sys.stderr)
                #species name from bullshit error
                #TypeError: species name '3515' maps to an ancestral name, not a leaf of the taxono
                species = full_error_message.split('species name ')[1].split(' ')[0].replace('\'','')
                print( 'trim tree'+species)
                tree = ete3.Tree(tree , format = 1)
                #select all nodes with name = species
                nodes = tree.search_nodes(name = species)
                #get the first node
                node = nodes[0]
                for c in node.children:
                    #delete all children
                    c.delete()
                #rerun with trimmed tree    
                ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string ) 
                #print(dir(ham_obj)) 
                ### Create tree profile for the top-level HOG
                tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0]) 
                ### save root name
                rootname = tp.treemap.name + '_0'
                ### get all subhogs
                subhogs  = tp.hog.get_all_descendant_hogs() 
                if hogid_for_all is None:
                    hogid_for_all = fam      
                hogs = { subhog.genome.name +'_' + str(hogid_for_all)+ '_' + str(i):  ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }

                ### If dataset_nodes are specified, avoid calculating unnecessary subhogs
                if dataset_nodes is not None:
                    #print('fam', fam)
                    #print('before',len(hogs.keys()))
                    hogs = {subhogname: hogs[subhogname] for subhogname in hogs if subhogname.split('_')[0] in dataset_nodes}
                    #print('after',len(hogs.keys()))
                    if len(hogs) == 0:
                        return {}

                ### first check rootHOG to see if there will be at least one hog returned
                ### if dataset_nodes is specified, this step cannot be done
                if dataset_nodes is None and not check_limits(hogs[rootname], limit_species, limit_events, 'root'):
                    return {}

                ### then check subhogs and remove the ones that do not meet the limits
                hogs = {subhogname: hogs[subhogname] for subhogname in hogs if check_limits(hogs[subhogname], limit_species, limit_events, subhogname)}

                return hogs            
            else:
                print('error' , full_error_message, file=sys.stderr)
        return {}#None
        #'''


def yield_families(h5file, start_fam):
    """
    Given a h5file containing OMA server, returns an iterator over the families
    (not sure if still in use)
    :param h5file: omafile
    :param start_fam: fam to start on
    :return: fam number
    """
    for row in h5file.root.OrthoXML.Index:
        if row[0] > start_fam:
            yield row[0]
def get_one_family(i, h5file):
    '''
    get one family from database
    Args:
        i : family number
        h5file : OMA server file
    Return :
        family
    Not sure if still in use
    '''
    return h5file.root.OrthoXML.Index[i][0]
