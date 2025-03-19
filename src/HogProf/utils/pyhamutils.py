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
                             limit_species = 10, limit_events = 0, dataset_nodes = None, taxmapper = {}):  
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
                                     limit_species = 10, limit_events = 0, dataset_nodes = None, hogid_for_all = None, taxmapper = {}):  
    fam, orthoxml = row
    format = 'newick_string'
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
        #### subfunctions here
        '''function to check if a treenode meets the limits of species and events'''
        def check_limits(treenode, limit_species, limit_events, subhogname):
            leaves_num = sum(1 for node in treenode.traverse() if node.is_leaf())
            if leaves_num < limit_species:
                return False
            total_dupl = sum(getattr(node, 'dupl', 0) for node in treenode.traverse())
            total_loss = sum(getattr(node, 'lost', 0) for node in treenode.traverse())
            return total_dupl > limit_events or total_loss > limit_events or (total_dupl + total_loss) > limit_events
        '''function to get the taxrange of a subhog'''    
        def get_subhog_taxrange(subhog):
            return taxmapper[subhog.genome.name]
        '''function to get all taxids from a treemap object'''
        def get_hog_taxids(treeprofile):
            toptaxid = treeprofile.hog.genome.name
            ids_list = treeprofile.treemap.search_nodes(name=toptaxid)[0]
            return [node.name for node in ids_list] 
        '''function to extract taxids from subhogids. Needs to be adjusted as subhogids change'''
        def get_taxid_from_subhogid(subhogid):
            ### assuming format HOG:E0707322.2i_8570_7_51
            #return subhogid.split("_")[-2]
            ### assuming format '1912_Protostomia_HOG:E0656112_0'
            return subhogid.split("_")[0] 
        '''main subfunction to return all subhogs''' 
        def return_hogs(tree, orthoxml, format, use_internal_name, orthoXML_as_string, fam, limit_species, limit_events, dataset_nodes, hogid_for_all):
            ### initialize the ham object
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string ) 
            #print(dir(ham_obj)) 
            ### Create tree profile for the top-level HOG
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0]) 
            ### if nothing in the tree profile belongs in the dataset, return empty dict (avoid calculating subhogs)
            if not any(taxid in dataset_nodes for taxid in get_hog_taxids(tp)):
                return {}
            ### save root name
            rootname = tp.treemap.name + '_0'
            ### get all subhogs
            subhogs  = tp.hog.get_all_descendant_hogs()
            hogid_for_all = hogid_for_all or fam
            #'''
            #hogs = { subhog.genome.name +'_' + str(hogid_for_all) + '_' + str(i):  ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }
            #'''
            ### working for non augmented orthoxmls but does not recover OMA browser subhog ID.
            hogs = { subhog.genome.name + '_' + get_subhog_taxrange(subhog) +'_' + str(hogid_for_all) + '_' + str(i):  
                    ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }
            ### working for (some) augmented orthoxmls. but we don't have those yet (planned for next OMA version)
            #hogs = { subhog.hog_id + '_' + subhog.genome.name + '_' + str(i):  ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }
            #print(hogs.keys())
            ### it wont be possible in some cases, so just use the genome name (taxnode )
            #except Exception as e:     
            #    hogs = { subhog.genome.name + '_' + str(i):  ham_obj.create_tree_profile(hog=subhog).treemap for i,subhog in enumerate(subhogs) }
            #print('hogs keys',hogs.keys())
            ### If dataset_nodes are specified, avoid calculating unnecessary subhogs
            if dataset_nodes is not None:
                hogs = {subhogname: hogs[subhogname] for subhogname in hogs if get_taxid_from_subhogid(subhogname) in dataset_nodes}
                if not hogs:
                    return {}
            #'''
            ### first check rootHOG to see if there will be at least one hog returned
            ### if dataset_nodes is specified, this step cannot be done
            if dataset_nodes is None and not check_limits(hogs[rootname], limit_species, limit_events, 'root'):
                return {}
            ### then check subhogs and remove the ones that do not meet the limits
            hogs = {subhogname: hogs[subhogname] for subhogname in hogs if check_limits(hogs[subhogname], limit_species, limit_events, subhogname)}
            #'''
            return hogs
        ### here try to execute the main function and catch the exception
        try:
            return return_hogs(tree, orthoxml, format, use_internal_name, orthoXML_as_string, fam, limit_species, limit_events, dataset_nodes, hogid_for_all)
        #'''
        except Exception as e:
            # Capture the exception and format the traceback
            full_error_message = str(e)
            if 'species name ' in full_error_message and 'maps to an ancestral name, not a leaf' in full_error_message:
                #print('error' , full_error_message, file=sys.stderr)
                #species name from bullshit error
                #TypeError: species name '3515' maps to an ancestral name, not a leaf of the taxono
                species = full_error_message.split('species name ')[1].split(' ')[0].replace('\'','')
                #print( 'trim tree'+species)
                tree_ete = ete3.Tree(tree , format = 1)
                #select all nodes with name = species
                nodes = tree_ete.search_nodes(name = species)
                if not nodes:
                    print(f"Node with name {node_name} not found in the tree.")
                    pass
                #get the first node
                node = nodes[0]
                print(node.children)
                #### WARNINGThis part is a tar pit. It is not deleting all children. But ete3 has a will of
                #### its own. Hours wasted trying to get it to do what it's supposed to: 3.
                for c in node.children[:]:
                    #delete all children -> does this delete only one child?
                    c.delete()
                '''
                ### for example, this causes infinite loops
                while len(node.children) > 0:
                    for c in node.children[:]:
                        #delete all children -> does this delete only one child?
                        c.delete()
                '''
                print(node.children)
                #rerun with trimmed tree  
                trimmed_tree = tree_ete.write(format=1)
                try: 
                    return return_hogs(trimmed_tree, orthoxml, format, use_internal_name, orthoXML_as_string, fam, limit_species, limit_events, dataset_nodes, hogid_for_all)
                except:
                    #pass
                    print('error' , full_error_message, file=sys.stderr)
                '''
                ham_obj = pyham.Ham(trimmed_tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string ) 
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
                    hogs = {subhogname: hogs[subhogname] for subhogname in hogs if get_taxid_from_subhogid(subhogname) in dataset_nodes}
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
                '''            
            else:
                #print('error' , full_error_message, file=sys.stderr)
                pass
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
