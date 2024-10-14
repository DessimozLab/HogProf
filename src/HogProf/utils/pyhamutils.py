import pyham
import xml.etree.cElementTree as ET
import ete3
import pickle
import traceback
import random

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
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'name' in child.attrib:
            child.attrib['name'] = mapper[child.attrib['name']]
    orthoxml = ET.tostring(root, encoding='unicode', method='xml')
    return orthoxml 

def get_ham_treemap_from_row(row, tree , levels = None , swap_ids = True , orthoXML_as_string = True , use_phyloxml = False , use_internal_name = True ,reformat_names= False, orthomapper = None ):  
    fam, orthoxml = row
    format = 'newick_string'
    #print('tree',tree)
    #print('orthoxml',orthoxml)
    if use_phyloxml:
        format = 'phyloxml'
    if orthoxml:
        if swap_ids == True and orthoXML_as_string == True:
            orthoxml = switch_name_ncbi_id(orthoxml)
            quoted = False
        elif reformat_names == True and orthoXML_as_string == True:
            orthoxml = orthoxml2numerical(orthoxml , orthomapper)
            #print('reformatted orthoxml', orthoxml)
            quoted = False
        else:
            quoted = True
        #print(orthoxml)
        try:
            # return multiple treemaps corresponding to slices at different levels
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string )            
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0]) 
            #check for losses / events and n leaves 
            return tp.treemap
        except Exception as e:
            # Capture the exception and format the traceback
            full_error_message = str(e)
            if 'species name ' in full_error_message and 'maps to an ancestral name, not a leaf' in full_error_message:
                print('error' , full_error_message)
                #species name from bullshit error
                #TypeError: species name '3515' maps to an ancestral name, not a leaf of the taxono
                species = full_error_message.split('species name ')[1].split(' ')[0].replace('\'','')
                #print( 'trim tree '+species)
                tree = ete3.Tree(tree , format = 1)
                #select all nodes with name = species
                nodes = tree.search_nodes(name = species)
                #print(nodes)
                #get the first node
                node = nodes[0]
                #print(type(getattr(node, 'name')))
                #print("node.features",node.features)
                #for f in node.features:
                #    print(getattr(node, f))
                '''
                print('children to be removed',node.children)
                for c in node.children:
                    print(c)
                    #delete all children
                    c.delete()
                '''
                '''
                # Instead of deleting the children, remove them from the parent node
                for c in node.children[:]:  # Iterate over a copy of the children list
                    node.remove_child(c)  # Remove the child from the parent node without deleting the node
                # Search for the node again after deleting its children
                nodes_after_deletion = tree.search_nodes(name=species)
                # Check if the node still exists
                if nodes_after_deletion:
                    print(f"Node {species} is still in the tree.")
                else:
                    print(f"Node {species} has been deleted from the tree.")
                # Check if the node is now a leaf
                if node.is_leaf():
                    print(f"Node {node.name} is now a leaf (it has no children).")
                else:
                    print(f"Node {node.name} is not a leaf (it still has children).")
                '''

                ##### Adrian: dont remove node, get node (as is here) and add child to it with new unique id (e.g. negative)
                ##### and then update orthoxml to have this new id on the species tag name attribute (exact search and replace 
                # in orthoxml string) e.g. species name="3891"
                def get_new_random():
                    ### generate random integer
                    rando = random.randint(0, 1000000)
                    ### check if it is already in the orthoxml string
                    if str(rando) in orthoxml:
                        return get_new_random()
                    return rando
                ### create the new leaf
                newchildname = f'-{get_new_random()}'
                node.add_child(name = newchildname, dist = 0.0, support = node.support)
                ### update orthoxml
                orthoxml = orthoxml.replace(f'species name="{species}"', f'species name="{newchildname}"')


                '''
                while node.children:
                    #print("Current children:", node.children)
                    # Delete each child in the current list of children
                    for c in node.children[:]:  # [:] ensures you're iterating over a copy of the list
                        #print("Deleting child:", c)
                        node.remove_child(c)
                '''
                #print(node.children)
                #rerun with trimmed tree
                #print(orthoxml)
                ### revursive way:
                #print('recursion')
                #return get_ham_treemap_from_row(row, tree.write(format=1) , levels = levels , swap_ids = swap_ids , orthoXML_as_string = orthoXML_as_string , use_phyloxml = use_phyloxml , use_internal_name = use_internal_name ,reformat_names= reformat_names, orthomapper = orthomapper )
                ### original way:
                #'''
                #print('rerunning with updated tree in format:', format, 'use internal names:', use_internal_name)
                #print('orthoxml after trimming:', orthoxml)
                ham_obj = pyham.Ham(tree.write(format=1), orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string )
                #ham_obj = pyham.Ham(tree.write(format=1), orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string, species_resolve_mode="OMA")
                #'''
                tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0])
                return tp.treemap
            else:
                print('error' , full_error_message)
                print('Warning: no fix was possible')
                return None

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
