import pyham
import xml.etree.cElementTree as ET
import ete3
import os
import pickle
import traceback

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

def get_ham_treemap_from_row(row, tree , levels = None , swap_ids = True , orthoXML_as_string = True , use_phyloxml = False , use_internal_name = True ,reformat_names= False, orthomapper = None , fallback = None):
    fam, orthoxml = row
    format = 'newick_string'
    if use_phyloxml:
        format = 'phyloxml'
    
    if os.path.exists('fallback.nwk'):
        tree = ete3.Tree('fallback.nwk' , format = 1 ).write(format = 1)

    if fallback:
        tree = ete3.Tree(fallback , format = 1 ).write(format = 1)

    if orthoxml:
        if swap_ids == True and orthoXML_as_string == True:
            orthoxml = switch_name_ncbi_id(orthoxml)
            quoted = False
        elif reformat_names == True and orthoXML_as_string == True:
            orthoxml = orthoxml2numerical(orthoxml , orthomapper)
            quoted = False
        else:
            quoted = True
        try:
            # return multiple treemaps corresponding to slices at different levels
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, orthoXML_as_string=orthoXML_as_string )            
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0]) 
            #check for losses / events and n leaves 
            return tp.treemap
        except Exception as e:
            # Capture the exception and format the traceback
            full_error_message = str(e)
            if 'maps to an ancestral name, not a leaf' in full_error_message:
                #species name from bullshit error
                #TypeError: species name '3515' maps to an ancestral name, not a leaf of the taxono
                species = full_error_message.split('species name ')[1].split(' ')[0].replace('\'','')
                if swap_ids == False and reformat_names == False:
                    species = ' '.join( full_error_message.split('species name ')[1].split(' ')[0:2] ).replace('\'','')

                #print( 'trim tree : '+species)
                if reformat_names == True and orthoXML_as_string == True:
                    species = str(species)


                tree = ete3.Tree(tree , format = 1 )
                #select all nodes with name = species

                nodes = tree.search_nodes(name = species)

                #print( 'nodes' , nodes) 
                #print( 'children' , nodes[0].get_children())
                #get the first node
                node = nodes[0]
                #get parent
                parent = node.up

                #create polytomy with children and internal node
                for child in node.get_children():
                    child.detach()
                    parent.add_child(child)
                #remove node
                tree.write(     outfile = 'fallback.nwk' , format = 1)
                
               
                #make sure node names are correctly formatted   
                try:
                    #rerun with trimmed tree
                    return get_ham_treemap_from_row(row, tree.write(format = 1) , levels = levels , swap_ids = swap_ids , orthoXML_as_string = orthoXML_as_string , use_phyloxml = use_phyloxml ,
                                                     use_internal_name = use_internal_name ,reformat_names= reformat_names, orthomapper = orthomapper , fallback= 'fallback.nwk' )
                except Exception as e:
                    print('error' , full_error_message)
                    return None
            else:
                print('error' , full_error_message)
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
