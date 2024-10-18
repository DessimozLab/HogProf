import pyham
import xml.etree.cElementTree as ET
import ete3
import pickle
import traceback
import random
import os
import re
from io import StringIO


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

'''Added on 17.10.2024 by Athina to try to adjust the orthoxml to the tree with correct format'''
def remove_namespace(xml_string):
    # Parse the XML while preserving its structure
    it = ET.iterparse(StringIO(xml_string))
    for _, el in it:
        # Remove the namespace from each tag
        if '}' in el.tag:
            el.tag = el.tag.split('}', 1)[1]  # Remove the namespace
        # Remove namespace declarations
        el.attrib = {key.split('}', 1)[-1]: value for key, value in el.attrib.items()}
    # Convert back to string (preserving formatting)
    return ET.tostring(it.root, encoding='unicode')


def adjust_orthoxml_to_tree(xml_string, newick_str):
    #print('adjusting orthoxml to tree')
    #print('xml_string before change',xml_string)

    # Helper function to get tree nodes
    def get_tree_nodes(newick_str):
        tree = ete3.Tree(newick_str, format=1)
        return [node.name for node in tree.traverse()]
    ### Helper function to get the taxon ID of the target group
    def find_outermost_ortholog_with_generef(root):
        # Traverse all orthologGroup elements
        for ortholog_group in root.findall(".//orthologGroup"):
            # Check if the orthologGroup contains a geneRef directly
            gene_refs = ortholog_group.findall(".//geneRef")
            if gene_refs:
                return ortholog_group  # Return this orthologGroup if it contains a geneRef directly

            # Check if the orthologGroup contains a geneRef inside a paralogGroup
            paralog_groups = ortholog_group.findall(".//paralogGroup")
            for paralog_group in paralog_groups:
                gene_refs_in_paralog = paralog_group.findall(".//geneRef")
                if gene_refs_in_paralog:
                    return ortholog_group  # Return the outermost orthologGroup containing a geneRef

        return None  # Return None if no orthologGroup contains a geneRef
    
    ### Helpefr function to move HOG id to the outermost orthologGroup
    def move_hog_id(xml_string, target_taxon_id):
        root = ET.fromstring(xml_string)

        # Find the orthologGroup that has an id (we assume there's only one)
        source_group = root.find(".//orthologGroup[@id]")
        
        if source_group is not None:
            # Get the HOG ID (we don't know its exact value, so we extract it from the source group)
            hog_id = source_group.attrib.get("id")
            
            # Remove the HOG ID from the source group
            source_group.attrib.pop('id', None)  # Safely remove the id attribute

            # Find the target orthologGroup by taxonId
            target_group = root.find(f".//orthologGroup[@taxonId='{target_taxon_id}']")
            
            if target_group is not None:
                # Assign the HOG ID to the target group
                target_group.set("id", hog_id)

        # Return the updated XML as a string
        updated_xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
        return updated_xml_string

    ### Helper function to remove HOG ID from an orthologGroup
    def remove_id_from_orthologgroup(root):
    # Find the orthologGroup with the id attribute
        ortholog_group = root.find(".//{http://orthoXML.org/2011/}orthologGroup[@id]")

        if ortholog_group is not None:
            # Remove the id attribute
            ortholog_group.attrib.pop('id', None)  # Safely remove the id attribute

        # Return the modified XML as a string
        updated_xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
        return updated_xml_string
    
    tree_nodes = get_tree_nodes(newick_str)
    root = ET.fromstring(xml_string)

    valid_species_names = {node for node in tree_nodes}
    gene_id_set = set()
    species_found_count = 0
    
    # Collect gene IDs from all species
    for species in root.findall("{http://orthoXML.org/2011/}species"):
        name = species.get("name")
        if name in valid_species_names:
            genes = species.findall(".//{http://orthoXML.org/2011/}gene")
            for gene in genes:
                gene_id_set.add(gene.get("id"))
            species_found_count += 1
        else:
            root.remove(species)
    #print('species found count', species_found_count)
    
    if species_found_count == 0:
        return None
    
    ### testing non-empty only  -- ATHINA CHECKPOINT - 18.10.2024. I am using a mask for Toxicofera
    ### and only one HOG so that lshbuilder runs fast. I am trying to adjust the orthoxml string to 
    ### the changed tree so that the pyHam object can be created correctly. So far it does not work.
    ### Also to consider: even if it does work, is there a chance that it will break lshbuilder later
    ### at the hashing step? if things dont have the same length? keep it in mind.
    #'''
    test_xml_file = '/work/FAC/FBM/DBC/cdessim2/default/agavriilidou/venom_project/2a_hogprof_testing/orthoxml_fixed_manual.xml'
    tree = ET.parse(test_xml_file)
    root = tree.getroot()
    updated_xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
    return updated_xml_string
    #'''
    
    parent_map = {child: parent for parent in root.iter() for child in parent}

    ### move the HOG ID and remove geneRefs when irelevant
    groups_element = root.find("{http://orthoXML.org/2011/}groups")
    if groups_element is not None:
        ortholog_groups_to_remove = []
        firstgroup=True
        checkiffirstwithgeneref = True
        for ortholog_group in groups_element.findall(".//{http://orthoXML.org/2011/}orthologGroup"):
            if firstgroup:
                print('first ortholog group', ortholog_group.find('.//{http://orthoXML.org/2011/}property').get('value'))
                hogid = ortholog_group.get('id')
                ortholog_group.attrib.pop('id', None)
                firstgroup=False
                #print(ortholog_group.get('id'))
            ### exact syntax for lookin in this group only - no children
            gene_refs = ortholog_group.findall("./{http://orthoXML.org/2011/}geneRef")
            
            # Remove invalid geneRefs
            for gene_ref in gene_refs:
                gene_id = gene_ref.get("id")
                #print('gene_id',gene_id)
                if gene_id not in gene_id_set:
                    parent = parent_map[gene_ref]
                    #print(f"Removing geneRef {gene_id} from orthologGroup {ortholog_group.find('.//{http://orthoXML.org/2011/}property').get('value')}")
                    parent.remove(gene_ref)
                else:
                    if checkiffirstwithgeneref:
                        ortholog_group.set('id', hogid)
                        checkiffirstwithgeneref = False
                        print('highest populated group', ortholog_group.find('.//{http://orthoXML.org/2011/}property').get('value'))
            #print('next')



    
    # Find the orthologGroup that has an id (we assume there's only one)
    #source_group = root.find(".//orthologGroup[@id]")
    #target_group = find_outermost_ortholog_with_generef(root)
    #if target_group is None:
    #    print('No target group found')
    #    print(ET.tostring(root, encoding='utf-8').decode('utf-8'))
    #target_group_taxonid = target_group.get("taxonId")
    #move_hog_id(xml_string, target_group_taxonid)



            # Collapse ortholog groups cautiously, keeping top-level groups intact
            #collapse_ortholog_groups(ortholog_group, gene_id_set)

        #for ortholog_group in ortholog_groups_to_remove:
        #    parent = parent_map.get(ortholog_group)
        #    if parent is not None and ortholog_group in parent:
        #        parent.remove(ortholog_group)

    filtered_xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
    filtered_xml_string = remove_namespace(filtered_xml_string)
    
    return filtered_xml_string


def collapse_ortholog_groups(ortholog_group, gene_id_set):
    # Track innermost group, deepest taxonId, and property
    current_group = ortholog_group
    last_valid_taxon_id = None
    last_valid_property = None

    while True:
        sub_groups = current_group.findall("{http://orthoXML.org/2011/}orthologGroup")
        print('sub_groups',sub_groups)
        if not sub_groups:
            break
        current_group = sub_groups[0]  # Collapse into the first sub-group
        
        # Update the deepest taxonId and property if they exist
        taxon_id = current_group.get("taxonId")
        if taxon_id:
            last_valid_taxon_id = taxon_id
        property_element = current_group.find("{http://orthoXML.org/2011/}property")
        if property_element is not None:
            last_valid_property = property_element

    # Set the outermost group's taxonId to the deepest valid taxonId
    if last_valid_taxon_id is not None:
        ortholog_group.set("taxonId", last_valid_taxon_id)

    # Remove all inner ortholog groups, only if no valid geneRefs exist
    for sub_group in ortholog_group.findall("{http://orthoXML.org/2011/}orthologGroup"):
        if not sub_group.findall(".//{http://orthoXML.org/2011/}geneRef"):
            ortholog_group.remove(sub_group)

    # Replace the property of the outer group with the deepest valid one
    if last_valid_property is not None:
        existing_property = ortholog_group.find("{http://orthoXML.org/2011/}property")
        if existing_property is not None:
            ortholog_group.remove(existing_property)
        ortholog_group.append(last_valid_property)














def get_ham_treemap_from_row(row, tree , levels = None , swap_ids = True , orthoXML_as_string = True , use_phyloxml = False , use_internal_name = True ,reformat_names= False, orthomapper = None ):  
    fam, orthoxml = row
    #print('fam',fam)
    format = 'newick_string'
    
    '''
     # Print all the arguments (both positional and keyword arguments)
    print("\nFunction arguments:")
    for arg_name, arg_value in locals().items():
        print(f"{arg_name}: {arg_value}")
    '''
    #print('tree',tree)
    
    '''
    treeete = ete3.Tree(tree, format=1)
    # Count the total number of nodes in the tree (including root)
    total_nodes = len([node for node in treeete.traverse()])
    print('Total number of nodes in the tree:', total_nodes)
    '''
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
        #print('orthoxml',orthoxml)
        #if fam == 706758:
        #    print('orthoxml before adjustment:', orthoxml)
        ### 14.10.2024: added update orthoxml to match the input tree
        orthoxml = adjust_orthoxml_to_tree(orthoxml, tree)
        ### check if the orthoxml is None
        if orthoxml is None:
            return None
        print('family', fam)
        if fam == 706758:
            print('orthoxml after adjustment:', orthoxml)
        #print(orthoxml)
        ### 17.10.2024: this was just for testing, to get the full error message
        #ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, 
        #                        orthoXML_as_string=orthoXML_as_string ) 
        try:
            #print('Creating HAM object')
            # return multiple treemaps corresponding to slices at different levels
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml" , tree_format = format  , use_internal_name=use_internal_name, 
                                orthoXML_as_string=orthoXML_as_string )   
            #print('Created HAM object')
            '''
            if not os.path.exists('ham_obj.html'):
                current_directory = os.getcwd()
                print(f"Current working directory: {current_directory}")
                print('creating ham obj at', 'ham_obj.html')
                print(ham_obj.get_list_top_level_hogs())
                ham_obj.create_iHam(hog=ham_obj.get_list_top_level_hogs()[0], outfile='ham_obj.html')
            #print('HAM leaves',len(ham_obj.taxonomy.leaves))
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0])  
            if not os.path.exists('tp.html'):
                print('creating tp at', 'tp.html')
                ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0],outfile='tp.html',as_html=True)
            '''
            if ham_obj.get_list_top_level_hogs() == []:
                #if fam == 706758:
                print(f'No top level hogs found for family {fam}.')
                print(orthoxml)
                return None
            else:
                print(f'Top level hogs found for family {fam}.')
            # Parse the tree_string back into an ete3 Tree object
            treeete = ete3.Tree(tree, format=1)
            # Count the total number of nodes in the tree (including root)
            total_nodes = len([node for node in treeete.traverse()])
            print('Total number of nodes in the tree:', total_nodes)        
            print('HAM obj',ham_obj)
            #print('orthoxml',orthoxml)
            
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
                tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0])  #### possibly where multi levels should be introduced!!!!!!
                
                return tp.treemap
            else:
                print('error' , full_error_message)
                print('Warning: no fix was possible')
                print(orthoxml)
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
