import pyham
import xml.etree.cElementTree as ET
import pickle
from utils import config_utils

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

def switch_name_ncbi_id(orthoxml , mapdict = None):
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



def get_ham_treemap_from_row(row, tree , level = None):
    fam, orthoxml = row
    orthoxml = switch_name_ncbi_id(orthoxml)
    try:
        if level is None:
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml", use_internal_name=True, orthoXML_as_string=True)
            tp = ham_obj.create_tree_profile(hog=ham_obj.get_list_top_level_hogs()[0])
            return tp.treemap
        else:
            ham_obj = pyham.Ham(tree, orthoxml, type_hog_file="orthoxml", use_internal_name=True, orthoXML_as_string=True)
            #return subHOGs at level
            slice = ham_obj.get_ancestral_genome_by_name(level)
            treeprofiles = [  ham_obj.create_tree_profile(hog=h) for h in ham_obj.get_list_top_level_hogs()[0].get_at_level(slice) ]

    except TypeError as err:
        print('Type error:', err)
        return None
    except AttributeError as err:
        print('Attribute error:', err)
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
