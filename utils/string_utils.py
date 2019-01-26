import time
import pandas as pd
import requests
from pandas.io.json import json_normalize
pd.set_option('max_colwidth', 200)


# should get someone a h5file, a dbobj, an omaidobj, and an xrefObj

def get_interactions(all_string_1, all_string_2):
    interactions = {}
    # all_string_# are dictionaries: {genome:[strings]}
    for genome_1, string_id_1 in all_string_1.items():
        for genome_2, string_id_2 in all_string_2.items():
            if genome_1 == genome_2:
                string_interactors = get_string_interactors([string_id_1, string_id_2], genome_1)
                interactions[''.join(sorted([string_id_1, string_id_2]))] = string_interactors

    return interactions

"""
def get_string_interactors(list_of_string_ids, genome):
    '''Gets the interactors based on stringids'''
    final_df = pd.DataFrame()

    method = "interaction_partners"
    string_api_url = "https://string-db.org/api"
    output_format = "json"
    caller_identity = "clement.train@unil.ch"
    echo_query = "1"
    ncbitaxid = genome # TODO correct

    print('inside interactors ... ', list_of_string_ids)


    request_url = construct_string_request(string_api_url, output_format, list_of_string_ids,
                                               ncbitaxid, method, 2, echo_query, caller_identity)

    response = requests.get(request_url)
    df = json_normalize(response.json())

    final_df = final_df.append(df)
    if len(final_df) > 0:
        return final_df
    else:
        return []"""

"""
def construct_string_request(string_api_url, output_format, my_genes,
                             ncbitaxid, method, limit, echo_query, caller_identity):
    ## Construct the request
    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=" + "%0d".join(my_genes)
    request_url += "&" + "species=" + str(ncbitaxid)
    request_url += "&" + "limit=" + str(limit)
    request_url += "&" + "echo_query=" + echo_query
    request_url += "&" + "caller_identity=" + caller_identity

    return request_url
"""
