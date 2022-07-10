"""
Purpose: To create functions that can 
compute statistics over a graph object
(whether it be a full neuron graph or a subgraph)

"""

import pandas_utils as pu
import networkx_utils as xu
import neuron_nx_utils as nxu
import pandas as pd

node_identifier = "u"

def summary_statistic_over_attributes(
    attributes,
    summary_statistic = "mean",
    summary_statisic_args = None,
    weight = None,
    G=None,
    node_df=None,
    verbose = False,
    return_df = False,
    
    ):
    if node_df is None:
        node_df = xu.node_df(G)
        
    node_df = node_df.query(f"{node_identifier} != 'S0'")
    
    return pu.summary_statistic_over_columns(
        columns = attributes,
        df = node_df,
        summary_statistic = summary_statistic,
        summary_statisic_args = summary_statisic_args,
        verbose = verbose,
        return_df = return_df,
        weight=weight,
    )

def n_branches(G):
    len(nxu.limb_branch_nodes(G))
    
def starting_coordinate(
    G,
    attribute = "endpoint_upstream"):
    """
    Purpose: To find the startng centroid of a graph

    Pseudocode: 
    1) Get the most upstream node
    """

    node = nxu.most_upstream_node(G)
    return G.nodes[node][attribute]

import time
def summary_statistic_over_dynamic_attribute(
    G,
    node,
    attribute,
    prefix = None,
    category_columns = None,
    attribute_summary_dicts = (dict(columns="volume",
                                   summary_statistic = "sum",
                                    summary_statisic_args = None),),
    verbose = False,
    default_value = None,
    debug_time = False,
    ):

    """
    Purpose: To get the dynamic attributes
    dictionary for a node in the graph
    """
    if default_value is None:
        default_value = {}
    
    if attribute_summary_dicts is not None:
        attribute_summary_dicts = list(attribute_summary_dicts)
    
    if debug_time:
        st = time.time()
        
    if prefix is None:
        prefix = attribute.replace("_data","")

    node_dict = G.nodes[node].copy()
    
    if len(node_dict[attribute]) == 0:
        return default_value
    
    
    df = pd.DataFrame.from_records(node_dict[attribute])
    
    if debug_time:
        print(f"Time for attribute dataframe = {time.time() - st}")
        st = time.time()

    return pu.summary_statistics_over_columns_by_category(
        df,
        prefix = prefix,
        category_columns = category_columns,
        attribute_summary_dicts = attribute_summary_dicts,
        add_counts_summary = True,
        verbose = verbose,
        special_count_name = True,
        debug_time = debug_time
    )

def summary_statistic_over_synapses(
    G,
    node,
    verbose = False,
    **kwargs):
    
    return nxst.summary_statistic_over_dynamic_attribute(
    G,
    node = node,
    attribute = "synapse_data",
    category_columns = ["head_neck_shaft","syn_type"],
    verbose = verbose,
    **kwargs
    )

def summary_statistic_over_spines(
    G,
    node,
    verbose = False,
    **kwargs):
    
    return nxst.summary_statistic_over_dynamic_attribute(
    G,
    node = node,
    attribute = "spine_data",
    category_columns = None,#["head_neck_shaft","syn_type"],
    verbose = verbose,
    **kwargs
    )

import time
def add_summary_statistic_over_dynamic_attributes_to_G(
    G,
    attributes = ("synapses","spines"),
    verbose = False,
    ):
    
    
    for node in nxu.limb_branch_nodes(G):
        st = time.time()
        #nxst = reload(nxst)
        curr_dict = {}
        for s in attributes:
            curr_dict.update(getattr(nxst,f"summary_statistic_over_{s}")(G,node))
            
        G.nodes[node].update(curr_dict)                         
                            
        if verbose:
            print(f"Node {node}: {time.time() - st}")



import neuron_nx_stats as nxst