"""
Purpose: tools that will help process a neuron that
is represented as a networkx graph

"""
import networkx_utils as xu
import numpy as np
import pandas as pd

soma_node_name_global = "S0"
node_label = "u"

dynamic_attributes_default = (
 'spine_data',
 'synapse_data',
 'width_data',
  'skeleton_data')

def delete_attributes(
    G,
    attributes=None,
    verbose = False
    ):
    
    if attributes is None:
        attributes = list(nxu.dynamic_attributes_default)
    
    return xu.delete_node_attributes(
        G,
        attributes=attributes,
        verbose = verbose)

def skeletal_length_on_path(G,path):
    return np.sum([G.nodes[k]["skeletal_length"] for k in path])

def calculate_soma_distance_to_data_objs(
    G,
    verbose = False,
    data_to_update = ("synapse_data","spine_data","width_data"),                     
    ):
    """
    Purpose: To set the soma distance for all attributes

    Pseudocode: 
    Iterate through all of the nodes with L in name
    1) Find the path back to the soma
    2) Calculate the sum of skeletal length back to the soma
    for each attrbute
    3) Add skeletal length to the upstream dist to get the soma distance

    """
    for node in G.nodes():
        if "L" not in node:
            continue
        if verbose:
            print(f"--working on node {node}--")
        inbetween_node_path = xu.shortest_path(G,"S0",node)[1:-1]
        data_to_update = ("synapse_data","spine_data","width_data")

        soma_dist_on_path = nxu.skeletal_length_on_path(G,inbetween_node_path)

        if verbose:
            print(f"inbetween_node_path = {inbetween_node_path}")
            print(f"soma_dist_on_path = {soma_dist_on_path}")

        for dtype in data_to_update:
            for obj in G.nodes[node][dtype]:
                obj["soma_distance"] = obj["upstream_dist"] + soma_dist_on_path
                
                
def draw_tree(
    G,
    draw_type = "dot",#"twopi" (the circular), "circo" (makes square like)
    node_size=4000,
    font_size = 20,
    font_color = "white",
    figsize=(32,32),
    **kwargs):
    
    xu.draw_tree(
        G,
        draw_type = draw_type,#"twopi" (the circular), "circo" (makes square like)
        node_size=node_size,
        font_size = font_size,
        font_color = font_color,
        figsize=figsize,
        **kwargs)
    
    
skeletal_length_min_starter_branches_global = 700
def small_starter_branches(
    G,
    skeletal_length_min = None,
    soma_node_name = None,
    verbose = False,
    ):
    """
    Purpose: To identify starting nodes may want to remove

    Pseuodocde: 
    1) Look for those with small lengths
    2) Look for those that neighbor the soma
    """
    if skeletal_length_min is None:
        skeletal_length_min = skeletal_length_min_starter_branches_global

    if soma_node_name is None:
        soma_node_name= soma_node_name_global

    small_nodes = xu.nodes_from_node_query(G,f"(skeletal_length < {skeletal_length_min})")
    if verbose:
        print(f"small_nodes = {small_nodes}")
    small_soma_neighbors= []
    if len(small_nodes) > 0:
        soma_partners = xu.downstream_nodes(G,soma_node_name)
        if verbose:
            print(f"soma_partners= {soma_partners}")
        small_soma_neighbors = list(set(small_nodes).intersection(set(soma_partners)))
        
    if verbose:
        print(f"small_soma_neighbors= {small_soma_neighbors}")
        
    return small_soma_neighbors

import neuron_nx_utils as nxu
import copy

def remove_small_starter_branches(
    G,
    skeletal_length_min= None,
    inplace = False,
    verbose = True,
    **kwargs):
    """
    Purpose: To remove small starter branches from 
    the graph
    
    Ex: 
    nxu.remove_small_starter_branches(G,verbose = True)
    """

    if not inplace:
        G = copy.deepcopy(G)

    sm_st_branches = nxu.small_starter_branches(G,verbose = verbose,**kwargs)

    for s_b in sm_st_branches:
        if verbose:
            print(f"Removing {s_b}")
        xu.remove_edge_reattach_children_di(G,s_b,inplace = True)

    return G

import numpy as np
def soma_radius(
    G,
    stat = "mean",
    verbose = False
    ):
    """
    Purpose: To find the [stat] distance of the soma from the neighbors
    """
    soma_neighbors = xu.downstream_nodes(G,soma_node_name_global)
    s_center = G.nodes[soma_node_name_global]["mesh_center"]
    endpoints_dist = np.linalg.norm([G.nodes[k]["endpoint_upstream"] - s_center for k in soma_neighbors],axis = 1)
    stat_dist = getattr(np,stat)(endpoints_dist)
    if verbose:
        print(f"soma_neighbors = {soma_neighbors}")
        print(f"endpoints_dist = {endpoints_dist}")
        print(f"{stat} dist = {stat_dist}")
    return stat_dist
    
    
# --------------- For exporting graph as differnt file types -------------
import networkx as nx
def export_swc_dicts(
    G,
    use_skeletal_coordinates = True,
    soma_node_name = soma_node_name_global,
    default_skeleton_pt = "mesh_center",
    default_compartment = "undefined",
    center_on_soma= True,
    coordinate_divisor = 1000,
    width_divisor = 1000,
    verbose = False,
    return_df = False,
    ):
    """
    Purpose: To convert a Graph Neuron Object
    into a SWC file dict objects by creating dictionaries of all of the points

    Pseudocode: 
    In a depth first search manner of searching
    For each node: 
        1) get the skeleton points
        2) Get the compartment
        3) Get the parent idx

    """

    index = 0

    compartment_index_map = dict(
        dendrite = 0,
        soma = 1,
        axon = 2,
        basal = 3,
        apical = 4,
        apical_tuft = 4,
        oblique = 4,
        apical_shaft = 4,
        custom = 5,
        undefined = 6,
        glia = 7
    )

    center_point = G.nodes[soma_node_name]["mesh_center"]

    node_order = nx.dfs_preorder_nodes(G,source = "S0")
    node_to_index_map = dict()
    swc_dicts = []
    for n in node_order:
        if verbose:
            print(f"\n---Working on {n}")
        node_to_index_map[n] = []

        if use_skeletal_coordinates and n != soma_node_name:
            skeleton_pts = G.nodes[n]["skeleton_data"][1:]
        else:
            skeleton_pts = np.array([G.nodes[n][default_skeleton_pt]])

        if n == soma_node_name:
            parent_idx = -1
            width_points = [nxu.soma_radius(G,verbose = False)]
        else:
            parent_name = xu.parent_node(G,n)
            parent_idx = node_to_index_map[parent_name][-1]
            if use_skeletal_coordinates:
                width_points = [k["width"] for k in G.nodes[n]["width_data"]]
            else:
                width_points = [G.nodes[n]["width_new"]["no_spine_median_mesh_center"]]

        curr_compartment  = G.nodes[n]["compartment"]
        if curr_compartment is None:
            curr_compartment = default_compartment

        if center_on_soma:
            skeleton_pts = skeleton_pts - center_point

        if verbose:
            print(f"skeleton_pts = {skeleton_pts}")
            print(f"parent_idx= {parent_idx}")
            print(f"curr_compartment= {curr_compartment}")

        for i,coord in enumerate(skeleton_pts):
            index += 1

            if i == 0:
                curr_parent_idx = parent_idx
            else:
                curr_parent_idx = index - 1

            curr_dict = {
                "n":index,
                "type":compartment_index_map[curr_compartment],
                "x": coord[0]/coordinate_divisor,
                "y":coord[1]/coordinate_divisor,
                "z":coord[2]/coordinate_divisor,
                "radius": width_points[i]/width_divisor,
                "parent":curr_parent_idx,
            }
            node_to_index_map[n].append(index)

            swc_dicts.append(curr_dict) 
    
    if return_df:
        pd.DataFrame.from_records(swc_dicts)
    return swc_dicts

def export_swc_df(G,**kwargs):
    return pd.DataFrame.from_records(export_swc_dicts(G,return_df=True,**kwargs))


import file_utils as fileu
def export_swc_file(
    G,
    filename=None,
    filename_append = None,
    directory="./",
    filepath = None,
    verbose = False,
    **kwargs):
    """
    Purpose: Create a SWC file from 
    graph object
    """
    
    swc_dicts= nxu.export_swc_dicts(G,verbose = False)
    if filename is None:
        filename = f"seg_{G.graph['segment_id']}_split_{G.graph['split_index']}_nucleus_{G.graph['nucleus_id']}"

    if filename_append is not None:
        filename = f"{filename}_{filename_append}"
        
    if filename[-4:] != ".swc":
        filename += ".swc"
        
    return fileu.file_from_dicts(
        swc_dicts,
        filename = filename,
        directory = directory,
        filepath = filepath,
        seperation_character=" ",
        verbose = verbose
    )

# -------------- computing morphology metrics using morphopy -----
import morphopy_utils as mpu

def swc_df(G,flip_z = True,**kwargs):
    swc = nxu.export_swc_df(G)
    swc_df = mpu.swc_rotated_resampled(swc,resample=False,flip_z = flip_z,**kwargs)
    return swc_df

def morphometrics(
    G=None,
    swc = None,
    **kwargs
    ):
    
    if swc is None:
        swc = nxu.swc_df(G)
    
    
    #swc_df = mpu.swc_rotated_resampled(swc,resample=False,flip_z = True)
        
    return mpu.morphometrics(swc=swc,**kwargs)


# ------------ for querying neuron graphs ----
def axon_dendrite_nodes(
    G,
    compartment="axon",
    return_node_df=False,):
    
    node_df = xu.node_df_from_node_query(G,query=f"axon_compartment=='{compartment}'")
    if return_node_df:
        return node_df
    else:
        return node_df[node_label].to_numpy()
    
def axon_nodes(
    G,
    return_node_df=False):
    
    return nxu.axon_dendrite_nodes(
    G,
    compartment="axon",
    return_node_df=return_node_df,)

def dendrite_nodes(
    G,
    return_node_df=False):
    
    return nxu.axon_dendrite_nodes(
    G,
    compartment="dendrite",
    return_node_df=return_node_df,)

    
    
import neuron_nx_utils as nxu