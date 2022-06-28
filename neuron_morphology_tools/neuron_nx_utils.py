"""
Purpose: tools that will help process a neuron that
is represented as a networkx graph

"""
import networkx_utils as xu
import numpy as np
import pandas as pd
import system_utils as su

soma_node_name_global = "S0"
node_label = "u"

auto_proof_filter_name_default = "auto_proof_filter"

dynamic_attributes_default = (
 'spine_data',
 'synapse_data',
 'width_data',
  'skeleton_data')
import copy


def delete_attributes(
    G,
    inplace = True,
    attributes=None,
    verbose = False
    ):
    
    if not inplace:
        G = copy.deepcopy(G)
    
    if attributes is None:
        attributes = list(nxu.dynamic_attributes_default)
    
    return xu.delete_node_attributes(
        G,
        attributes=attributes,
        verbose = verbose)


def name_from_G(G,append=None):
    graph_attr = xu.graph_attr_dict(G)
    segment_id = graph_attr.get("segment_id",None)
    split_index = graph_attr.get("split_index",None)
    nucleus_id = graph_attr.get("nucleus_id",None)
    
    name = f"seg_{segment_id}_{split_index}_nucleus_{nucleus_id}"
    if append is not None:
        name += f"_{append}"
    return name
    
def save_G(
    G,
    filepath=None,
    filename_append = None,
    delete_dynamic_attributes = True,
    verbose = False
    ):
    
    """
    Purpose: To save the graph
    nxu.save_G(G_ax,filename_append="axon_high_fid")
    """
    
    if delete_dynamic_attributes:
        G = nxu.delete_attributes(
            G,
            inplace = False,
            verbose = False)
    if filepath is None:
        filepath = f"./{nxu.name_from_G(G,append=filename_append)}"
    
    return su.compressed_pickle(G,filepath)

def load_G(filepath):
    return su.decompress_pickle(filepath)
    

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
        
    if nxu.soma_only_graph(G):
        return []

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

import copy
import networkx_utils as xu

def remove_node(
    G,
    node,
    inplace = False,
    verbose = False,
    maintain_skeleton_connectivity = True,
    remove_all_downstream_nodes = False,
    **kwargs
    ):
    """
    Purpose: To remove a node from the graph
    and to reconnect the downstream children

    Pseudocode: 
    1) Get any downstream nodes 
        a. (if none then just delete and return)
    2) For each downstream node:

    Attributes that need to be changed: 

    # if ask to alter skeleton
    endpoint_upstream --> parent endpoint_upstream
    skeleton_coordinates --> add parent endpoint_upstream

    #the soma information
    'soma_start_vec': array([-0.375219  , -0.91207943, -0.16529313]),
     'soma_start_angle': 24.21}
     
     
    Ez: 
    new_G = nxu.remove_node(
    G,
    node="L0_22",
    inplace = False,
    verbose = True,
    maintain_skeleton_connectivity = True,
        )

    nxu.soma_connected_nodes(new_G)
    new_G.nodes["L0_19"]["skeleton_data"],new_G.nodes["L0_20"]["skeleton_data"]
    
    
    Ex 2: Deleting all downstream nodes
    nodes = ["L1_10","L1_2","L1_8","L0_20"]
    G_del = nxu.remove_node(
        G,
        node=nodes,
        inplace=False,
        verbose = True,
        remove_all_downstream_nodes = True
    )

    """

    if not inplace:
        G = copy.deepcopy(G)

    nodes = nu.convert_to_array_like(node)
    for node in nodes:
        if verbose:
            print(f"--Working on removing node {node}")
        
        if node not in G:
            if verbose:
                print(f"node {node} wasn't in graph nodes so continuing")
            continue
            
        downstream_nodes = xu.downstream_nodes(G,node)

        if verbose:
            print(f"downstream_nodes = {downstream_nodes}")

        if len(downstream_nodes) > 0 and not remove_all_downstream_nodes:
            if node in nxu.soma_connected_nodes(G):
                soma_vals = ["soma_start_vec","soma_start_angle"]
                if verbose:
                    print(f"Adding {soma_vals} to the downstream nodes")
                for att in soma_vals:
                    for n in downstream_nodes:
                        G.nodes[n][att]  =  G.nodes[node][att]

            if maintain_skeleton_connectivity:
                endpoint_upstream = G.nodes[node]["endpoint_upstream"]
                width_upstream = G.nodes[node]["width_new"]["no_spine_median_mesh_center"]
                skeletal_length_upstream = G.nodes[node]["skeletal_length"]
                if verbose:
                    print(f"Adding endpoint_upstream {endpoint_upstream} to the downstream nodes skeleton")
                for n in downstream_nodes:
                    G.nodes[n]["endpoint_upstream"]  = endpoint_upstream
                    G.nodes[n]["skeleton_data"] = np.concatenate([[endpoint_upstream],G.nodes[n]["skeleton_data"]])
                    for i,(k) in enumerate(G.nodes[n]["width_data"]):
                        G.nodes[n]["width_data"][i]["upstream_dist"] = G.nodes[n]["width_data"][i]["upstream_dist"] + width_upstream
                    G.nodes[n]["width_data"] = [dict(upstream_dist=skeletal_length_upstream,
                                                     width = width_upstream)] + G.nodes[n]["width_data"]
                    
        if not remove_all_downstream_nodes:
            xu.remove_node_reattach_children_di(G,node,inplace = True)
        else:
            all_downstream_nodes = xu.all_downstream_nodes(G,node)
            total_nodes_to_delete = [node] 
            if len(all_downstream_nodes) > 0:
                total_nodes_to_delete += list(all_downstream_nodes)
                
            if verbose:
                print(f"Removing all downstream nodes along with node {node}: {total_nodes_to_delete}")
            G = xu.remove_nodes_from(G,total_nodes_to_delete)
            
    return G
    

def remove_small_starter_branches(
    G,
    skeletal_length_min= None,
    inplace = False,
    verbose = True,
    maintain_skeleton_connectivity = True,
    **kwargs):
    """
    Purpose: To remove small starter branches from 
    the graph
    
    Ex: 
    nxu.remove_small_starter_branches(G,verbose = True)
    
    Ex 2:
    new_G = nxu.remove_small_starter_branches(
        G,
        skeletal_length_min = 100000000000,
    )
    """

    if not inplace:
        G = copy.deepcopy(G)

    sm_st_branches = nxu.small_starter_branches(
        G,
        verbose = verbose,
        skeletal_length_min = skeletal_length_min,
        **kwargs)

#     for s_b in sm_st_branches:
# #         if verbose:
# #             print(f"Removing {s_b}")
    nxu.remove_node(
        G,
        sm_st_branches,
        verbose= verbose,
        inplace = True,
        maintain_skeleton_connectivity=maintain_skeleton_connectivity,
        **kwargs)

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

compartment_index_swc_map = dict(
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

def export_swc_dicts(
    G,
    use_skeletal_coordinates = True,
    soma_node_name = soma_node_name_global,
    default_skeleton_pt = "mesh_center",
    default_compartment = "basal",#"undefined",
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
            
#         if curr_compartment == "dendrite":
#             curr_compartment = default_compartment

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
                "type":compartment_index_swc_map[curr_compartment],
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
    apply_basal_dendrite_swap = True,
    **kwargs
    ):
    
    if swc is None:
        swc = nxu.swc_df(G)
        
    if apply_basal_dendrite_swap:
        # --- fixed so dendrite compartment will be changed to basal and recognized as dendrite ---
        import pandas_utils as pu
        def basal_rename(row):
            if row["type"] == compartment_index_swc_map["dendrite"]:
                return compartment_index_swc_map["basal"]
            else:
                return row["type"]

        swc["type"] = pu.new_column_from_row_function(swc,basal_rename).astype('int')
        #return swc
    
    #swc_df = mpu.swc_rotated_resampled(swc,resample=False,flip_z = True)
        
    return mpu.morphometrics(swc=swc,**kwargs)


# ------------ for querying neuron graphs ----
def axon_dendrite_nodes(
    G,
    compartment="axon",
    return_node_df=False,):
    
    node_df = xu.node_df(G)
    try:
        node_df = node_df.query(f"axon_compartment=='{compartment}'")
    except:
        pass 
    
    if return_node_df:
        return node_df
    else:
        return node_df[node_label].to_numpy()
    
def axon_nodes(
    G,
    return_node_df=False,
    **kwargs):
    
    return nxu.axon_dendrite_nodes(
    G,
    compartment="axon",
    return_node_df=return_node_df,**kwargs)

def dendrite_nodes(
    G,
    return_node_df=False):
    
    return nxu.axon_dendrite_nodes(
    G,
    compartment="dendrite",
    return_node_df=return_node_df,)

def axon_dendrite_subgraph(
    G,
    compartment,
    include_soma = True,
    verbose = False,
    ):
    
    compartment_nodes =  list(axon_dendrite_nodes(
    G,
    compartment=compartment,
    return_node_df=False,))
    
    if verbose:
        print(f"compartment_nodes = {compartment_nodes}")
        
    if include_soma:
        compartment_nodes.append(nxu.soma_node_name_global)
        
    return G.subgraph(compartment_nodes).copy()
        
def axon_subgraph(
    G,
    include_soma = True,
    verbose = False,
    ):
    return nxu.axon_dendrite_subgraph(
    G,
    compartment="axon",
    include_soma = include_soma,
    verbose = verbose,
    )

def dendrite_subgraph(
    G,
    include_soma = True,
    verbose = False,
    ):
    return nxu.axon_dendrite_subgraph(
    G,
    compartment="dendrite",
    include_soma = include_soma,
    verbose = verbose,
    )

import time
def plot(
    G,
    verbose = False,
    xlim = None,
    ylim=None,
    zlim = None,
    **kwargs):
    
    
    #st = time.time()
    swc_df = nxu.swc_df(G)
    ntree_obj = mpu.ntree_obj_from_swc(swc=swc_df)
    mpu.plot_ntree(
        ntree_obj,
        xlim = xlim,
        ylim=ylim,
        zlim = zlim,
        **kwargs)

# --------- mapping errors to nodes ----------
def limb_from_node_name(name):
    return int(name[1:name.find("_")])

def branch_from_node_name(name):
    return int(name[name.find("_")+1:])

def limb_branch_nodes(G):
    return [k for k in G.nodes() if "S" not in k]

def limb_branch_subgraph(G):
    return G.subgraph(nxu.limb_branch_nodes(G)).copy()

def all_limb_idxs_in_G(G):
    return np.sort(np.unique([nxu.limb_from_node_name(k)
                     for k in limb_branch_nodes(G)]))

def all_limb_graphs(G,return_idxs = False):
    limb_idxs = nxu.all_limb_idxs_in_G(G)
    return_graphs = [limb_graph(G,k) for k in limb_idxs]
    if return_idxs:
        return return_graphs,limb_idxs
    else:
        return return_graphs

def limb_graph(
    G,
    limb_idx=None,
    most_upstream_node = None,
    branches_idx = None,
    verbose = False,
    ):
    """
    Purpose: To return a graph of a certain limb
    
    Ex: nxu.limb_graph(G_ax,limb_idx = 3,verbose = True)
    """
    
    if limb_idx is not None:
        if verbose:
            print(f"Using the limb_idx method")
        if type(limb_idx) == str:
            limb_idx = int(limb_idx[1:])

        if limb_idx == -1:
            if verbose:
                print(f"Returning whole graph")
            return G

        subgraph_nodes = [k for k in nxu.limb_branch_nodes(G)
                         if nxu.limb_from_node_name(k) == limb_idx]
        
        if branches_idx is not None:
            subgraph_nodes = [k for k in subgraph_nodes
                             if nxu.branch_from_node_name(k) in branches_idx]

            if verbose:
                print(f"subgraph_nodes after branch restriction: {subgraph_nodes}")
    elif most_upstream_node is not None:
        if verbose:
            print(f"Using the most upstream method")
        subgraph_nodes = xu.all_downstream_nodes(
            G,
            most_upstream_node,
            include_self=True)
    else:
        raise Exception("")
    if verbose:
        print(f"subgraph_nodes after restriction: {subgraph_nodes}")
        
    return G.subgraph(subgraph_nodes).copy()


def node_match_by_dict(
    dict1,
    dict2,
    attributes = (
        "endpoint_upstream",
        #"endpoint_downstream",
        "skeleton_vector_upstream",
    ),
    verbose = False,
    ):
    """
    Purpose: To compare whether two
    nodes are the same or not
    
    Ex: 
    node_match_by_dict(
    G_ax.nodes["L0_0"],
    G_ax.nodes["L0_1"],
    verbose = True
    )
    """
    
    match_flag = True
    
    for a in attributes:
        val1 = dict1[a]
        val2 = dict2[a]
        
        if "array" in str(type(val1)):
            comp_result = np.array_equal(val1,val2)
        elif type(val1) == list:
            comp_result = np.all([k1==k2 for k1,k2 in zip(val1,val2)])
        else:
            comp_result  = val1==val2
            
        if comp_result == False:
            match_flag = False
            if verbose:
                print(f"Did not have equal {a}: {val1},{val2}")
            
    
    return match_flag


def node_map(
    G_source,
    G_target,
    verbose = False,
    ):
    """
    Purpose: To find the matches of all the nodes
    of one graph to another:

    Pseudocode: 
    For each node in the source graph
        Iterate through the target graph nodes
            Compare the certain properties, and if all match then assign as the correct mapping
            
            
    nxu.node_map(
        G_source = G_ax,
        G_target = G_proof,
        verbose = True
    )
    """
    source_to_target_map = dict()

    for n in nxu.limb_branch_nodes(G_source):

        limb_idx = nxu.limb_from_node_name(n)
        matching_nodes = []
        for nt in nxu.limb_branch_nodes(G_target):

            limb_idx_t = nxu.limb_from_node_name(nt)
            if limb_idx != limb_idx_t:
                continue

            match_result = nxu.node_match_by_dict(
                G_source.nodes[n],
                G_target.nodes[nt],
            )

            if match_result:
                matching_nodes.append(nt)

        if verbose:
            print(f"Node {n} --> {matching_nodes}")

        if len(matching_nodes) > 1:
            raise Exception("More than 1 match")

        if len(matching_nodes) == 0:
            source_to_target_map[n] = None
            continue

        match = matching_nodes[0]

        if nxu.limb_from_node_name(n) != nxu.limb_from_node_name(match):
            raise Exception("Not on same limb")

        source_to_target_map[n] = match

    #check to see if unique
    final_matches = list([k for k in source_to_target_map.values()
                          if k is not None])
    if len(final_matches) != len(np.unique(final_matches)):
        raise Exception("Some repeats")
            
    return source_to_target_map
            
def nodes_without_match(
    G_source,
    G_target,
    verbose = False):
    
    curr_map = nxu.node_map(
        G_source,
        G_target)
    
    no_map = np.array([k for k,v in curr_map.items() if v is None])
    
    if verbose:
        print(f"Nodes with no match = {len(no_map)}")
        
    return no_map

def skeleton_from_node(
    G,
    n):
    """
    Ex: nxu.skeleton_from_node(G_obj,n = "L0_16")
    """
    skeleton_points = G.nodes[n]["skeleton_data"]
    sk_idx = np.vstack([
        np.arange(0,len(skeleton_points)-1),
        np.arange(1,len(skeleton_points))]).T
    return skeleton_points[sk_idx]
    
import numpy_utils as nu
def skeleton_coordinates_from_G(
    G,
    include_upstream_node = False,
    verbose = False,
    return_node_names = False,
    ):
    """
    Purpose: To output the skeleton points of a graph
    and the node names those skeleton points came from
    
    Ex: 
    sk_pts,sk_names = nxu.skeleton_coordinates_from_G(
        G_ax,
    return_node_names=True,
    verbose = True)
    """

    total_skeleton_points = []
    total_node_names = []

    for n in nxu.limb_branch_nodes(G):
        skeleton_points = G.nodes[n]["skeleton_data"]
        if not include_upstream_node:
            upstream_coordinate = G.nodes[n]["endpoint_upstream"].reshape(-1,3)
            skeleton_points = nu.setdiff2d(skeleton_points,upstream_coordinate)

        node_names = np.array([n]*len(skeleton_points))

        if verbose:
            print(f"Node {n} has {len(skeleton_points)} skeleton points")

        total_skeleton_points.append(skeleton_points)
        total_node_names.append(node_names)

    if len(total_skeleton_points) > 0:
        total_skeleton_points = np.vstack(total_skeleton_points)
        total_node_names = np.concatenate(total_node_names)
    else:
        total_skeleton_points = np.array([])
        total_node_names = np.array([])

    if verbose:
        print(f"Total skeleton points = {len(total_skeleton_points)}")

    if return_node_names:
        return total_skeleton_points,total_node_names
    else:
        return total_skeleton_points
    
    
from pykdtree.kdtree import KDTree
import pandas as pd

def split_location_node_map_df(
    G,
    split_locations,
    G_target = None,
    nodelist = None,
    distance_max = 5000,
    error_if_no_one_match = True,
    error_on_non_unique_node_names = True,
    verbose = False,
    verbose_loop = False,
    ):
    """
    Purpose: Want to map split locaitons to the nodes that
    were cut off due to proofreading with what rule was used to cut off

    Things can leverage:
        The limb name
        Nodes that shouldn't be mapped

    Pseudocode: 
    Iterate through all filters and limbs to work
    with the split locations and possible nodes to match with



    """

    split_info = []

    if G_target is not None and nodelist is None:
        nodelist = nxu.nodes_without_match(
            G_source = G,
            G_target=G_target,
            verbose = verbose )

    if nodelist is not None: 
        G_search = G.subgraph(nodelist).copy()
    else:
        G_search = G

    split_counter = 0
    for filter_name,filter_splits in split_locations.items():
        if type(filter_splits) != dict:
            filter_splits = {"L-1":filter_splits}

        for limb_name,split_loc in filter_splits.items():
            G_limb = nxu.limb_graph(G_search,limb_name)

#             sk_pts,sk_names = nxu.skeleton_coordinates_from_G(
#                 G_limb,
#                 return_node_names=True,
#                 verbose = False)

            if len(G_limb) == 0:
                raise Exception("Empty graph")

            for i,s in enumerate(split_loc):
                
                local_split = dict(limb_split_idx = i,
                                   limb_name = limb_name,
                                   coord = s,
                                  filter_name = filter_name)
                
                """
                4/13: Want to determine the node that has the lowest average distance
                
                Pseudocode:
                1) Iterate through all of the nodes
                2) Build a KDTree on their skeleton points
                3) Find the average distance
                4) Add to name and distance to list to find min later
                
                """
                sk_names = []
                dist_from_split = []
                for n in nxu.limb_branch_nodes(G_limb):
                    '''    
                    dist_from_split = np.linalg.norm(sk_pts - s,axis = 1)
                    '''
                    n_kd = KDTree(G_limb.nodes[n]["skeleton_data"])
                    dist,_ = n_kd.query(s.reshape(-1,3))
                    sk_names.append(n)
                    dist_from_split.append(np.mean(dist))
                    
                sk_names = np.array(sk_names)
                dist_from_split = np.array(dist_from_split)

                min_dist = np.min(dist_from_split)
                min_arg = np.where(dist_from_split == min_dist)[0]
                min_nodes = sk_names[min_arg]

                if len(min_nodes) != 1 and error_if_no_one_match:
                    raise Exception(f"Split {i} ({s}) had min_dist ({min_dist}) had more than one match ({min_nodes})")

                if error_if_no_one_match:
                    min_nodes = min_nodes[0]

                if min_dist > distance_max:
                    raise Exception(f"Split {i} ({s}) had min_dist ({min_dist}) that was greater than distance_max ({distance_max})")

                local_split["node"] = min_nodes
                local_split["min_distance"] =  min_dist

                split_info.append(local_split)
                split_counter += 1

                if verbose_loop:
                    print(f"local_split = {local_split}")

    split_df = pd.DataFrame.from_records(split_info)
    
    
    if error_on_non_unique_node_names and len(split_df) > 0:
        if len(split_df["node"].unique()) != len(split_df):
            raise Exception("Not unique nodes")
    
    if verbose:
        print(f"Total split mappings: {len(split_df)}")
    
    return split_df

def nodes_with_auto_proof_filter(
    G,
    return_filter_names=False,
    verbose = False,
    ):
    return xu.nodes_with_non_none_attributes(
        G,
        attribute_name=nxu.auto_proof_filter_name_default,
        return_attribute_value=return_filter_names,
        verbose=verbose)

def set_auto_proof_filter_attribute(
    G,
    split_df = None,
    split_locations = None,
    inplace = True,
    filter_attribute_name = None,
    default_value = None,
    verbose = False,
    error_on_non_unique_node_names=True,
    filter_axon_on_dendrite_splits = True,
    **kwargs
    ):
    """
    Purpose: To take the split locations df and
    to label the nodes with the right error
    """
    
    if split_df is None:
        split_df = nxu.split_location_node_map_df(
            G = G,
            split_locations=split_locations,
            verbose = verbose,
            error_on_non_unique_node_names=error_on_non_unique_node_names,
            **kwargs
            )

    if filter_attribute_name is None:
        filter_attribute_name = nxu.auto_proof_filter_name_default

    if not inplace:
        G = copy.deepcopy(G)

    xu.set_node_attribute(
        G,
        attribute_name = filter_attribute_name,
        attribute_value = default_value
    )

    if len(split_df) > 0:
        filter_names = split_df["filter_name"]
        nodes_to_label = split_df["node"]

        for n,f in zip(nodes_to_label,filter_names):
            if error_on_non_unique_node_names:
                G.nodes[n][filter_attribute_name] = f
            else:
                if G.nodes[n][filter_attribute_name] is None:
                    G.nodes[n][filter_attribute_name] = []
                G.nodes[n][filter_attribute_name].append(f)

    if verbose:
        print(f"nodes with filter attributes")
        nxu.nodes_with_auto_proof_filter(
            G,
            verbose = True)
        
    if filter_axon_on_dendrite_splits:
        if verbose:
            print(f"filtering the axon on dendrite splits")
        G = nxu.filter_axon_on_dendrite_splits_to_most_upstream(G)

    return G

def filter_axon_on_dendrite_splits_to_most_upstream(G):
    return G

def segment_id_from_G(G,return_split_index = True):
    seg_id = G.graph["segment_id"]
    if return_split_index:
         return seg_id,G.graph["split_index"]
    else:
        return seg_id


# ---------- For filtering parts of the graph ---------
def soma_connected_nodes(
    G,
    ):
    return list(G[nxu.soma_node_name_global].keys())


import numpy as np
def most_upstream_nodes(
    G,
    nodes,
    verbose = False,
    return_downstream_count = True,
    ):
    """
    Purpose: To get a count of the number of
    downstream nodes
    
    Ex: 
    nxu.most_upstream_nodes(
        G,
        nodes = ["L1_10","L1_2","L1_8","L0_20"],
        verbose = True
    )
    """
    nodes = np.array(nodes)

    down_count = np.array([xu.n_all_downstream_nodes(G,k) for k in nodes])
    down_count_idx_sorted = np.flip(np.argsort(down_count))
    nodes_sorted = nodes[down_count_idx_sorted]
    down_count_sorted = down_count[down_count_idx_sorted]

    if verbose:
        print(f"nodes_sorted = {nodes_sorted}")
        print(f"down_count_sorted = {down_count_sorted}")
    
    if return_downstream_count:
        return nodes_sorted,down_count_sorted
    else:
        return nodes_sorted
    
def limb_graphs_from_soma_connected_nodes(
    G,
    verbose = False,
    plot = False,
    ):
    """
    Purpose: To get all of the limb graphs
    defined by those boardering the soma

    Pseucode: 
    1) Get all of the soma connected nodes
    2) For each soma connected node: get the connected subgraph
    
    Ex: 
    limb_graphs = nxu.limb_graphs_from_soma_connected_nodes(
    G,
    verbose = True,
    plot = False,
    )
    """
    soma_conn_nodes = nxu.soma_connected_nodes(G)
    if verbose:
        print(f"soma_conn_nodes = {soma_conn_nodes}")

    total_limb_graphs = []
    for node in soma_conn_nodes:
        if verbose:
            print(f"-- Working on soma node {node}")
        limb_graph = nxu.limb_graph(G,most_upstream_node=node)

        if plot:
            most_up_node = xu.most_upstream_node(limb_graph)
            print(f"Plotting limb with most upstream node {most_up_node}")
            nxu.draw_tree(limb_graph)
            print(f"\n\n")

        total_limb_graphs.append(limb_graph)

    return total_limb_graphs


def distance_from_node_to_soma(
    G,
    node,
    include_self_distance = False,
    distance_attribute = "skeletal_length",
    destination_node = "S0",
    verbose = False,
    return_path= False,
    ):
    """
    Purpose: To find the distance of a node from soma
    (both with inclusion and without of its own skeletal length)

    Pseudocode: 
    1) 
    """


    total_lengths = []
    node_path = []
    
    if verbose:
        print(f"Finding distance from {node} to {destination_node}")

    if include_self_distance:
        total_lengths.append(G.nodes[node][distance_attribute])
        node_path.append(node)

    upstream_node = xu.upstream_node(G,node)

    while upstream_node != destination_node:
        if verbose:
            print(f"Working on upstream node {upstream_node}")

        node_path.append(upstream_node)
        total_lengths.append(G.nodes[upstream_node][distance_attribute])
        upstream_node = xu.upstream_node(G,upstream_node)

    path_length = np.sum(total_lengths)
    if verbose:
        print(f"path_length = {path_length}")
        print(f"Path to soma: {node_path}")
        print(f"Path {distance_attribute}: {total_lengths}")

    if return_path:
        return path_length,node_path
    else:
        return path_length
    
def distance_upstream_from_soma(
    G,
    node,
    verbose = False,
    **kwargs
    ):
    """
    nxu.distance_upstream_from_soma(
    G,
    "L0_19",
    verbose = True,   
    )
    """
    
    return distance_from_node_to_soma(
    G,
    node,
    include_self_distance = False,
    verbose = verbose,
    return_path= False,
    **kwargs
    )

def distance_downstream_from_soma(
    G,
    node,
    verbose = False,
    **kwargs
    ):
    """
    Ex: 
    nxu.distance_downstream_from_soma(
        G,
        "L0_19",
        verbose = True,   
    )
    """
    
    return distance_from_node_to_soma(
    G,
    node,
    include_self_distance = True,
    verbose = verbose,
    return_path= False,
    **kwargs
    )
    
import pandas as pd
def distance_from_soma_df(
    G,
    nodes = None,
    distance_type = "upstream",):
    """
    Purpose: Find all the soma distances of 
    all the nodes
    """
    if nodes is None:
        nodes = nxu.limb_branch_nodes(G)

    dist_dict = [{"node":n,"soma_distance":getattr(nxu,f"distance_{distance_type}_from_soma")(G,n)}
                for n in nodes]

    dist_dict = pd.DataFrame.from_records(dist_dict)

    return dist_dict

def distance_upstream_from_soma_df(
    G,
    nodes = None,):
    
    return nxu.distance_from_soma_df(
    G,
    nodes = nodes,
    distance_type = "upstream",)

def distance_downstream_from_soma_df(
    G,
    nodes = None,):
    
    return nxu.distance_from_soma_df(
    G,
    nodes = nodes,
    distance_type = "downstream",)

def nodes_distance_query_from_soma(
    G,
    distance_threshold,
    within_distance = True,
    distance_type = "upstream",
    nodes=None,
    return_subgraph=False,
    return_soma_with_sugraph = True,
    verbose = False,
    maintain_skeleton_connectivity = True,
    **kwargs
    ):
    """
    Purpose: Find all the nodes within
    a certain distance or farther away than
    a certain distance from soma
    """
    dist_df = getattr(nxu,f"distance_{distance_type}_from_soma_df")(G,nodes=nodes)
    
    
    
    if within_distance:
        query = f"soma_distance <= {distance_threshold}"
        query_descriptor = "closer"
    else:
        query = f"soma_distance > {distance_threshold}"
        query_descriptor = "farther"

    filt_df = dist_df.query(query)
    
    filt_nodes = filt_df["node"].to_list()
    if verbose:
        print(f"Nodes {distance_type} dist {query_descriptor} than {distance_threshold} from soma ({len(filt_nodes)}):\n{np.array(filt_nodes)} ")
        
        
    if return_subgraph:
        
        if return_soma_with_sugraph:
            filt_nodes += [nxu.soma_node_name_global]
            
        # Old way:
        #G_sub = G.subgraph(filt_nodes).copy()
        
        node_to_delete = np.setdiff1d(list(G.nodes()),filt_nodes)
        if len(node_to_delete) > 0:
            G_sub = nxu.remove_node(
            G,
            node_to_delete,
            verbose= False,
            inplace = False,
            maintain_skeleton_connectivity=maintain_skeleton_connectivity,
            **kwargs)
        else:
            G_sub = G
       
        
        return G_sub
    else:
        return filt_nodes

def nodes_within_distance_from_soma(
    G,
    distance_threshold,
    distance_type = "upstream",
    nodes=None,
    return_subgraph=False,
    verbose = False,
    ):
    
    """
    Ex: 
    distance_threshold = 1000
    curr_nodes = nxu.nodes_within_distance_from_soma(G,distance_threshold = distance_threshold,verbose = True)
    """
    
    return_G =  nxu.nodes_distance_query_from_soma(
    G,
    distance_threshold,
    within_distance = True,
    distance_type = distance_type,
    nodes=nodes,
    return_subgraph=return_subgraph,
    verbose = verbose,
    )
    
    
    return return_G
    
def nodes_farther_than_distance_from_soma(
    G,
    distance_threshold,
    distance_type = "upstream",
    nodes=None,
    return_subgraph=False,
    verbose = False,
    ):
    
    return nxu.nodes_distance_query_from_soma(
    G,
    distance_threshold,
    within_distance = False,
    distance_type = distance_type,
    nodes=nodes,
    return_subgraph=return_subgraph,
    verbose = verbose,
    )

def nodes_within_distance_upstream_from_soma(
    G,
    distance_threshold,
    nodes=None,
    return_subgraph=False,
    verbose = False,
    ):
    
    return nxu.nodes_within_distance_from_soma(
    G,
    distance_threshold,
    distance_type = "upstream",
    nodes=nodes,
    return_subgraph=return_subgraph,
    verbose = verbose,
    )


import numpy as np
def soma_filter_by_complete_graph(
    G,
    inplace = False,
    verbose = False,
    connect_previous_touching_soma_nodes = False,
    plot = False):
    """
    Problem: Want to resolve the soma so it does
    not affect...

    Solution 1: Delete the soma network and then 
    connect all the 

    Pseudocode:
    1) Connect all the soma connecting nodes
    2) Get the subgraph of only the limb branches
    
    Ex: 
    G_no_soma = nxu.soma_filter_by_complete_graph(G_filt,plot=True)
    nxu.draw_tree(G_no_soma)
    """
    if not inplace:
        G = copy.deepcopy(G)
    
    if connect_previous_touching_soma_nodes:
        soma_conn_nodes = nxu.soma_connected_nodes(G)
        if verbose:
            print(f"soma_conn_nodes ({len(soma_conn_nodes)})= {soma_conn_nodes}")

        all_conn = np.concatenate([[(k,v) for v in soma_conn_nodes if v != k] for k in soma_conn_nodes])

        if verbose:
            print(f"New connections ({len(all_conn)}) = {all_conn}")

        G.add_edges_from(all_conn)
        
    return_G = nxu.limb_branch_subgraph(G)
    
    if plot:
        nx.draw(nx.Graph(return_G),with_labels = True)
        
    return return_G


def nodes_with_auto_proof_filter_type(
    G,
    filter_type,
    filter_feature = "auto_proof_filter",
    verbose = False,
    ):
    """
    Purpose: To get all nodes with a certain filter marking 
    
    Ex: 
    nxu.nodes_with_auto_proof_filter_type(G,filter_type="axon_on_dendrite",verbose = True)
    """
    filter_df = xu.node_df(G).query(f"{filter_feature} == {filter_feature}")
    nodes_with_filter_type = filter_df[filter_df[filter_feature].str.contains(filter_type)][xu.upstream_name].to_numpy()
    if verbose:
        print(f"nodes with {filter_type} = {nodes_with_filter_type}")
        
    return nodes_with_filter_type

def clear_nodes_auto_proof_filter_feature(
    G,
    nodes,
    filter_feature = "auto_proof_filter",
    verbose=False,
    ):
    for n in nodes:
        curr_filter_name = G.nodes[n][filter_feature]
        if verbose:
            print(f"Clearning {curr_filter_name} for node {n} ")
        G.nodes[n][filter_feature] = None
    return G

def filter_axon_on_dendrite_splits_to_most_upstream(
    G,
    verbose=False,
    inplace = False,
    filter_out_only_if_parent_in_split = False,
    **kwargs):
    """
    Purpose: To reduce the axon on dendrite merges to only those
    that are the most upstream of the group
    
    Pseudocode:
    1) Get all of the nodes that are in axon on dendrite mergers
    2) For each node:
    a. Get either just the parent or all of the upstream node
    b. add node to list ot be cleared if has upstream axon on dendrite
    
    Ex: 
    G_filt = filter_axon_on_dendrite_splits_to_most_upstream(G,verbose = True)
    nxu.nodes_with_auto_proof_filter(G_filt,"axon_on_dendrite")
    
    """
    axon_on_dendrite_nodes = nxu.nodes_with_auto_proof_filter_type(G,"axon_on_dendrite")
    
    if not inplace:
        G = copy.deepcopy(G)
    
    if verbose:
        print(f"axon_on_dendrite_nodes ({len(axon_on_dendrite_nodes)}) = {axon_on_dendrite_nodes}")
        
    nodes_to_clear = []
    for n in axon_on_dendrite_nodes:
        if filter_out_only_if_parent_in_split:
            up_nodes = [xu.upstream_node(G,n)]
        else:
            up_nodes = xu.all_upstream_nodes(G,n)
            
        ax_on_dendr_upstream = np.intersect1d(axon_on_dendrite_nodes,up_nodes)
        if len(ax_on_dendr_upstream) > 0:
            if verbose:
                print(f"Removing node {n} because had upstream axon on dendrite: {ax_on_dendr_upstream}")
            nodes_to_clear.append(n)
            
    G = nxu.clear_nodes_auto_proof_filter_feature(G,nodes_to_clear,verbose = verbose)
        
    
    if verbose:
        axon_on_dendrite_nodes = nxu.nodes_with_auto_proof_filter_type(G,"axon_on_dendrite")
        print(f"axon_on_dendrite_nodes AFTER FILTERING ({len(axon_on_dendrite_nodes)}) = {axon_on_dendrite_nodes}")
        
    return G
    
    
import numpy as np
def fix_flipped_skeletons(
    G,
    verbose = False,):
    """
    Purpose: To fix the skeleton data
    in Graph objects if they are
    not aligned correctly
    
    Ex: 
    segment_id,split_index = 864691135162621741,0

    G_obj = hdju.graph_obj_from_auto_proof_stage(
        segment_id=segment_id,
        split_index=split_index,
    )
    
    G_obj = fix_flipped_skeletons(G_obj,verbose = True)
    """
    
    node_flipped = []
    for n in G.nodes():
        if "s" in n.lower():
            continue
        curr_dict = G.nodes[n]
        endpt_up = curr_dict["endpoint_upstream"]
        sk_up = curr_dict["skeleton_data"][0]
        if not np.array_equal(endpt_up,sk_up):
            G.nodes[n]["skeleton_data"] = np.flip(curr_dict["skeleton_data"],axis=0)
            node_flipped.append(n)
            
    if verbose:
        print(f"Nodes with skeleton flipped: {node_flipped}")
        
    return G

def soma_only_graph(
    G,
    soma_node_name=soma_node_name_global):
    """
    Purpose: To check if only a soma node is in the nodes
    """
    node_names = [n for n in G.nodes()]
    
    if len(node_names) == 1:
        if node_names[0] == soma_node_name:
            return True
    else:
        return False
    
# --------- For outputing graph attributes --------
def compartment_from_node(G,n=None):
    """
    Purpose: To get the compartment from
    a node in a graph
    """
    if "graph" in str(type(G)).lower():
        G = G.nodes[n]
        
    return G["compartment"].replace("_","")

def width_from_node(
    G,
    n,
    verbose = False):
    """
    Purpose: To get the width of a certain node
    
    Ex: 
    import neuron_nx_utils as nxu
    nxu.width_from_node(
     G = G_obj,
     n = "L0_5",
        verbose = True

    )
    """
    curr_key = G.nodes[n]
    width = None
    if len(curr_key["width_new"]) > 0:
        width = curr_key["width_new"].get("no_spine_median_mesh_center",None)
    if width is None:
        if verbose:
            print(f"Getting width from default width (not new width)")
        width = key.get("width",None)

    return width


def attribute_graph_from_graph_obj(
    G,
    attribute,
    ids_name = None,
    verbose = False,
    return_attribute_nodes = False,
    return_attribute_nodes_to_branch_dict = True,
    return_upstream_dist = False,
    ):
    """
    Purpose: To convert a graph object into 
    a graph where an attribute of the graph 
    object are the nodes along with the branch points.
    The edges between nodes will be the upstream distance.

    Application: Will help find the distances between the attributes
    """

    attribute_name = f"{attribute}_data"

    if ids_name is None:
        if attribute == "synapse":
            ids_name = f"syn_id"
        else:
            ids_name = f"{attribute}_id"

    graph_edges = []
    graph_weights = []
    ids_counter = 0
    attribute_ids = []
    attribute_ids_branch_dict = dict()
    upstream_dist = {v:None for v in G.nodes()}
    nodes = list(G.nodes())

    for n in nodes:
        key = G.nodes[n]
        curr_dict = dict()

        upstream_node = xu.upstream_node(G,n)
        if upstream_node is None:
            upstream_node = f"L-1"

        comp = compartment_from_node(key)
        attribute_data = key[attribute_name]

        #if len(attribute_data) > 0:
        upstream_data = np.array([k["upstream_dist"] for k in attribute_data])

        # generates the ids for the current branch
        ids_data = np.array([k[ids_name] for k in attribute_data])
        if None in ids_data:
            ids_data = np.arange(ids_counter,ids_counter + len(attribute_data)).astype('int')
            ids_counter += len(attribute_data)
            
        
        attribute_ids.append(ids_data)
        attribute_ids_branch_dict.update({str(k):n for k in ids_data})
        upstream_dist.update({str(k):upstream_data[j] for j,k in enumerate(ids_data)})

        order_idx = np.argsort(upstream_data)
        upstream_data_sort = upstream_data[order_idx]
        ids_data_sort = ids_data[order_idx]

        upstream_data_sort = np.hstack([[0],upstream_data_sort,[key["skeletal_length"]]])
        ids_data_sort = np.hstack([[upstream_node],ids_data_sort,[n]])

        edge_weights = upstream_data_sort[1:] - upstream_data_sort[:-1]
        edge_weights[edge_weights <= 0] = 0
        edges = np.vstack([ids_data_sort[:-1],ids_data_sort[1:]]).T
        #weighted_edges = np.vstack([edges,edge_weights]).T

        if verbose:
            print(f" --> node {n}: # of edges = {len(edges)}")
            #print(f"weighted_edges = {weighted_edges}")

        graph_edges.append(edges)
        graph_weights.append(edge_weights)


    output_G = nx.Graph()

    if len(graph_edges) > 0:
        graph_edges = np.vstack(graph_edges)
        graph_weights = np.hstack(graph_weights)
        output_G = xu.edges_and_weights_to_graph(graph_edges,
                              weights_list=graph_weights)
        
    if len(attribute_ids) > 0:
        attribute_ids = np.hstack(attribute_ids)
        
    attribute_ids = np.array(attribute_ids).astype('int').astype("str")
        
    if verbose:
        print(f"Total Graph stats:")
        xu.print_node_edges_counts(output_G)

    if ((not return_attribute_nodes) and (not return_attribute_nodes_to_branch_dict)
        and (not return_upstream_dist)):
        return output_G
    return_value =[output_G]
    
    if return_attribute_nodes:
        return_value.append(attribute_ids)
    if return_attribute_nodes_to_branch_dict:
        return_value.append(attribute_ids_branch_dict)
        
    if return_upstream_dist:
        return_value.append(upstream_dist)
    
    return return_value
    
    
import matplotlib.pyplot as plt
def plot_inter_attribute_intervals(
    inter_attribute_info,
    title = None,
    bins=50,
    attribute = "attribute",
    verbose = False):
    """
    Purpose: To plot the histograms of the
    closest attribute arrays
    """
    if type(inter_attribute_info) != dict:
        inter_attribute_info = {1:inter_attribute_info}
        
    if type(inter_attribute_info[1]) == dict:
        if verbose:
            print(f"Combining the branch specific dicts")
        inter_attribute_info = {k:np.concatenate(list(v.values())) for k,v in inter_attribute_info.items()}
        
    fig,ax = plt.subplots(1,1,)
    for i in range(len(inter_attribute_info)):
        ax.hist(
            np.array(inter_attribute_info[i+1])/1000,
            label=f"{i+1} Hop",
            alpha = 0.5,
            bins = bins)

        if title is None:
            title = f"Closest {attribute}"
        ax.set_title(title)
        ax.set_xlabel(f"Distance (um)")
        ax.set_ylabel(f"Count")
        ax.legend()
        
    return ax

import networkx_utils as xu
import time
from tqdm_utils import tqdm
def inter_attribute_intervals_from_G(
    G,
    attribute,
    n_closest_neighbors = 1,
    default_value = -1,
    debug_time = False,
    plot = False,
    verbose = False,
    separate_branches = True,
    return_upstream_dist = True,
    ):
    """
    Purpose: To find the k inter-attribute
    distances for the attributes on the graph

    1) Turn the graph into an attribute graph
    2) For every attribute id:
        Remove the attribute id from the total list
        a. For 1 to k (the number of attributes away):
           if list empty then add distance of -1 and continue
           Calculate the closest neighbors from remaining nodes and get distance

           Save the distance
           Remove that closest neighbor from the list

    3) Return the lists of closest distances
    """

    (G_disc,
     att_nodes,
     att_to_branch_dict,
     att_to_upstream_dist)= nxu.attribute_graph_from_graph_obj(
        G,
        attribute = attribute,
        verbose = verbose,
        return_attribute_nodes = True,
        return_attribute_nodes_to_branch_dict = True,
        return_upstream_dist = True,
    )

    upstream_dists = {v:[] for v in G.nodes()}
    if not separate_branches:
        closest_neighbors = {k:[] for k in range(1,n_closest_neighbors + 1)}
        
    else:
        closest_neighbors = {k:{v:[] for v in G.nodes()} for k in range(1,n_closest_neighbors + 1)}
        

    for n in tqdm(att_nodes):
        local_nodes = att_nodes.copy()
        local_nodes = local_nodes[local_nodes != n]
        previous_closest = None
        
        upstream_dists[att_to_branch_dict[n]].append(att_to_upstream_dist[n])
        
        for i in closest_neighbors:
            if previous_closest is not None:
                local_nodes = local_nodes[local_nodes != previous_closest]
            if len(local_nodes) == 0:
                default_value
                short_path_dist = default_value
                n2 = None

            else:
                if debug_time:
                    st = time.time()
                short_path_dist,n1,n2 = xu.shortest_path_between_two_sets_of_nodes(
                    G_disc,
                    node_list_1 = [n],
                    node_list_2 = local_nodes,
                    return_node_pairs=True,
                    return_path_distance=True
                )

                if debug_time:
                    print(f"Time for calculating shortest path = {time.time() - st}")
                    st = time.time()


            if verbose:
                print(f"{i}th Closest neighbor of {n} was {n2} with path distance {short_path_dist}")
            previous_closest = n2
            
            if not separate_branches:
                closest_neighbors[i].append(short_path_dist)
            else:
                closest_neighbors[i][att_to_branch_dict[n]].append(short_path_dist)

    if not separate_branches:
        closest_neighbors = {k:np.array(v) for k,v in closest_neighbors.items()}
    else:
        for k in closest_neighbors:
            for v in closest_neighbors[k]:
                closest_neighbors[k][v] = np.array(closest_neighbors[k][v])
        

    if plot: 
        nxu.plot_inter_attribute_intervals(
            closest_neighbors,
            attribute=attribute)
        plt.show()
        
    if return_upstream_dist:
        return closest_neighbors,upstream_dists
    else:
        return closest_neighbors


def n_data_attribues(G,attribute,n=None):
    """
    Purpose: To get the number of data attributes
    belonging to a branch
    """
    if n is None:
        n = [k for k in list(G.nodes()) if "S" not in k]
    else:
        n = [n]
    return np.sum([len(G.nodes[n1][f"{attribute}_data"]) for n1 in n])



import networkx_utils as xu
import networkx as nx
import numpy_utils as nu

def inter_attribute_intervals_dict_from_neuron_G(
    G,
    attribute,
    dendrite_only = True,
    remove_starter_branches = True,
    #return_attribute_density = True,
    #return_compartment = True
    verbose = True,
    n_closest_neighbors = 3,
    branch_features_to_add = (
        "mesh_volume",
        'n_synapses_head',
        'n_synapses_neck',
        'n_synapses_no_head',
        'n_synapses_post',
        'n_synapses_pre',
        'n_synapses_shaft',
        'n_synapses_spine',
        "total_spine_volume",
        "soma_distance_euclidean",
        'parent_skeletal_angle',
        'siblings_skeletal_angle_max',
        'width_upstream',
         'width_downstream',
        ),
    ):
    """
    Purpose: Want to build an inter attribute 
    distance for all branches in a neuron
    """
    if type(attribute) == tuple:
        attribute = list(attribute)
    attribute= nu.convert_to_array_like(attribute)
    data_dicts = []


    if not nxu.soma_only_graph(G):
        if dendrite_only:
            G = nxu.dendrite_subgraph(G)

        if remove_starter_branches:
            G = nxu.remove_small_starter_branches(
                    G,
                    verbose = False,
                    maintain_skeleton_connectivity = True)


        limb_graphs,limb_idxs = nxu.all_limb_graphs(G,return_idxs=True)

        if verbose:
            print(f"# of limb graphs: {len(limb_graphs)}")

        for idx,G_limb in zip(limb_idxs,limb_graphs):

            if verbose:
                print(f"\n--------Working on Limb {idx}---------")

            limb_dicts = dict()
            for n in G_limb.nodes():
                curr_dt = dict(branch=n,
                                   compartment=nxu.compartment_from_node(G_limb,n),
                                     skeletal_distance_to_soma = nxu.distance_upstream_from_soma(G,node=n),
                                   skeletal_length = G_limb.nodes[n]["skeletal_length"],
                                    width = nxu.width_from_node(G_limb,n))
                if branch_features_to_add is not None:
                    for k in branch_features_to_add:
                        curr_dt[k] = G_limb.nodes[n][k]
                        
                limb_dicts[n] = curr_dt
                              
            
            
            for att in attribute:
                for n in limb_dicts:
                    limb_dicts[n][f"n_{att}"] = nxu.n_data_attribues(G_limb,attribute=att,n=n)
                (inter_dict_att,
                 upstream_dist)= nxu.inter_attribute_intervals_from_G(
                    G_limb,
                    attribute=att,
                    n_closest_neighbors=n_closest_neighbors,
                    separate_branches = True,
                    plot = False,
                     return_upstream_dist = True
                )

                # Add the attribute data to the dicts
                for i,data in inter_dict_att.items():
                    for node_name,node_array in data.items():
                        limb_dicts[node_name][f"{att}_intervals_{i}"] = node_array
                        
                for n in limb_dicts.keys():
                    limb_dicts[n][f"{att}_upstream_dist"] = np.array(upstream_dist[n])

            data_dicts += list(limb_dicts.values())


    return data_dicts

    

import neuron_nx_utils as nxu