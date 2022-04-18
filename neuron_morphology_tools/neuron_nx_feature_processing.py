import neuron_nx_utils as nxu
import numpy_utils as nu
import networkx as nx

import networkx_utils as xu
import numpy as np
import copy

def features_list(
    G,
    limb_branch_features = True,
    features_to_ignore = ("u",),
    verbose = False):
    """
    Purpose: Find all of the current features
    
    Ex: 
    import neuron_nx_feature_processing as nxf
    nxf.features_list(G)
    """
    if limb_branch_features:
        G = nxu.limb_branch_subgraph(G)
    node_df = xu.node_df(G)
    current_features = list(node_df.columns)
    if features_to_ignore is not None:
        current_features = list(np.setdiff1d(current_features,features_to_ignore))

    if verbose:
        print(f"current_features = {current_features}")
        
    return current_features

def add_node_feature(
    G,
    feature_func,
    feature_name = None,
    nodes=None,
    inplace = True,
    verbose=False,
    default_value = None,
    feature_value_dict = None,
    verbose_loop = False,
    #skip_if_exists = True,
    ):
    
    """
    Puprose: Will apply an feature_func
    to a node based on the current node values
    """
    
    if nodes is None:
        nodes = nxu.limb_branch_nodes(G)
        
    if not inplace:
        G = copy.deepcopy(G)
    
    feature_func = nu.convert_to_array_like(feature_func)
    
#     if verbose:
#         print(f"feature functions")
#         print(f"{[k.__name__ for k in feature_func]}")
        
    for att_func in feature_func:
        if "str" in str(type(att_func)):
            att_func = getattr(nxf,att_func)
            
        if feature_name is None:
            curr_name = str(att_func.__name__)
        else:
            curr_name = feature_name

        if verbose:
            print(f"\n\n---Setting {curr_name}, att_func ={att_func}")

        for n in nodes:
            if feature_value_dict is not None:
                attr_value = feature_value_dict.get(n,curr_name,default_value)
            else:
                try:
                    attr_value = att_func(G.nodes[n])
                except:
                    if default_value is not None:
                        attr_value = default_value
                    else:
                        raise Exception("")
            if verbose_loop:
                print(f"For node {n}, {curr_name} = {attr_value}")
            G.nodes[n][curr_name] = attr_value
        
    return G


def skeleton_vector_upstream_x(node_dict):
    return node_dict[f"skeleton_vector_upstream"][0]
def skeleton_vector_upstream_y(node_dict):
    return node_dict[f"skeleton_vector_upstream"][1]
def skeleton_vector_upstream_z(node_dict):
    return node_dict[f"skeleton_vector_upstream"][2]
def skeleton_vector_downstream_x(node_dict):
    return node_dict[f"skeleton_vector_downstream"][0]
def skeleton_vector_downstream_y(node_dict):
    return node_dict[f"skeleton_vector_downstream"][1]
def skeleton_vector_downstream_z(node_dict):
    return node_dict[f"skeleton_vector_downstream"][2]

def skeleton_vector_downstream_theta(node_dict):
    return nu.polar_3D_from_cartesian(*node_dict[f"skeleton_vector_downstream"])[1]
def skeleton_vector_downstream_phi(node_dict):
    return nu.polar_3D_from_cartesian(*node_dict[f"skeleton_vector_downstream"])[2]

def skeleton_vector_upstream_theta(node_dict):
    return nu.polar_3D_from_cartesian(*node_dict[f"skeleton_vector_upstream"])[1]
def skeleton_vector_upstream_phi(node_dict):
    return nu.polar_3D_from_cartesian(*node_dict[f"skeleton_vector_upstream"])[2]

def width_no_spine(node_dict):
    return node_dict["width_new"]["no_spine_median_mesh_center"]

def axon_label(node_dict):
    return int(node_dict["axon_compartment"] == "axon")

def add_skeleton_vector_features(
    G,
    use_polar_coords = True,
    verbose = False,
    **kwargs
    ):
    
    

        
    if not use_polar_coords:
        attr_funcs = [skeleton_vector_upstream_x,
        skeleton_vector_upstream_y,
        skeleton_vector_upstream_z,
        skeleton_vector_downstream_x,
        skeleton_vector_downstream_y,
        skeleton_vector_downstream_z,]
    else:
        if verbose:
            print(f"Using polar coordinates")
        attr_funcs = [skeleton_vector_downstream_theta,
        skeleton_vector_downstream_phi,
        skeleton_vector_upstream_theta,
        skeleton_vector_upstream_phi,]
    G_new_feats = nxf.add_node_feature(
        G,
        feature_func=attr_funcs,
        **kwargs
        )
    
    return G_new_feats

def add_any_missing_node_features(
    G,
    features,
    verbose = False,
    inplace = False,
    ):
    """
    Purpose:
    1) Check that all the features are requested
    2) Generate the features that are not
    """

    curr_features = nxf.features_list(G)
    features_not_computed = np.setdiff1d(features,curr_features)

    if verbose:
        print(f"features_not_computed = {features_not_computed}")

    G_new = nxf.add_node_feature(
        G,
        feature_func=features_not_computed,
        verbose = verbose,
        inplace = inplace
    )

    return G_new


def filter_G_features(
    G,
    features,
    inplace = False,
    verbose = False,
    ):
    """
    Purpose: To reduce a networkx graph to a 
    certain number of features (and all other superflous features are deleted)
    
    Ex: 
    axon_vs_dendrite_features = [
    "mesh_volume",
    "n_spines",
    "total_spine_volume",
    "n_synapses_post",
    "n_synapses_pre",
    #"n_synapse_head",
    #"parent_skeletal_angle",
    "skeletal_length",
    "skeleton_vector_upstream_theta",
    "skeleton_vector_upstream_phi",
    "skeleton_vector_downstream_theta",
    "skeleton_vector_downstream_phi",
    "width_upstream",
    "width_no_spine",
    "width_downstream",
    ]
    
    G_feat_filt = nxf.filter_G_features(
        G,
        features=axon_vs_dendrite_features,
        inplace = False,
        verbose = True,
    )

    xu.node_df(G_feat_filt)
    """


    if not inplace:
        G = copy.deepcopy(G)


    G_ret = nxf.add_any_missing_node_features(
        G,
        features = features,
        verbose = verbose
    )

    if verbose:
        print(f"Number of features after adding missing ones = {len(xu.node_df(G_ret).columns)}")

    G_ret = xu.delete_node_attributes(G_ret,attributes_not_to_delete=features)

    if verbose:
        print(f"Number of features after adding missing ones = {len(xu.node_df(G_ret).columns)}")
    
    return G_ret





import neuron_nx_feature_processing as nxf