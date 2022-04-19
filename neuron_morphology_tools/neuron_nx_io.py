# ----------- For outputing for use in GNN
import system_utils as su
import neuron_nx_feature_processing as nxf
import neuron_nx_utils as nxu
import networkx as nx
import networkx_utils as xu
from pathlib import Path
import time
def export_GNN_info_dict(
    G,
    features_to_output,
    axon_dendrite = None,
    remove_starter_branches = True,
    distance_threshold = 100_000,
    divide_into_limbs = True,
    

    label_name = None,#"axon_label",
    graph_label = None,

    feature_matrix_dtype = "float",
    
    #output file features
    folder = "./",
    filename = None,
    description = None,
    verbose = False,
    
    return_filepaths = False,
    
    
    
    ):
    """
    To process a neuron object before output
    in dictionary format to be used by a GNN
    """
    print(f"return_filepaths =- {return_filepaths}")
    
    
    st = time.time()
    if label_name is None:
        if verbose:
            print(f"*** Warning label_name is None")
            
    if graph_label is None:
        if verbose:
            print(f"*** Warning graph_label is None")
    
    
    G_dict = xu.graph_attr_dict(G)
    
    if verbose:
        print(f"G_dict = {G_dict}")
    
    filename = nxu.name_from_G(G)
    
    if description is not None:
        filename += f"_{description}"
        
        
        
    filepaths = []

    # ----------- Does a lot of the preprocessing before outputting----------
    if axon_dendrite is not None:
        if verbose:
            print(f"Filtering for {axon_dendrite}")
        G = nxu.axon_dendrite_subgraph(
            G,
            compartment=axon_dendrite,
            include_soma = True,
            verbose = verbose,
            )
        
        
    if remove_starter_branches:
        G_filt = nxu.remove_small_starter_branches(
            G,
            verbose = verbose,
            maintain_skeleton_connectivity = True)
    else:
        G_filt = G
        

    if distance_threshold is not None:
        G_dist_filt = nxu.nodes_within_distance_upstream_from_soma(
            G_filt,
            verbose = verbose,
            distance_threshold = distance_threshold,
            return_subgraph = True,
        )
    else:
        G_dist_filt = G_filt


    G_with_feats = nxf.filter_G_features(
                G_dist_filt,
                features=features_to_output,
                inplace = False,
                verbose = verbose,
            )

    # ----------- Dividing up and outputting the files -----------
    if divide_into_limbs:
        limb_graphs_for_axon = nxu.limb_graphs_from_soma_connected_nodes(G_with_feats)
        for j,G_limb in enumerate(limb_graphs_for_axon):
            if verbose:
                print(f"Outputing limb {j}---")


            most_starting_branch = xu.most_upstream_node(G_limb)
            #converting to non-directional
            G_limb = nx.Graph(G_limb)

            G_limb = nxu.limb_branch_subgraph(G_limb)
            limb_info = xu.adjacency_feature_info(
                G = G_limb,
                return_df_for_feature_matrix = False,
                feature_matrix_dtype = feature_matrix_dtype,
                dense_adjacency=True

            )

            limb_info["label_name"] = label_name
            limb_info["graph_label"] = graph_label

            curr_filename = f"{filename}_limb_{j}_starting_branch_{most_starting_branch}"

            output_path = str((Path(folder)/Path(curr_filename)).absolute())

            if return_filepaths:
                ret_filepath = su.compressed_pickle(
                    limb_info,
                    output_path,
                    return_filepath=True,
                    verbose = verbose)
            else:
                ret_filepath = limb_info

            filepaths.append(ret_filepath)
    else:
        """
        Pseudocode: 
        1) Remove the soma from the graph
        2) Attach the label of the graph
        """
        if verbose:
            print()

        G_no_soma = nxu.soma_filter_by_complete_graph(G_with_feats,plot=False)
        G_no_soma = nx.Graph(G_no_soma)

        G_info = xu.adjacency_feature_info(
                G = G_no_soma,
                return_df_for_feature_matrix = False,
                feature_matrix_dtype = "float",
                dense_adjacency=True

            )

        G_info["label_name"] = label_name
        G_info["graph_label"] = graph_label

        curr_filename = f"{filename}"
        output_path = str((Path(folder)/Path(curr_filename)).absolute())

        if return_filepaths:
            ret_filepath = su.compressed_pickle(
                    G_info,
                    output_path,
                    return_filepath=True,
                    verbose = verbose)
        else:
            ret_filepath = G_info

        filepaths.append(ret_filepath)

    if verbose:
        print(f"\n\n---Total time = {time.time() - st}")
        
        
    return filepaths
    
    

import system_utils as su
import networkx_utils as xu
import pandas as pd

import pandas as pd
import numpy_utils as nu
import pandas_utils as pu
def feature_df_from_gnn_info(
    gnn_info,
    return_data_labels_split = True,
    inf_fill_value = 1000):
    df = pd.DataFrame(gnn_info["feature_matrix"])
    df.columns = gnn_info["features"]

    df=df.replace([np.inf],inf_fill_value)

    label_name = gnn_info["label_name"]

    if return_data_labels_split:
        if label_name is not None:
            label_name = list(nu.convert_to_array_like(label_name))
            y = df[list(label_name)].to_numpy()
            x = pu.delete_columns(df,label_name)
        else:
            y = None

        return x.to_numpy(),y
    else:
        return df

def feature_df_from_adj_feature_dict(adj_feature_dict):
    if "pandas" not in str(type(adj_feature_dict["feature_matrix"])):
        df = pd.DataFrame(adj_feature_dict["feature_matrix"])
        df.columns = adj_feature_dict["features"]
        df["node"] = adj_feature_dict["nodelist"]
    else:
        df = adj_feature_dict["feature_matrix"]
        
    return df

def G_from_adj_feature_dict(
    adj_feature_dict=None,
    filepath = None,
    plot = False,
    verbose = False
    ):
    """
    Purpose: To recover the original graph
    stored in the adjacency dict information
    
    Ex: 
    import neuron_nx_io as nxio
    G_rec = nxio.G_from_adj_feature_dict(
        filepath = filepaths[1],
        plot = True,
        verbose = True
        )
    """
    
    if adj_feature_dict is None:
        if verbose:
            print(f"Reading from {filepath}")
        adj_feature_dict = su.decompress_pickle(filepath)

    df = nxio.feature_df_from_adj_feature_dict(adj_feature_dict)
        
    
    G = xu.G_from_adjacency_matrix(
    matrix = adj_feature_dict["adjacency"],
    nodelist = adj_feature_dict["nodelist"],
    plot = True,
    )
    
    G = xu.set_node_attributes_from_df(G,df,index_name="node")
    
    if verbose:
        print("label_name,graph_label = ",(adj_feature_dict["label_name"],adj_feature_dict["graph_label"]))
    return G



# ----------------- exporting different types of graph attributes for GNNs -----------
def GNN_info_axon_vs_dendrite(
    G,
    distance_threshold = 100_000,
    
    remove_starter_branches = True,
    divide_into_limbs = False,
    label_name = "axon_label",
    graph_label = None,
    
    return_filepaths = False,
    folder = "./Axon_vs_Dendrite/",
    description = "ax_vs_dendr",
    verbose = False,
    
    **kwargs
    ):
    
    
    features_to_output = [
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
        "axon_label"
        ]
    
    #print(f"return_filepaths =- {return_filepaths}")
    
    filepaths = nxio.export_GNN_info_dict(
            G,
            features_to_output=features_to_output,
            remove_starter_branches = remove_starter_branches,
            divide_into_limbs = divide_into_limbs,

            label_name = label_name,#"axon_label",
            graph_label = graph_label,

            distance_threshold = distance_threshold,

            #output file features
            folder = folder,
            description = description,
            verbose = verbose,
        
            return_filepaths = return_filepaths,
            **kwargs

            )
    
    
    return filepaths
    
    
def GNN_info_compartment_proof(
    G,
    distance_threshold = None,
    
    remove_starter_branches = True,
    divide_into_limbs = False,
    label_name = ("axon_label","basal_label","apical_label","dendrite_label"),
    
    graph_label = None,
    
    return_filepaths = False,
    folder = "./Compartments_Proof/",
    description = "compartment_proof",
    verbose = False,
    
    **kwargs
    ):
    
    
    features_to_output = [
        "mesh_volume",
        "n_spines",
        "total_spine_volume",
        "n_synapses_post",
        "n_synapses_pre",
        "n_synapses_head",
        "n_synapses_neck",
        #"parent_skeletal_angle",
        "skeletal_length",
        "skeleton_vector_upstream_theta",
        "skeleton_vector_upstream_phi",
        "skeleton_vector_downstream_theta",
        "skeleton_vector_downstream_phi",
        "width_upstream",
        "width_no_spine",
        "width_downstream",
        "axon_label",
        "dendrite_label",
        "basal_label",
        "apical_label",
        ]
    
    filepaths = nxio.export_GNN_info_dict(
    G,
    features_to_output=features_to_output,
    remove_starter_branches = remove_starter_branches,
    divide_into_limbs = divide_into_limbs,

    label_name = label_name,#"axon_label",
    graph_label = graph_label,
    
    distance_threshold = distance_threshold,
    
    #output file features
    folder = folder,
    description = description,
    verbose = verbose,
        
    return_filepaths = return_filepaths,
    
    )
    
    return filepaths


def GNN_info_merge_errors(
    G,
    distance_threshold = None,

    remove_starter_branches = True,
    divide_into_limbs = False,
    #label_name = "auto_proof_filter_label",
    label_name = ("merge_high_degree_branching_label",
        "merge_low_degree_branching_label",
        "merge_width_jump_up_axon_label",
        "merge_axon_on_dendrite_label",
        "merge_high_degree_branching_dendrite_label",
        "merge_width_jump_up_dendrite_label",
        "merge_double_back_dendrite_label",),
    graph_label = None,

    axon_dendrite = None,
    
    return_filepaths = False,
    folder = "./Merge_Errors/",
    description = "merge_errors",
    
    
    verbose = False,
    
    **kwargs
    ):
    
    features_to_output = [
        "mesh_volume",
        "n_spines",
        "total_spine_volume",
        "n_synapses_post",
        "n_synapses_pre",
        "n_synapses_head",
        "n_synapses_neck",
        #"parent_skeletal_angle",
        "skeletal_length",
        "skeleton_vector_upstream_theta",
        "skeleton_vector_upstream_phi",
        "skeleton_vector_downstream_theta",
        "skeleton_vector_downstream_phi",
        "width_upstream",
        "width_no_spine",
        "width_downstream",
        "min_dist_synapses_pre_downstream_clip",
        "min_dist_synapses_pre_upstream_clip",
        
    ]
    
    if type(label_name) == str:
        features_to_output.append(label_name)
    else:
        features_to_output += label_name
        # -------- labels ------------
#         "merge_high_degree_branching_label",
#         "merge_low_degree_branching_label",
#         "merge_width_jump_up_axon_label",
#         "merge_axon_on_dendrite_label",
#         "merge_high_degree_branching_dendrite_label",
#         "merge_width_jump_up_dendrite_label",
#         "merge_double_back_dendrite_label",
        
#         ]
    
    filepaths = nxio.export_GNN_info_dict(
        G,
        features_to_output=features_to_output,
        remove_starter_branches = remove_starter_branches,
        divide_into_limbs = divide_into_limbs,

        label_name = label_name,#"axon_label",
        graph_label = graph_label,

        distance_threshold = distance_threshold,

        #output file features
        folder = folder,
        description = description,
        verbose = verbose,

        axon_dendrite = axon_dendrite,
        return_filepaths = return_filepaths
        )
    
    return filepaths


def GNN_info_cell_type_fine(
    G,
    distance_threshold = None,

    remove_starter_branches = True,
    divide_into_limbs = False,
    label_name = None,
    graph_label = None,

    axon_dendrite = None,
    
    return_filepaths = False,
    folder = "./Cell_Type_Fine/",
    description = "cell_type_fine",
    
    
    verbose = False,
    
    **kwargs
    ):
    
    features_to_output = [
        "mesh_volume",
        "n_spines",
        "total_spine_volume",
        "n_synapses_post",
        "n_synapses_pre",
        "n_synapses_head",
        "n_synapses_neck",
        #"parent_skeletal_angle",
        "skeletal_length",
        "skeleton_vector_upstream_theta",
        "skeleton_vector_upstream_phi",
        "skeleton_vector_downstream_theta",
        "skeleton_vector_downstream_phi",
        "width_upstream",
        "width_no_spine",
        "width_downstream",
        "axon_label",
        "dendrite_label",
        "basal_label",
        "apical_label"
        ]
    
    filepaths = nxio.export_GNN_info_dict(
        G,
        features_to_output=features_to_output,
        remove_starter_branches = remove_starter_branches,
        divide_into_limbs = divide_into_limbs,

        label_name = label_name,#"axon_label",
        graph_label = graph_label,

        distance_threshold = distance_threshold,

        #output file features
        folder = folder,
        description = description,
        verbose = verbose,

        axon_dendrite = axon_dendrite,
        return_filepaths = return_filepaths
        )
    
    return filepaths

import neuron_nx_io as nxio