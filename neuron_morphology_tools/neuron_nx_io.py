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
    remove_starter_branches = True,
    divide_into_limbs = True,

    label_name = None,#"axon_label",
    graph_label = None,
    
    distance_threshold = 50_000,
    
    #output file features
    folder = "./",
    filename = None,
    description = None,
    verbose = False,
    
    
    ):
    """
    To process a neuron object before output
    in dictionary format to be used by a GNN
    """
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
                feature_matrix_dtype = "float",

            )

            limb_info["label_name"] = label_name

            curr_filename = f"{filename}_limb_{j}_starting_branch_{most_starting_branch}"

            output_path = str((Path(folder)/Path(curr_filename)).absolute())

            ret_filepath = su.compressed_pickle(
                limb_info,
                output_path,
                return_filepath=True,
                verbose = verbose)

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

            )

        G_info["label"] = graph_label

        curr_filename = f"{filename}"
        output_path = str((Path(folder)/Path(curr_filename)).absolute())

        ret_filepath = su.compressed_pickle(
                G_info,
                output_path,
                return_filepath=True,
                verbose = verbose)

        filepaths.append(ret_filepath)

    if verbose:
        print(f"\n\n---Total time = {time.time() - st}")
    return filepaths
    
    

import system_utils as su
import networkx_utils as xu
import pandas as pd


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

    df = pd.DataFrame(adj_feature_dict["feature_matrix"])
    df.columns = adj_feature_dict["features"]
    df["node"] = adj_feature_dict["nodelist"]
    
    G = xu.G_from_adjacency_matrix(
    matrix = adj_feature_dict["adjacency"],
    nodelist = adj_feature_dict["nodelist"],
    plot = True,
    )
    
    G = xu.set_node_attributes_from_df(G,df,index_name="node")
    return G

import neuron_nx_io as nxio