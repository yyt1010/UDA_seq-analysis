import matplotlib.pyplot as plt
import scanpy as sc
import dill
import os
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *
from scenicplus.cistromes import *
from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS
from scenicplus.enhancer_to_gene import export_to_UCSC_interact 
from scenicplus.utils import format_egrns, export_eRegulons
from scenicplus.eregulon_enrichment import *
from scenicplus.TF_to_gene import *
from scenicplus.dimensionality_reduction import *
from scenicplus.RSS import *
from scenicplus.diff_features import *
from scenicplus.loom import *
from typing import Dict, List, Mapping, Optional, Sequence
import os
import dill
import time

plt.rcParams["pdf.fonttype"]=42

tmp_dir = '/xtdisk/jiangl_group/huangzh/tmp5'
work_dir= '/p300s/jiangl_group/huangzh/Post_index/Kidney/Merge/T3/RNA/results/Scissor1/2_proteinuria_3g/All_cell_0.8/POD/SCENICplus4/SCENICplus'


scplus_obj=dill.load(open(os.path.join(work_dir, 'scenicplus2/scplus_obj.pkl'), 'rb'))

from scenicplus.RSS import *

calculate_TF_eGRN_correlation=True
#var="ACC_Cluster"
variable=["ACC_Cluster"]
if True:
    log.info('Calculating TF-eGRNs AUC correlation')
    for var in variable:
        generate_pseudobulks(scplus_obj, 
                             variable = var,
                             auc_key = 'eRegulon_AUC',
                             signature_key = 'Gene_based',
                             nr_cells = 5,
                             nr_pseudobulks = 100,
                             seed=555)
        generate_pseudobulks(scplus_obj, 
                                     variable = var,
                                     auc_key = 'eRegulon_AUC',
                                     signature_key = 'Region_based',
                                     nr_cells = 5,
                                     nr_pseudobulks = 100,
                                     seed=555)
        TF_cistrome_correlation(scplus_obj,
                            variable = var, 
                            auc_key = 'eRegulon_AUC',
                            signature_key = 'Gene_based',
                            out_key = var+'_eGRN_gene_based')
        TF_cistrome_correlation(scplus_obj,
                                    variable = var, 
                                    auc_key = 'eRegulon_AUC',
                                    signature_key = 'Region_based',
                                    out_key = var+'_eGRN_region_based')
log.info('Calculating eRSS')
for var in variable:
    regulon_specificity_scores(scplus_obj, 
                         var,
                         signature_keys=['Gene_based'],
                         out_key_suffix='_gene_based',
                         scale=False)
    regulon_specificity_scores(scplus_obj, 
                         var,
                         signature_keys=['Region_based'],
                         out_key_suffix='_region_based',
                         scale=False)
regulon_specificity_scores(
        scplus_obj,
        variable = 'ACC_Cluster',
        auc_key = 'eRegulon_AUC_filtered',
        signature_keys = ['Region_based'],
        selected_regulons = [x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
        out_key_suffix = '_filtered')

flat_list = lambda t: [item for sublist in t for item in sublist]
selected_markers = list(set(flat_list(
    [scplus_obj.uns['RSS']['ACC_Cluster_filtered'].loc[celltype].sort_values(ascending = False).head(10).index.to_list()
    for celltype in scplus_obj.uns['RSS']['ACC_Cluster_filtered'].index])))
plot_rss(scplus_obj, 'ACC_Cluster_filtered', num_columns=2, top_n=10, figsize = (5, 10))
plt.show()
plt.savefig("step7_rss.pdf")
from scenicplus.plotting.correlation_plot import *

region_intersetc_data, Z = jaccard_heatmap(
        scplus_obj,
        method = 'intersect',
        gene_or_region_based = 'Region_based',
        use_plotly = False,
        selected_regulons = selected_markers,
        signature_key = 'eRegulon_signatures_filtered',
        figsize = (10, 10), return_data = True, vmax = 0.5, cmap = 'plasma')

from scenicplus.plotting.dotplot import heatmap_dotplot
a=heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'], #specify what to plot as dot sizes, target region enrichment in this case
        color_matrix = scplus_obj.to_df('EXP'), #specify  what to plot as colors, TF expression in this case
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'ACC_Cluster',
        subset_eRegulons = scplus_obj.uns['selected_eRegulon']['Gene_based'],
        index_order = ['cluster_0','cluster_1','cluster_2','cluster_3','cluster_4','cluster_5'],
        orientation = 'vertical')
print(a)
plt.show()
plt.savefig("step7_heatmap.pdf")

from pycisTopic.diff_features import find_highly_variable_features
hvr = find_highly_variable_features(scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata']['Region']))], n_top_features=250, plot = False)
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata']['Gene']))].T, n_top_features=250, plot = False)

print("+++")
print(hvr)
print(hvg)
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
print(scplus_obj.uns['eRegulon_metadata_filtered'])
nx_tables = create_nx_tables(
    scplus_obj = scplus_obj,
    eRegulon_metadata_key ='eRegulon_metadata_filtered',
    subset_eRegulons = ['PLAG1','ZNF398',"E2F3"],
    subset_regions = hvr,
    subset_genes = hvg,
    add_differential_gene_expression = True,
    add_differential_region_accessibility = True,
    differential_variable = ['ACC_Cluster'])

print(nx_tables)
G, pos, edge_tables, node_tables = create_nx_graph(nx_tables,
                   use_edge_tables = ['TF2R','R2G'],
                   color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {'PLAG1': 'Orange', 'ZNF398': 'Purple', 'E2F3': 'Red'}},
                                    'R2G': {'variable' : 'R2G_rho', 'continuous_color' : 'viridis', 'v_min': -1, 'v_max': 1}},
                   transparency_edge_by =  {'R2G': {'variable' : 'R2G_importance', 'min_alpha': 0.1, 'v_min': 0}},
                   width_edge_by = {'R2G': {'variable' : 'R2G_importance', 'max_size' :  1.5, 'min_size' : 1}},
                   color_node_by = {'TF': {'variable': 'TF', 'category_color' : {'PLAG1': 'Orange', 'ZNF398': 'Purple', 'E2F3': 'Red'}},
                                    'Gene': {'variable': 'ACC_Cluster_Log2FC_cluster_1', 'continuous_color' : 'bwr'},
                                    'Region': {'variable': 'ACC_Cluster_Log2FC_cluster_1', 'continuous_color' : 'viridis'}},
                   transparency_node_by =  {'Region': {'variable' : 'ACC_Cluster_Log2FC_cluster_1', 'min_alpha': 0.1},
                                    'Gene': {'variable' : 'ACC_Cluster_Log2FC_cluster_1', 'min_alpha': 0.1}},
                   size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 30},
                                    'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
                                    'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
                   shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
                                    'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
                   label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 20.0},
                                    'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 10.0},
                                    'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
                   layout='kamada_kawai_layout',
                   scale_position_by=250)
plt.figure(figsize=(10,10))
plot_networkx(G, pos)

plt.show()
plt.savefig("step7_new.pdf")

print("PTGIS" in hvg)
