import os, sys
import polars as pl, polars.selectors as cs
import pandas as pd, numpy as np
import matplotlib.pylab as plt, seaborn as sns

sys.path.append('/home/karbabi/projects/def-wainberg/karbabi/utils')
from single_cell import SingleCell, Pseudobulk, DE
from utils import Timer, print_df, savefig, get_coding_genes, debug
from ryp import r, to_r

debug(third_party=True)
data_dir = 'projects/def-wainberg/single-cell/Maze/Neptune'
working_dir = 'projects/def-wainberg/karbabi/maze-neptune'

# Cell type annotation with reference atlas ####################################
# Load the NEPTUNE human kidney single cell 
with Timer ('[Neptune] Loading and QCing single cell data'):
    additional_genotypes = pl.read_csv(
        f'{data_dir}/additional_neptune_genotypes.csv',
        null_values='N/A')
    sc_query = SingleCell(
        f'{data_dir}/neptune_10x.h5ad', num_threads=None)\
        .cast_obs({'ID': pl.String})\
        .join_obs(additional_genotypes, on='ID')\
        .with_columns_obs(
            pl.coalesce('APOL_Alleles', 'APOL_Alleles_right')
                .alias('APOL_Alleles'))\
        .drop_obs(['orig.ident', 'APOL_Allele_Number', 
                   'APOL_Alleles_right', 'N264K'])\
        .qc(cell_type_confidence_column=None, 
            doublet_column='IsDoublet',
            max_mito_fraction=0.05,
            min_genes=200,
            allow_float=True)
    sc_query_raw = sc_query 
    
# Load the Lake 2023 atlas of healthy and injured human kidney
# doi.org/10.1038/s41586-023-05769-3
with Timer ('[Lake] Loading and QCing single cell data'):
    sc_ref = SingleCell(
        f'{data_dir}/Lake_2023_integrated_human_kidney.h5ad',
        num_threads=None)\
        .set_var_names('feature_name')\
        .cast_var({'feature_name': pl.String})\
        .cast_obs({'suspension_type': pl.String, 'subclass.l1': pl.String, 
                   'subclass.l3': pl.String})\
        .with_columns_obs(
            pl.coalesce(
                pl.when(pl.col('subclass.l1').eq('PT')).then(None)
                    .otherwise(pl.col('subclass.l1')),
                pl.col('subclass.l3')).alias('subclass_combined'))\
        .qc(cell_type_confidence_column=None,
            doublet_column=None, 
            max_mito_fraction=0.05,
            min_genes=200,
            allow_float=True)

# Find highly variable genes, normalize expression, then run PCA
with Timer('Highly variable genes'):
    sc_query, sc_ref = sc_query.hvg(sc_ref, min_cells=3, allow_float=True)

with Timer('Normalize'):
    sc_query = sc_query.normalize(allow_float=True)
    sc_ref = sc_ref.normalize(allow_float=True)
    
with Timer('PCA'):
    sc_query, sc_ref = sc_query.PCA(sc_ref, verbose=True)

# Plot PC1 vs PC2
with Timer('Plot PCA'):
    plt.scatter(
        sc_query.obsm['PCs'][:, 0], sc_query.obsm['PCs'][:, 1],
        label='Neptune', s=1, alpha=0.05, rasterized=True)
    plt.scatter(
        sc_ref.obsm['PCs'][:, 0], sc_ref.obsm['PCs'][:, 1],
        label='Lake Atlas', s=1, alpha=0.05, rasterized=True)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    savefig(f'{working_dir}/figures/compare_pcs.png')

# Harmonize the principal components between the two datasets with Harmony:
# github.com/slowkow/harmonypy
# Note: the order of the two datasets doesn't matter   
with Timer('Harmony'):
    sc_query, sc_ref = sc_query.harmonize(
        sc_ref,  pytorch=True, num_threads=None)

# Generate new PCA plots after harmonization
with Timer('Plot PCA'):
    plt.scatter(
        sc_query.obsm['Harmony_PCs'][:, 0], 
        sc_query.obsm['Harmony_PCs'][:, 1],
        label='Neptune', s=1, alpha=0.05, rasterized=True)
    plt.scatter(
        sc_ref.obsm['Harmony_PCs'][:, 0],
        sc_ref.obsm['Harmony_PCs'][:, 1],
        label='Lake Atlas', s=1, alpha=0.05, rasterized=True)
    plt.xlabel('Harmony PC1')
    plt.ylabel('Harmony PC2')
    savefig(f'{working_dir}/figures/compare_pcs_harmony.png')
    
# Transfer cell-type labels from Lake et al. to Neptune
with Timer('Label transfer'):
    sc_query = sc_query\
    .label_transfer_from(
        sc_ref, cell_type_column='subclass_combined',
        cell_type_confidence_column='subclass_combined_confidence',
        min_cell_type_confidence=0.80,
        num_neighbors=20,
        num_index_neighbors=50)
    print_df(sc_query.obs.group_by('subclass_combined')
          .agg(pl.col('subclass_combined_confidence').mean())
          .sort('subclass_combined_confidence'))
    sns.ecdfplot(data=sc_query.obs, x="subclass_combined_confidence")
    savefig(f'{working_dir}/figures/sc_query_subclass_combined_ecdf.png')

# Generate and plots UMAP
with Timer('UMAP plot'):
    sc_query = sc_query.UMAP(seed=None, num_threads=24)
    sc_query.plot_UMAP(
        'subclass_combined', 
        f'{working_dir}/figures/sc_query_subclass_combined_umap2.png',
        label=True, label_kwargs={'size': 6},
        legend=True, legend_kwargs={'fontsize': 'x-small', 'ncols': 1})
    sc_query.plot_UMAP(
        'subclass_combined_confidence',
        f'{working_dir}/figures/sc_query_subclass_combined_confidence_umap.png')
    sc_query.plot_UMAP(
        'ID', f'{working_dir}/figures/sc_query_sample_umap.png', legend=False)
    
# Save labelled single cell data 
with Timer('[Neptune] Saving single cell'):
    sc_query.X = sc_query_raw.X
    sc_query.obsm['X_umap'] = sc_query.obsm['UMAP'] # for CellXGene
    sc_query.save(f'{data_dir}/neptune_10x_labelled.h5ad', overwrite=True)

# Peudobulk differential expression testing ####################################
# Pseudobulk the data 
with Timer('[Neptune] Pseuduobulking'):
    sc_query = SingleCell(f'{data_dir}/neptune_10x_labelled.h5ad')
    pb = sc_query\
        .cast_obs({'ID': pl.String, 'subclass_combined': pl.String})\
        .pseudobulk(ID_column='ID', cell_type_column='subclass_combined')\
        .filter_var(pl.col._index.is_in(get_coding_genes()['gene']))\
        .with_columns_obs(
            pl.col.APOL_Alleles.cast(pl.String)
                .str.count_matches(r'[1-9]')
                .fill_null(strategy='mean')
                .alias('APOL_dosage'))\
        .with_columns_obs(
            pl.col.Age.fill_null(strategy='mean'),
            Additive_Genotype=pl.col.APOL_dosage,
            Dominant_Genotype=pl.col.APOL_dosage>0,
            Recessive_Genotype=pl.col.APOL_dosage==2)
    drop_cell_types = []
    drop_cell_types.extend([
        cell_type for cell_type, (_, obs, _) in pb.items() 
        if obs.shape[0] < 30])
    pb = pb.drop_cell_types(drop_cell_types)
    if not os.path.exists(f'{data_dir}/subclass_combined_pseudobulk'):
        pb.save(f'{data_dir}/subclass_combined_pseudobulk')

# Limma-voom differential expression 
pb = Pseudobulk(f'{data_dir}/subclass_combined_pseudobulk')

label_column = [
    'Additive_Genotype', 'Dominant_Genotype', 'Recessive_Genotype']
case_control_column = {
    'Additive_Genotype': None, 'Dominant_Genotype': 'Dominant_Genotype', 
    'Recessive_Genotype': 'Recessive_Genotype'}
case_control = {
    'Additive_Genotype': False, 'Dominant_Genotype': True, 
    'Recessive_Genotype': True}

for label in label_column:
    with Timer(f'[Neptune] Differential expression for {label}'):
        de = pb\
            .qc(case_control_column=case_control_column[label], 
                custom_filter=pl.col(label).is_not_null())\
            .DE(label_column=label, 
                case_control=case_control[label],
                covariate_columns=['eGFR_Bx', 'Age', 'Sex'],
                include_library_size_as_covariate=True,
                include_num_cells_as_covariate=True,
                verbose=False)
        de.plot_voom(f'{working_dir}/figures/voom/{label}', 
                     overwrite=True, PNG=True)
        de.save(f'{working_dir}/output/DE_{label}', overwrite=True)   
        de.table.write_csv(f'{working_dir}/output/DE_{label}/table.csv')
        print(label)
        print_df(de.get_num_hits(threshold=0.1).sort('cell_type'))

''' 
Additive_Genotype
 cell_type  num_hits 
 EC         474      
 FIB        262      
 IMM        92       
 PC         5        
 aPT        3       
Dominant_Genotype
 cell_type  num_hits 
 EC         403      
 FIB        140      
 IMM        25       
 TAL        3        
 aPT        1     
Recessive_Genotype  
 cell_type  num_hits 
 FIB        59       
 IMM        23       
 aPT        8      
'''

# Plotting #####################################################################

de_a = DE(f'{working_dir}/output/DE_Additive_Genotype').table
de_d = DE(f'{working_dir}/output/DE_Dominant_Genotype').table
de_r = DE(f'{working_dir}/output/DE_Dominant_Genotype').table

# plotting to compare models 
plt.figure(figsize=(8, 5.5))
plt.subplot(2, 1, 1)
sns.kdeplot(np.log(de_a['SE']), 
            fill=True, color='blue', label='Additive model')
sns.kdeplot(np.log(de_d['SE']), 
            fill=True, color='green', label='Dominant model')
sns.kdeplot(np.log(de_r['SE']), 
            fill=True, color='red', label='Recessive model')
plt.title('Model Residual Error')
plt.xlabel('log(Genewise Error)')
plt.ylabel('Density')
plt.legend()

plt.subplot(2, 1, 2) 
sns.histplot(de_a['P'], 
             color='blue', fill=False, label='Additive model')
sns.histplot(de_d['P'], 
             color='green', fill=False, label='Dominant model')
sns.histplot(de_r['P'], 
             color='red', fill=False, label='Recessive model')
plt.title('Genewise P')
plt.xlabel('P')
plt.ylabel('Frequency')
plt.legend()
plt.tight_layout()
savefig(f'{working_dir}/figures/model_comparisons.png', dpi=300)
plt.clf()
