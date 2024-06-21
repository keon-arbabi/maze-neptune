import sys, polars as pl, numpy as np
import matplotlib.pylab as plt, seaborn as sns

sys.path.append('/home/karbabi/projects/def-wainberg/karbabi/utils')
from single_cell import SingleCell
from utils import Timer, print_df, savefig, debug

debug(third_party=True)
data_dir = 'projects/def-wainberg/single-cell/Maze/Neptune'
working_dir = 'projects/def-wainberg/karbabi/maze-neptune'

# Cell type annotation with reference atlas ####################################

# Load and QC the NEPTUNE human kidney single cell 
# Output from `preprocessing.R`
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
        .with_columns_obs(
            pl.col.APOL_Alleles.cast(pl.String)
                .str.count_matches(r'[1-9]')
                .fill_null(strategy='mean')
                .alias('APOL_dosage'))\
        .with_columns_obs(
            pl.col.Age.fill_null(strategy='mean'),
            Additive_Genotype=pl.col.APOL_dosage,
            Dominant_Genotype=pl.col.APOL_dosage>0,
            Recessive_Genotype=pl.col.APOL_dosage==2)\
        .drop_obs(['orig.ident', 'APOL_Allele_Number', 
                   'APOL_Alleles_right', 'N264K'])\
        .qc(cell_type_confidence_column=None, 
            doublet_column='IsDoublet',
            allow_float=True)
    sc_query_raw = sc_query.copy()

'''
Starting with 301,859 cells.
Filtering to cells with ≤10.0% mitochondrial counts...
301,859 cells remain after filtering to cells with ≤10.0% mitochondrial counts.
Filtering to cells with ≥100 genes detected (with nonzero count)...
301,859 cells remain after filtering to cells with ≥100 genes detected.
Removing doublets...
288,313 cells remain after removing doublets.
'''

# Load and QC Lake 2023 single cell atlas of healthy and injured human kidney
# doi.org/10.1038/s41586-023-05769-3
with Timer ('[Lake] Loading and QCing single cell data'):
    sc_ref = SingleCell(
        f'{data_dir}/Lake_2023_integrated_human_kidney.h5ad',
        num_threads=None)\
        .set_var_names('feature_name')\
        .cast_var({'feature_name': pl.String})\
        .rename_obs({'subclass.l1': 'subclass_l1'})\
        .cast_obs({'suspension_type': pl.String, 'subclass_l1': pl.String, 
                   'subclass.l3': pl.String})\
        .with_columns_obs(
            pl.coalesce(pl.when(pl.col('subclass_l1').eq('PT')).then(None)
                        .otherwise(pl.col('subclass_l1')),
                        pl.col('subclass.l3')).alias('subclass_combined'))\
        .qc(cell_type_confidence_column=None,
            doublet_column=None, 
            max_mito_fraction=0.15,
            allow_float=True)    
    sc_ref_raw = sc_ref.copy()

'''
Starting with 304,652 cells.
Filtering to cells with ≤15.0% mitochondrial counts...
224,577 cells remain after filtering to cells with ≤15.0% mitochondrial counts.
Filtering to cells with ≥100 genes detected (with nonzero count)...
224,577 cells remain after filtering to cells with ≥100 genes detected.
'''

# Find highly variable genes, normalize expression, then run PCA
with Timer('Highly variable genes'):
    sc_query, sc_ref = sc_query.hvg(sc_ref, allow_float=True)
with Timer('Normalize'):
    sc_query = sc_query.normalize(allow_float=True, num_threads=None)
    sc_ref = sc_ref.normalize(allow_float=True, num_threads=None)
with Timer('PCA'):
    sc_query, sc_ref = sc_query.PCA(sc_ref, verbose=True)

# Harmonize the principal components between the two datasets with Harmony:
# github.com/slowkow/harmonypy
with Timer('Harmony'):
    sc_query, sc_ref = sc_query.harmonize(
        sc_ref, pytorch=True, num_threads=None)
    
# Transfer cell-type labels from Lake et al. to Neptune
with Timer('Label transfer'):
    sc_query = sc_query\
    .label_transfer_from(
        sc_ref, cell_type_column='subclass_l1',
        cell_type_confidence_column='subclass_l1_confidence')\
    .label_transfer_from(
        sc_ref, cell_type_column='subclass_combined',
        cell_type_confidence_column='subclass_combined_confidence')\
    .with_columns_obs(
        passed_subclass_combined=pl.col.subclass_combined_confidence.ge(0.8),
        passed_subclass_l1=pl.col.subclass_l1_confidence.ge(0.8))
    
with Timer('Label transfer stats'):
    total_cells = sc_query.obs['passed_subclass_combined'].sum()
    print(f'{total_cells: ,} cells remain after filtering to cells with ≥80% '
        f'cell-type confidence.')
    
    # Mean confidence per cell type 
    print_df(sc_query.obs.group_by('subclass_combined')
            .agg(pl.col('subclass_combined_confidence').mean())
            .sort('subclass_combined_confidence'))
    
    # Barplot of percent change in cells filtered by cell type confidence
    # for each APOL dosage
    df = sc_query.obs.group_by([
        'APOL_dosage', 'subclass_combined', 'passed_subclass_combined'])\
        .count()\
        .filter(pl.col.passed_subclass_combined.is_not_null())\
        .pivot(values='count', 
               index=['APOL_dosage', 'subclass_combined'], 
               columns='passed_subclass_combined')\
        .with_columns(((pl.col.false / (pl.col.true + pl.col.false)) * 100)
                    .alias('percent_false'))
    sns.barplot(data=df, x='subclass_combined', y='percent_false', 
                order=df.sort(
                    'percent_false', descending=True)['subclass_combined'],
                hue='APOL_dosage')
    plt.ylabel('Percent change in cells failing annotation confidence')
    plt.xticks(rotation=90)
    savefig(f'{working_dir}/figures/sc_query_dropped_cells_barplot.png')

    # ECDF plot of cell type confidence 
    sns.ecdfplot(data=sc_query.obs, x='subclass_combined_confidence')
    kde = sns.kdeplot(data=sc_query.obs['subclass_combined_confidence'],
                      cumulative=True)
    plt.axvline(x=0.8, linestyle='--', color='r')
    plt.axvline(x=0.9, linestyle='--', color='b')
    x_kde, y_kde = kde.get_lines()[0].get_data()
    plt.text(0.8, y_kde[np.abs(x_kde - 0.8).argmin()],
            f'{y_kde[np.abs(x_kde - 0.8).argmin()]:.2f}', color='r')
    plt.text(0.9, y_kde[np.abs(x_kde - 0.9).argmin()],
            f'{y_kde[np.abs(x_kde - 0.9).argmin()]:.2f}', color='b')
    savefig(f'{working_dir}/figures/sc_query_subclass_combined_ecdf.png')

'''
200,508 cells remain after filtering to cells with ≥80% cell-type confidence.

subclass_combined  subclass_combined_confidence 
PapE               0.505556                     
NEU                0.605682                     
cycPT              0.68056                      
dPT                0.686882                     
PT-S1/2            0.704228                     
ATL                0.7875                       
PT-S3              0.788915                     
DTL                0.84083                      
aPT                0.853119                     
CNT                0.859353                     
VSM/P              0.900235                     
PC                 0.916042                     
FIB                0.924451                     
DCT                0.930926                     
TAL                0.946639                     
PEC                0.946858                     
IC                 0.950194                     
IMM                0.968465                     
EC                 0.974136                     
POD                0.977655  
'''

# Generate and plots UMAP
with Timer('UMAP plot'):
    sc_query = sc_query\
        .with_columns_obs(passed_tmp=pl.col.passed_QC & 
                          pl.col.passed_subclass_combined)\
        .embed(QC_column='passed_tmp', embedding_key='PaCMAP_subclass_combined',
               num_threads=24)
    sc_query.plot_embedding(
        'subclass_combined', 
        f'{working_dir}/figures/sc_query_subclass_combined_umap.png',
        embedding_key='PaCMAP_subclass_combined',
        cells_to_plot_column='passed_tmp',
        label=True, label_kwargs={'size': 6},
        legend=True, legend_kwargs={'fontsize': 'x-small', 'ncols': 1})
    sc_query = sc_query\
        .with_columns_obs(passed_tmp=pl.col.passed_QC & 
                          pl.col.passed_subclass_l1)\
        .embed(QC_column='passed_tmp', embedding_key='PaCMAP_subclass_l1',
               num_threads=24)
    sc_query.plot_embedding(
        'subclass_l1', 
        f'{working_dir}/figures/sc_query_subclass_l1_umap.png',
        embedding_key='PaCMAP_subclass_l1',
        cells_to_plot_column='passed_tmp',
        label=True, label_kwargs={'size': 6},
        legend=True, legend_kwargs={'fontsize': 'x-small', 'ncols': 1})
    sc_query = sc_query.drop_obs('passed_tmp')
    
# Save labelled single cell data 
with Timer('[Neptune] Saving single cell'):
    sc_query.X = sc_query_raw.X
    sc_query.save(f'{data_dir}/neptune_10x_labelled.h5ad', overwrite=True)
    
    sc_query_shareable = sc_query\
        .filter_obs(pl.col.passed_QC & pl.col.passed_subclass_combined)\
        .drop_obs(['passed_subclass_combined', 'passed_subclass_l1'])
    del sc_query_shareable.obsm['PaCMAP_subclass_l1']
    sc_query_shareable.save(f'{data_dir}/neptune_10x_labelled_shareable.h5ad',
                            overwrite=True)

    sc_ref.X = sc_ref_raw.X
    sc_ref.save(f'{data_dir}/Lake_2023_integrated_human_kidney_processed.h5ad', 
                overwrite=True)