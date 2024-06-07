import os, sys
import polars as pl, polars.selectors as cs
import pandas as pd, numpy as np
import matplotlib.pylab as plt, seaborn as sns

sys.path.append('/home/karbabi/projects/def-wainberg/karbabi/utils')
from single_cell import SingleCell, Pseudobulk, DE
from utils import Timer, print_df, savefig, get_coding_genes, debug

debug(third_party=True)
data_dir = 'projects/def-wainberg/single-cell/Maze/Neptune'
working_dir = 'projects/def-wainberg/karbabi/maze-neptune'

# Load the NEPTUNE human kidney single cell 
with Timer ('[Neptune] Loading and QCing single cell data'):
    additional_genotypes = pl.read_csv(
        f'{data_dir}/additional_neptune_genotypes.csv',
        null_values='N/A')
    sc_query = SingleCell(f'{data_dir}/neptune_10x.h5ad')\
        .cast_obs({'ID': pl.String})\
        .join_obs(additional_genotypes, on='ID', validate='m:1')\
        .with_columns_obs(
            pl.col.ID.alias('Batch'),
            pl.coalesce('APOL_Alleles', 'APOL_Alleles_right')
                .alias('APOL_Alleles'))\
        .drop_obs(['orig.ident', 'APOL_Allele_Number', 
                   'APOL_Alleles_right', 'N264K'])\
        .qc(cell_type_confidence_column=None, 
            doublet_column='IsDoublet',
            max_mito_fraction=0.05,
            min_genes=200,
            allow_float=True)

# Load the Lake 2023 atlas of healthy and injured human kidney
# https://cellxgene.cziscience.com/collections/bcb61471-2a44-4d00-a0af-ff085512674c
with Timer ('[Lake] Loading and QCing single cell data'):
    sc_ref = SingleCell(f'{data_dir}/Lake_2023_integrated_human_kidney.h5ad')\
        .set_var_names('feature_name')\
        .cast_var({'feature_name': pl.String})\
        .cast_obs({'specimen': pl.String, 'subclass.l1': pl.String, 
                   'subclass.l3': pl.String})\
        .with_columns_obs(
            pl.col.specimen.alias('Batch'),
            pl.coalesce(
                pl.when(pl.col('subclass.l1').eq('PT')).then(None)
                    .otherwise(pl.col('subclass.l1')),
                pl.col('subclass.l3')).alias('subclass.combined'))\
        .qc(cell_type_confidence_column=None,
            doublet_column=None, 
            # custom_filter=pl.col.Batch.len().over('Batch') > 500,
            max_mito_fraction=0.05,
            min_genes=200,
            allow_float=True)

# Find highly variable genes
with Timer('Highly variable genes'):
    sc_query, sc_ref = sc_query.hvg(
        sc_ref, min_cells=None, allow_float=True)
    
# Run PCA
with Timer('PCA'):
    sc_query, sc_ref = sc_query.PCA(
        sc_ref, allow_float=True, verbose=True)

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
    savefig(f'{working_dir}/figures/PC1_vs_PC2.png')

# Harmonize the principal components between the two datasets with Harmony:
# github.com/slowkow/harmonypy
# Note: the order of the two datasets doesn't matter   
with Timer('Harmony'):
    sc_query, sc_ref = sc_query.harmonize(
        sc_ref, pytorch=True, num_threads=None)

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
    savefig(f'{working_dir}/figures/Harmony_PC1_vs_PC2.png')
    
# Transfer cell-type labels from Lake et al. to Neptune
with Timer('Label transfer'):
    sc_query = sc_query\
    .label_transfer_from(
        sc_ref, cell_type_column='subclass.combined',
        cell_type_confidence_column='subclass.combined_confidence',
        min_cell_type_confidence=0.80,
        num_neighbors=5,
        num_index_neighbors=30)\
    .label_transfer_from(
        sc_ref, cell_type_column='condition.l2',
        cell_type_confidence_column='condition.l2_confidence')
    print_df(sc_query.obs.group_by('subclass.combined')
          .agg(pl.col('subclass.combined_confidence').mean())
          .sort('subclass.combined_confidence'))

# Generate and plots UMAP
with Timer('UMAP plot'):
    sc_query = sc_query.UMAP(seed=None, num_threads=24)
    sc_query.plot_UMAP('subclass.combined', 
                       'figures/sc_query_subclass.combined_umap.png')
    sc_query.plot_UMAP('subclass.combined_confidence',
                       'figures/sc_query_subclass.combined_confidence_umap.png')

# Save labelled single cell data 
with Timer('[Neptune] Saving single cell'):
    sc_query.save(f'{data_dir}/neptune_10x_labelled.h5ad', overwrite=False)
    
# Pseudobulk the data 
with Timer('[Neptune] Pseuduobulking'):
    # sc_query = SingleCell(f'{data_dir}/neptune_10x_labelled.h5ad')
    pb = sc_query\
        .cast_obs({'ID': pl.String})\
        .pseudobulk(ID_column='ID', cell_type_column='subclass.combined')\
        .filter_var(pl.col._index.is_in(get_coding_genes()['gene']))\
        .with_columns_obs(
            pl.col.APOL_Alleles
                .str.count_matches(r'[1-9]')
                .fill_null(strategy='mean')
                .alias('APOL_dosage'))\
        .with_columns_obs(
            pl.col.Age.fill_null(strategy='mean'),
            Additive_Genotype=pl.col.APOL_dosage,
            Dominant_Genotype=pl.col.APOL_dosage>0,
            Recessive_Genotype=pl.col.APOL_dosage==2)
    # drop cell types missing a lot of cells across people 
    pb = pb.drop_cell_types(['NEU', 'ATL'])
    if not os.path.exists(f'{data_dir}/subclass.combined_pseudobulk'):
        pb.save(f'{data_dir}/subclass.combined_pseudobulk')

# Limma-voom differential expression 
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
        # pb = Pseudobulk(f'output/pseudobulk/neptune_pseudobulk')
        de = pb\
            .qc(case_control_column=case_control_column[label], 
                custom_filter=pl.col(label).is_not_null())\
            .DE(label_column=label, 
                case_control=case_control[label],
                covariate_columns=['eGFR_Bx', 'Age', 'Sex'],
                include_library_size_as_covariate=True,
                include_num_cells_as_covariate=True)
        de.plot_voom(save_to=f'figures/voom/{label}', overwrite=True)
        de.save(f'output/DE_{label}', overwrite=True)    
        print(label)
        print_df(de.get_num_hits(threshold=0.05).sort('cell_type'))

for obs in pb.iter_obs():
    print(obs.shape[0])


meta = sc_query.get_sample_covariates(ID_column='ID')
print_df(meta)
print_df(pl.DataFrame({
    "Column Name": meta.columns,
    "Data Type": [str(dtype) for dtype in meta.dtypes]
}))


pb.qc(case_control_column='Recessive_Genotype',
        min_nonzero_percent=80)

pb.obs['TAL']['APOL_Allele_Number']
pb.obs['TAL'].group_by('Recessive_Genotype').count()

# differential expression 
de_a = pb\
    .DE(DE_label_column='Additive_Genotype', 
        case_control=False,
        covariate_columns=['eGFR_Bx', 'Age', 'Sex', 'Cohort'],
        return_voom_weights=False,
        voom_plot_directory=None)
de_d = pb\
    .DE(DE_label_column='Dominant_Genotype', 
        covariate_columns=['eGFR_Bx', 'Age', 'Sex', 'Cohort'],
        return_voom_weights=False,
        voom_plot_directory=None)
de_r = pb\
    .DE(DE_label_column='Recessive_Genotype', 
        covariate_columns=['eGFR_Bx', 'Age', 'Sex', 'Cohort'],
        return_voom_weights=False,
        voom_plot_directory='figures/voom')
de_r.write_csv('results/differential-expression/DE_results_broad_recessive.csv')

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
plt.savefig("figures/model_comparisons.png", dpi=300)
plt.clf()

Pseudobulk.get_num_DE_hits(de_r, threshold=0.05)

'''
cell_type  num_hits 
DCT        755      
TAL        4        
POD        1782     
VSM/P      1        
PEC        1        
IC         3222     
PT         1682     
CNT        10  
'''

print_df(Pseudobulk.get_DE_hits(de_r, threshold=0.1, num_top_hits=50)\
        .filter(pl.col.cell_type=='POD')) 

def plot_volcano(df, logFC_threshold=(-1, 1), P_threshold=0.05,
                colors=('blue', 'grey', 'red'), alpha=1, size=6,
                plot_threshold_lines=True, label_top_genes=True, 
                num_top_genes=20, adjust=True,
                plot_as_pdf=False, plot_directory=None):
    import matplotlib.patches as mpatches
    from adjustText import adjust_text  

    df = df.to_pandas()
    sig_down = (df['logFC'] <= logFC_threshold[0]) & (df['P'] < P_threshold)
    sig_up = (df['logFC'] >= logFC_threshold[1]) & (df['P'] < P_threshold)
    df['color'] = colors[1]  
    df.loc[sig_down, 'color'] = colors[0]
    df.loc[sig_up, 'color'] = colors[2]
    df['logP'] = -np.log10(df['P'])  
    
    plt.figure(figsize=(5, 4.5))
    scatter = plt.scatter(df['logFC'], df['logP'], c=df['color'],
                          alpha=alpha, s=size, marker="o")
    if label_top_genes:
        top_genes = pd.concat([
            df[sig_down].assign(score=lambda x: x['logP']*abs(x['logFC']))\
                .nlargest(num_top_genes, 'score'),
            df[sig_up].assign(score=lambda x: x['logP']*abs(x['logFC']))\
                .nlargest(num_top_genes, 'score')])
        texts = [plt.text(row['logFC'], row['logP'], row['gene'], fontsize=8)
                 for _, row in top_genes.iterrows()]
        if adjust:
            adjust_text(texts, expand=(1.3, 1.5),
                        arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
    if plot_threshold_lines:
        plt.axhline(y=-np.log10(P_threshold), color='grey', linestyle='--')
        plt.axvline(x=logFC_threshold[0], color='grey', linestyle='--')
        plt.axvline(x=logFC_threshold[1], color='grey', linestyle='--')
    cell_type = df['cell_type'].unique()[0]
    plt.title(f"{cell_type}")
    plt.xlabel(r'$\log_2(Fold Change)$')
    plt.ylabel(r'$-\log_{10}(P-value)$')
    legend_patches = [
        mpatches.Patch(color=colors[0], label='significant down'),
        mpatches.Patch(color=colors[1], label='not significant'),
        mpatches.Patch(color=colors[2], label='significant up')
    ]
    plt.legend(handles=legend_patches, loc='best')
    if plot_directory is not None:
        plot_file = os.path.join(
            plot_directory,
            f'{cell_type.replace("/", "-")}.'
            f'{"pdf" if plot_as_pdf else "png"}')
        savefig(plot_file)

for cell_type in de_rg['cell_type'].unique():
    df = de_rg.filter(pl.col.cell_type==cell_type)
    plot_volcano(df, plot_directory='figures/volcano')

with Timer('[Neptune] UMAP'):
    sc_query = sc_query.UMAP(seed=None, num_threads=24)
    sc_query.plot_UMAP(None, 'figures/sc_query_umap.pdf')

with Timer('[Lake] UMAP'):
    sc_ref = sc_ref.UMAP(seed=None, num_threads=24)
    sc_ref.plot_UMAP('subclass.l1', 'figures/sc_ref_broad_umap.pdf', 
                     label=True, legend=False)
    sc_ref.plot_UMAP('subclass.full', 'figures/sc_ref_fine_umap.pdf')