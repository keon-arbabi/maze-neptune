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

# Peudobulk differential expression testing ####################################
# Pseudobulk
with Timer('[Neptune] Pseuduobulking'):
    sc_query = SingleCell(f'{data_dir}/neptune_10x_labelled.h5ad',
                          num_threads=None)\
        .cast_obs({'ID': pl.String, 'subclass_combined': pl.String})
    pb = sc_query.with_columns_obs(passed_tmp=pl.col.passed_QC & 
                                   pl.col.passed_subclass_combined)\
        .pseudobulk(ID_column='ID', cell_type_column='subclass_combined',
                    QC_column='passed_tmp') | \
        sc_query.with_columns_obs(passed_tmp=pl.col.passed_QC & 
                                  pl.col.passed_subclass_l1)\
            .filter_obs(pl.col.subclass_l1.eq('PT'))\
        .pseudobulk(ID_column='ID', cell_type_column='subclass_l1',
                    QC_column='passed_tmp')
    pb = pb.filter_var(pl.col._index.is_in(get_coding_genes()['gene']))
    drop_cell_types = []
    drop_cell_types.extend([
        cell_type for cell_type, (_, obs, _) in pb.items() 
        if obs.shape[0] < 30])
    pb = pb.drop_cell_types(drop_cell_types)

# Limma-voom 
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
                custom_filter=pl.col(label).is_not_null(),
                verbose=False)\
            .DE(label_column=label, 
                case_control=case_control[label],
                covariate_columns=['eGFR_Bx', 'Age', 'Sex'],
                include_library_size_as_covariate=True,
                include_num_cells_as_covariate=True,
                verbose=False)
        de.plot_voom(f'{working_dir}/figures/voom/{label}', 
                     overwrite=True, PNG=True)
        de.save(f'{working_dir}/output/DE_{label}', overwrite=True)   
        de.table.write_csv(f'{working_dir}/output/DE_{label}_table.csv')
        print(label)
        print_df(de.get_num_hits(threshold=0.1).sort('cell_type'))

''' 
Additive_Genotype
 cell_type  num_hits 
 EC         390      
 FIB        218      
 IMM        69       
 PC         1        
 aPT        21         
Dominant_Genotype
 cell_type  num_hits 
 EC         369      
 FIB        106      
 IMM        3        
 TAL        2        
 aPT        1 
Recessive_Genotype
 cell_type  num_hits 
 FIB        4        
 IMM        4        
 aPT        10       
'''

# Plotting #####################################################################

de_a = DE(f'{working_dir}/output/DE_Additive_Genotype').table
de_d = DE(f'{working_dir}/output/DE_Dominant_Genotype').table
de_r = DE(f'{working_dir}/output/DE_Recessive_Genotype').table

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

def plot_volcano(df, significance_column='FDR', threshold=0.05, 
                 label_top_genes=True, num_top_genes=30, min_distance=0.1,
                 colors=('blue', 'grey', 'red'), alpha=1, size=6,
                 plot_as_pdf=False, plot_directory=None):
    
    import matplotlib.patches as mpatches
    from scipy.spatial.distance import pdist, squareform
    
    df = df.to_pandas()
    df['logP'] = -np.log10(df['P'])
    df['rank'] = abs(df['logFC']) * df['logP']
    df['color'] = colors[1]
    sig = df[significance_column] < threshold
    df.loc[sig & (df['logFC'] > 0), 'color'] = colors[2]
    df.loc[sig & (df['logFC'] < 0), 'color'] = colors[0]

    plt.figure(figsize=(6, 8))
    plt.scatter(df['logFC'], df['logP'], c=df['color'], alpha=alpha, s=size)
    if label_top_genes:
        for direction in ['up', 'down']:
            top_genes = \
                df[(sig & ((df['logFC'] > 0) if direction == 'up' 
                    else (df['logFC'] < 0)))].nlargest(num_top_genes, 'rank')
            if not top_genes.empty:
                distances = squareform(pdist(top_genes[['logFC', 'logP']]))
                np.fill_diagonal(distances, np.inf)
                filter_idx = [np.min(distances[i]) > min_distance 
                              for i in range(len(distances))]
                filtered_genes = top_genes.iloc[filter_idx]
                for _, gene in filtered_genes.iterrows():
                    plt.text(gene['logFC'], gene['logP'], gene['gene'],
                             fontweight='semibold', fontsize=8)
    plt.title(f"{df['cell_type'].unique()[0]}")
    plt.xlabel(r'$\log_2(Fold Change)$')
    plt.ylabel(r'$-\log_{10}(P-value)$')
    legend_patches = [
        mpatches.Patch(color=color, label=label) 
        for color, label in zip(colors, ['down', '', 'up'])]
    plt.legend(handles=legend_patches, 
               title=f"{significance_column} < {threshold}", loc='best')
    if plot_directory is not None:
        filename = (
            f"{df['cell_type'].unique()[0].replace('/', '-')}"
            f".{'pdf' if plot_as_pdf else 'png'}"
        )
        plot_file = os.path.join(plot_directory, filename)
        savefig(plot_file)
        
for cell_type in de_a['cell_type'].unique():
   plot_volcano(de_a.filter(cell_type=cell_type), threshold=0.10,
                plot_directory=f'{working_dir}/figures/volcano')