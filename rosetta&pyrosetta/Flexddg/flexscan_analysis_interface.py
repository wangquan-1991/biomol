# setup rstoolbox:
import rstoolbox
# ====================
# hotfix rstoolbox, for metal protein;
# script_dir = os.path.dirname(os.path.realpath("__file__"))
# os.system(f'mv {script_dir}/rosetta.py /root/miniconda3/lib/python3.7/site-packages/rstoolbox/io/rosetta.py')

# load reference and design data.
import os
import rstoolbox as rs
import pandas as pd
import seaborn as sns
from pyrosetta import init
from pyrosetta.io import poses_from_silent
from multiprocessing import Pool


def extract_pdb_by_tag(tag_label, silent_file):
    """
    Summary: extract pdbs from silent by tag labels;
    Args:
        tag_label (str): tag list of the description;
        silent_file (str): raw silent file;
    """

    os.system('extract_pdbs.mpi.macosclangrelease \
              -in::file::silent %s \
              -in:file:tags %s' % (silent_file, tag_label))


def extract_pdb_from_pose(pose, top_decoy_list):
    '''
    分离pose，并检测是否存在于top rank list中.
    '''
    description = pose.pdb_info().name()
    if description in top_decoy_list:
        pose.dump_pdb(pose.pdb_info().name())
    else:
        pass


def analyze_flexddg(silent):
    # read df
    df = rs.io.parse_rosetta_file(silent)
    df_wt = df[df['description'].str.contains('wt')]
    df_mut = df[df['description'].str.contains('design')]

    # identify_mutants;
    avg_df = pd.DataFrame()
    for index, mut_info in enumerate(set(df_wt['mut_info'])):
        # selection & avg;
        df_wt_ = (df_wt[df_wt['mut_info'] == mut_info])
        df_wt_mean = df_wt_.mean()
        df_mut_ = (df_mut[df_mut['mut_info'] == mut_info])
        df_mut_mean = df_mut_.mean()

        # delta:
        delta_df = df_mut_mean - df_wt_mean
        delta_df.loc['mut_info'] = mut_info
        for terms in delta_df.index:
            avg_df.loc[index, terms] = delta_df[terms]
        avg_df.loc[index, 'description'] = ','.join(df_wt_['description'])

    # save data to csv;
    avg_df.to_csv('avg.csv')

    # Analyze decompose scores to fig;
    sns.set_style('whitegrid')

    # decomposed all energy terms;
    avg_df = avg_df.set_index('mut_info')
    for index in avg_df.index:
        fig = avg_df.loc[index, ['mut_ddG', 'fa_atr', 'fa_rep', 'fa_sol', 'fa_elec', 'fa_dun_rot',
                                 'fa_dun_semi', 'rama_prepro', 'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc',
                                 'hbond_sc', 'dG_cross', 'dSASA_int', 'dSASA_polar', 'dSASA_hphobic',
                                 'dG_separated/dSASAx100', 'hbonds_int', 'delta_unsatHbonds', 'sc_value']
                             ].plot.bar(title=f'{index}', fontsize=8, figsize=(10, 7))
        # plot;
        for tick in fig.get_xticklabels():
            tick.set_rotation(90)
        fig = fig.get_figure()
        fig.savefig('%s_decomposed.png' % index, dpi=300)
        fig.clear()

    # ddG vs Mutation Distribution;
    positions = [index.split('_')[1] for index in avg_df.index]
    set(positions)

    for position in positions:
        index_list = [index for index in avg_df.index if index.split('_')[1] == position]
        tmp_df = pd.DataFrame()
        for index in index_list:
            tmp_df.loc[index, 'ddG'] = avg_df.loc[index, 'mut_ddG']
        fig = tmp_df.plot.barh(title=f'Position {position} ddG Distribution', figsize=(10, len(index_list)+10))
        fig = fig.get_figure()
        fig.savefig('Position %s ddG Distribution' % position, dpi=300)
        fig.clear()


if __name__ == '__main__':
    silent = 'result.silent'
    analyze_flexddg(silent)
