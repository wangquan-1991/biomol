import os
from pyrosetta import *
from pyrosetta.rosetta.core.simple_metrics.metrics import SelectedResiduesMetric
from pyrosetta.rosetta.core.select.residue_selector import InterGroupInterfaceByVectorSelector, AndResidueSelector, NotResidueSelector, ChainSelector
from pyrosetta.rosetta.core.simple_metrics import get_sm_data

from warnings import simplefilter
from Bio import PDB
from Bio import BiopythonWarning
from Bio.PDB.Polypeptide import *
simplefilter('ignore', BiopythonWarning)
pdb_io = PDB.PDBIO()
pdb = PDB.PDBParser()

std_aa = ['ASP', 'ALA', 'ARG', 'ASN', 'GLU', 'HIS', 'ILE', 'CYS', 'LEU',
         'LYS', 'MET', 'PRO', 'PHE', 'GLN', 'SER', 'THR', 'TRP', 'TYR',
         'VAL']

std_daa = ['DPH', 'DTR', 'DVA', 'DCS', 'DAS', 'DGU', 'DHI', 'DIL', 'DLY',
             'DLE', 'DME', 'DAN', 'DGN', 'DAR', 'DSE', 'DTH', 'DTY', 'DAL', 'DPR']

std_aa_3to1 = {
    'ALA': 'A', 
    'ARG': 'R', 
    'ASN': 'N', 
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y', 
    'VAL': 'V'  
}
# std_aa2daa = {'ALA': "DAL", 'CYS': "DCS", 'ASP': "DAS", 'GLU': "DGU", "PHE": "DPH", "GLY": "G", "HIS": "DHI", "ILE": "DIL",
#               "LYS": "DLY", "LEU": "DLE", "MET": "DME", "ASN": "DAN", "PRO": "DPR", "GLN": "DGN", "ARG": "DAR",
#               "SER": "DSE", "THR": "DTH", "VAL": "DVA", "TRP": "DTR", "TYR": "DTY"}

std_aa2daa = {'ALA': "DALA", 'CYS': "DCYS", 'ASP': "DASP", 'GLU': "DGLU", "PHE": "DPHE", "GLY": "GLY", "HIS": "DHIS", "ILE": "DILE",
              "LYS": "DLYS", "LEU": "DLEU", "MET": "DMET", "ASN": "DASN", "PRO": "DPRO", "GLN": "DGLN", "ARG": "DARG",
              "SER": "DSER", "THR": "DTHR", "VAL": "DVAL", "TRP": "DTRP", "TYR": "DTYR"}

xrotlib_aa = ['L00', 'L01', 'L02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P12', 'P18', 'P24', 'P25', 'P26', 'P27', 'P31', 'P32', 'P36', 'P37', 'P43', 'P44', 'P45', 'C20', 'P46', 'P47', 'P48', 'P49', 'P50', 'P55', 'P56', 'P57', 'P58', 'P59', 'P60', 'P61', 'P62', 'B12', 'P63', 'P64', 'P65', 'P66', 'P67', 'P68', 'P69', 'P70', 'P71', 'P72', 'P73', 'P74', 'P79', 'P80', 'P81', 'P82', 'P83', 'P85', 'P86', 'P87', 'P88', 'P89', 'P90', 'P91', 'P93', 'P94', 'P95', 'P96', 'P97', 'P98', 'P99', 'R00', 'R01', 'R02', 'R03', 'R04', 'R05', 'R06', 'R07', 'R09', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18', 'R19', 'R20', 'R21', 'R23', 'R24', 'R25', 'R26', 'V01', 'R27', 'R28', 'R29', 'R30', 'R31', 'R32', 'R33', 'R34', 'R35', 'R36', 'R37', 'R38', 'R39', 'R40', 'R41', 'R42', 'R43', 'R44', 'R45', 'R46', 'R47', 'R48', 'R49', 'R50', 'R51', 'R52', 'R53', 'R54', 'R55', 'R56', 'R57', 'R58', 'R59', 'R60', 'R62', 'R63', 'R64', 'R65', 'R66', 'R67', 'R68', 'R69', 'R70', 'R71', 'R72', 'R73', 'R74', 'R75', 'R76', 'R77', 'R78', 'R79', 'R80', 'R81', 'R82', 'R83', 'R84', 'R85', 'R86', 'R87', 'R88', 'R89', 'C87', 'R91', 'R92', 'R93', 'R94', 'R95', 'R97', 'R99', 'X00', 'X01', 'X02', 'X04', 'X05', 'X06', 'X07', 'X08', 'X09', 'X13', 'X14', 'X15', 'X16', 'X17', 'X18', 'X19', 'X20', 'X21', 'B19', 'X23', 'X24', 'X26', 'X27', 'X30', 'X31', 'X32', 'X33', 'X34', 'X35', 'X36', 'X37', 'X38', 'X39', 'X40', 'X41', 'A30', 'X42', 'X43', 'X44', 'X45', 'X46', 'X47', 'X50', 'X52', 'X53', 'X54', 'B50', 'C05', 'X55', 'X57', 'X58', 'NLU', 'A92', 'X60', 'X61', 'X66', 'X71', 'X72', 'X73', 'X74', 'X75', 'X76', 'X77', 'X78', 'X79', 'X80', 'X81', 'X82', 'X83', 'X84', 'B63', 'X85', 'X87', 'X88', '0TD', 'X91', 'X93', 'X95', 'A83', 'A07', 'X96', 'X97', 'X98', 'X99', 'P00', 'P01', 'P02', 'P51', 'P52', 'A04', 'A98', 'B02', 'B06', 'C15', 'C16', 'C89', 'A24']

def make_saturation_resfile(pdbfile, partner_chain, cutoff, use_ncaa=False, savedpath=os.getcwd()):

    # load complex pose;
    init()
    pose = pose_from_pdb(pdbfile)

    # for cyclic-peptide:
    #if cyclization_type:
     #   connect_termini(pose, cyclization_type)

    # check pdbnum in partner chain;
    # partner_chain = ','.join(partner_chain)
    group1 = ChainSelector(partner_chain)
    group2 = NotResidueSelector(group1)
    interface_seletor = InterGroupInterfaceByVectorSelector(group1, group2)
    interface_seletor.vector_dist_cut(float(cutoff))
    interface_seletor.nearby_atom_cut(8.0)
    designable_selector = AndResidueSelector(interface_seletor, group1)

    #
    selected_index = SelectedResiduesMetric(designable_selector)
    selected_index.set_output_in_rosetta_num(True)  # output as pose id.
    selected_index.apply(pose)
    sm_data = get_sm_data(pose)
    selected_pose_id_list = sm_data.get_string_metric_data()['selection']
    print(f"ppi residue: {selected_pose_id_list}")

    # get peptide start/end res:
    start_res = pose.chain_begin(pose.num_chains())
    end_res = pose.chain_end(pose.num_chains())

    for pose_id in selected_pose_id_list.split(','):
       # if 1==2 :
        #    # 不允许突变:
         #   continue
        info = pose.pdb_info().pose2pdb(int(pose_id))
        pdb_num = info.split(' ')[0]
        chain_id = info.split(' ')[1]
        phi_angles = pose.phi(int(pose_id))
        native_residue = pose.residue(int(pose_id)).name3()
        write_resfile(pdb_num, chain_id, native_residue, phi_angles, use_ncaa, savedpath)


def write_resfile(pdb_num, chain_id, wt, phi_angles, use_ncaa, cyclization_type, savedpath=os.getcwd()):
    '''
    wt: name3; / mut: name3;
    '''
    for mut in std_aa:
        mut_1letter = std_aa_3to1[mut]
        with open(os.path.join(savedpath, '%s%s%s_%s.resfile' % (wt, pdb_num, mut_1letter, chain_id)), 'w') as f:
            f.write('NATAA\nstart\n') 
            f.write('%s %s PIKAA %s\n' % (pdb_num, chain_id, mut_1letter))
        #with open(os.path.join(savedpath, '%s%s%s_%s.resfile' % (wt, pdb_num, mut, chain_id)), 'w') as f:
         #   f.write('NATAA\nstart\n')  # header
          #  f.write('%s %s PIKAA %s\n' % (pdb_num, chain_id, mut))

    # FOR NCAA, 默认使用D型氨基酸，只能处理有环肽的结构;
    #if use_ncaa:
     #   for mut in xrotlib_aa:
      #      if phi_angles > 0: # D residue type mut here;
       #         mut = 'D'+mut
        #    with open(os.path.join(savedpath, '%s%s%s_%s.resfile' % (wt, pdb_num, mut, chain_id)), 'w') as f:
         #       f.write('NATAA\nstart\n')  # header
          #      f.write('%s %s PIKAA X[%s]\n' % (pdb_num, chain_id, mut))

if __name__ == '__main__':
    pdbfile = '/data/wq/antibody/integrin/3v4p_50ns_0001_0001.pdb'
    chain_id = 'C,D'
    neighbor_cutoff = '10'
    make_saturation_resfile(pdbfile, chain_id, neighbor_cutoff)