# Copyright (c) 2020 by weikun wu. All Rights Reserved.
#
from pyrosetta.io import poses_to_silent
import pyrosetta
import os
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, OperateOnResidueSubset, PreventRepackingRLT, ReadResfile, RestrictToRepackingRLT, UseMultiCoolAnnealer
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector,
    NeighborhoodResidueSelector, NotResidueSelector,
    PrimarySequenceNeighborhoodSelector, ResidueIndexSelector)
from pyrosetta.rosetta.protocols import rosetta_scripts
from pyrosetta.rosetta.protocols.constraint_movers import ClearConstraintsMover
from pyrosetta.rosetta.protocols.minimization_packing import (
    MinMover, PackRotamersMover)
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.pack.palette import CustomBaseTypePackerPalette
from pyrosetta.rosetta.core.pose import get_chain_from_chain_id, get_chain_id_from_chain
from pyrosetta.rosetta.protocols.simple_moves import DeleteChainMover
from pyrosetta.rosetta.core.pose import renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.utility import vector1_int
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.core.pose import setPoseExtraScore
import random
from pyrosetta.rosetta.protocols.moves import MoverStatus
from pyrosetta.rosetta.protocols.generalized_kinematic_closure.filter import filter_type
from pyrosetta.rosetta.core.pose import remove_lower_terminus_type_from_pose_residue, remove_upper_terminus_type_from_pose_residue
from pyrosetta.rosetta.protocols.cyclic_peptide import *
from pyrosetta.rosetta.protocols.generalized_kinematic_closure import GeneralizedKIC
from pyrosetta.rosetta.protocols.cyclic_peptide import TryDisulfPermutations
from pyrosetta.rosetta.protocols.rosetta_scripts import ParsedProtocol
from pyrosetta.rosetta.protocols.generalized_kinematic_closure.selector import selector_type

def delete_chain_from_pose(pose, pdb_chain):
    # get to_delete pose chain id;
    to_delete_chain = get_chain_id_from_chain(pdb_chain, pose)
    print(f'****DeleteChainsFromPose****: {pdb_chain}(chain_id:{to_delete_chain})')

    # assign pose;
    pose_ = Pose().assign(pose.clone())

    # delete chains;
    delete_chain = DeleteChainMover()
    delete_chain.chain_num(to_delete_chain)
    delete_chain.apply(pose_)
    return pose_


def extract_chain_from_pose(pose, pdb_chain):
    # get to_delete pose chain id;
    to_extract_chain = get_chain_id_from_chain(pdb_chain, pose)

    # split_by_chain;
    pose_ = pose.split_by_chain()[to_extract_chain]
    print(f'****SplitChainsFromPose****: {pdb_chain}(chain_id:{to_extract_chain})')

    return pose_


def get_partner(pose1, pose2):
    '''
    setup docking_foldtree by pose1_pose2;
    '''
    # get partner:
    pose1_chains = []
    for i in range(1, pose1.num_chains()+1):
        pose1_chains.append(get_chain_from_chain_id(i, pose1))
    pose1_chains = ''.join(pose1_chains)

    pose2_chains = []
    for j in range(1, pose2.num_chains()+1):
        pose2_chains.append(get_chain_from_chain_id(j, pose2))
    pose2_chains = ''.join(pose2_chains)

    return pose1_chains + '_' + pose2_chains


def megre_to_complex_pose(pose1, pose2):
    '''
    add pose2 to pose1 by the last residue of pose1;
    setup the docking foldtree;
    '''
    # add pose;
    complex_pose = Pose()
    complex_pose.assign(pose1)
    complex_pose.append_pose_by_jump(pose2, pose1.total_residue())

    # update infos;
    complex_pose.update_residue_neighbors()
    complex_pose.update_pose_chains_from_pdb_chains()
    complex_pose.conformation().detect_disulfides()  # 可能会破坏二硫键的定义，导致genKIC再找不到SG connect。
    # renumber_pdbinfo_based_on_conf_chains(complex_pose)  # renumber, 否则编号会导致后续程序出错!
    complex_pose.pdb_info().obsolete(False)

    # no jump moveable;
    partner = get_partner(pose1, pose2)
    v = vector1_int(1)
    setup_foldtree(complex_pose, partner, v)

    return complex_pose, partner

def get_mutation(pose, resfile):
    icode_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                  'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    mut_info = []
    with open(resfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'PIKAA' in line:
                info = line.strip('\n').split(' ')
                # for antibody:
                if info[0][-1] in icode_list:
                    icode_ = info[0][-1]
                else:
                    icode_ = None
                mut_pdb_id = info[0]+info[1]
                mut_res = info[3]
                if icode_:
                    wt_res = pose.residue(pose.pdb_info().pdb2pose(info[1], int(info[0][:-1]), icode_)).name1()
                else:
                    wt_res = pose.residue(pose.pdb_info().pdb2pose(info[1], int(info[0]))).name1()
                mut_info.append(wt_res+"_"+mut_pdb_id+"_"+mut_res)
    mut_info
    return mut_info

class FlexDdgDesignMover(object):
    def __init__(self, pdbfile, partner_chain, interface_mode=False, pssm_resfile=None,
                 position_list=None, resfile=None, ntrials=35000, scorefnx='beta_nov16',
                 extra_res_fa_args=None, extra_res_list=None, cyclization_type=None, debug=False):

        if extra_res_fa_args:
            init(f'-ex1 -ex2 -input_ab_scheme Chothia_Scheme -corrections::beta_nov16 -extrachi_cutoff 0 \
                   -in:auto_setup_metals \
                   -extra_res_fa crosslinker/1.4_bisbromomethylbenzene.params crosslinker/1.3.5_trisbromomethylbenzene.params sidechain_conjugation/CYX.params {extra_res_fa_args}')
        else:
            init('-ex1 -ex2 -input_ab_scheme Chothia_Scheme -corrections::beta_nov16 -extrachi_cutoff 0 \
                  -in:auto_setup_metals \
                  -extra_res_fa crosslinker/1.4_bisbromomethylbenzene.params crosslinker/1.3.5_trisbromomethylbenzene.params sidechain_conjugation/CYX.params')

        self.chain_to_move = partner_chain
        self.pose = pose_from_pdb(pdbfile)
        self.debug = debug
        # support for cyclic_peptide:
        self.cyclization_type = cyclization_type
        if self.cyclization_type:
            connect_termini(self.pose, self.cyclization_type)

        self.resfile = resfile
        self.pssm_resfile = pssm_resfile
        self.extra_res_fa_args = extra_res_fa_args
        self.extra_res_list = extra_res_list

        # mutation:
        self.mut_info = ','.join(get_mutation(self.pose, resfile))

        # setup scorefuntion.
        self.scorefnx = create_score_function(scorefnx)
        self.scorefnx_cst = create_score_function(scorefnx+'_cst')

        # setup foldtree:
        from pyrosetta.rosetta.core.pose import get_chain_from_chain_id
        all_chains = [get_chain_from_chain_id(chain_id, self.pose) for chain_id in range(1, self.pose.num_chains()+1)]
        print(all_chains)

        # partner1(receptor)
        partner1 = self.chain_to_move.split(',')
        print(partner1)
        # partner2(cal_chains)
        partner2 = [chain_id for chain_id in all_chains if chain_id not in partner1]
        print(partner2)

        # receptor
        receptor_pose = Pose().assign(self.pose.clone())
        for chain_id in partner1:
            receptor_pose = delete_chain_from_pose(receptor_pose, chain_id)
            receptor_pose.pdb_info().obsolete(False)

        # cal_chain pose;
        cal_pose = Pose().assign(self.pose.clone())
        for chain_id in partner2:
            cal_pose = delete_chain_from_pose(cal_pose, chain_id)
            cal_pose.pdb_info().obsolete(False)

        # set docking foldtree;
        self.pose, self.partners = megre_to_complex_pose(receptor_pose, cal_pose)
        print(self.pose.pdb_info())
        print(self.pose.fold_tree())
        print(self.partners)
        if self.cyclization_type:
            connect_termini(self.pose, self.cyclization_type)

        if self.cyclization_type == 'n_to_c_amide_bond':
            start_res = self.pose.chain_begin(self.pose.num_chains())
            end_res = self.pose.chain_end(self.pose.num_chains())
            add_design_NC_constraint(self.pose, start_res, end_res)

        if interface_mode:
            # enable all interface residue in patner_chain to be designed.
            group1 = ChainSelector(self.chain_to_move)
            group2 = NotResidueSelector(group1)
            interface_seletor = InterGroupInterfaceByVectorSelector(group1, group2)
            interface_seletor.nearby_atom_cut(6)
            interface_seletor.vector_dist_cut(8)
            self.designable_selector = AndResidueSelector(interface_seletor, group1)

        elif position_list:
            # manually mode;
            pose_id_list = ','.join([str(self.pose.pdb_info().pdb2pose(
                pdbnumber[-1], int(pdbnumber[:-1]))) for pdbnumber in position_list.split(',')])
            self.designable_selector = ResidueIndexSelector(pose_id_list)
            if pssm_resfile:
                raise Exception('PSSM profile can only use in Interface Mode')

        elif resfile:
            # user define designable residue manually
            self.designable_selector = self.parse_resfile(resfile)
            if pssm_resfile:
                raise Exception('PSSM profile can only use in Interface Mode')
        else:
            raise Exception('No Running Mode Defined!')

        # set backrup trail number.
        self.ntrials = ntrials

        # add D residue-type
        ncaa_list = 'DALA,DCYS,DASP,DGLU,DPHE,DHIS,DILE,DLYS,DLEU,DMET,DASN,DPRO,DGLN,DARG,DSER,DTHR,DVAL,DTRP,DTYR'
        if self.extra_res_list:  # load user ncaa
            ncaa_list = ncaa_list + ',' + self.extra_res_list
        print(ncaa_list)
        self.pp = CustomBaseTypePackerPalette()
        self.pp.parse_additional_residue_types(ncaa_list)
    def apply(self, prefix, job_id):

        # 产生新的随机数.
        # init('-ex1 -ex2 -corrections::beta_nov16 -extrachi_cutoff 0 -in:auto_setup_metals')
        if self.extra_res_fa_args:
            init(f'-ex1 -ex2 -corrections::beta_nov16 -input_ab_scheme Chothia_Scheme -extrachi_cutoff 0 -in:auto_setup_metals -extra_res_fa {self.extra_res_fa_args}')
        else:
            init('-ex1 -ex2 -corrections::beta_nov16 -input_ab_scheme Chothia_Scheme -extrachi_cutoff 0 -in:auto_setup_metals')

        # setup taskop & packer & movers;
        self.setup()
        print('setup done')

        # Sampling Flexibility
        pose_ = Pose().assign(self.pose)
        self.backrub_protocol.apply(pose_)

        # Repack to generate wildtype pose
        wt_bound_pose = Pose()
        wt_bound_pose.assign(pose_)
        # support for cyclic_peptide:
        if self.cyclization_type:
            connect_termini(self.pose, self.cyclization_type)
        self.repacker.apply(wt_bound_pose)
        if self.debug:
            wt_bound_pose.dump_pdb('debug_repack.pdb')
        self.minimization(wt_bound_pose)
        if self.debug:
            wt_bound_pose.dump_pdb('debug_min.pdb')

        # Design mutation for better affinity
        mut_bound_pose = Pose()
        mut_bound_pose.assign(pose_)
        # support for cyclic_peptide:
        if self.cyclization_type:
            connect_termini(self.pose, self.cyclization_type)
        self.designer.apply(mut_bound_pose)
        self.minimization(mut_bound_pose)

        # dump to silent;
        ddg = self.interface_analysis(self.partners, mut_bound_pose)
        setPoseExtraScore(mut_bound_pose, "mut_ddG", ddg)
        setPoseExtraScore(mut_bound_pose, "mut_info", self.mut_info)
        mut_bound_pose.pdb_info().name('designed_%s_%s.pdb' % (prefix, job_id))  # rename decoy.
        poses_to_silent(mut_bound_pose, 'designed_%s_%s.silent' % (prefix, job_id))

        # dump wt data to silent;
        ddg = self.interface_analysis(self.partners, wt_bound_pose)
        setPoseExtraScore(wt_bound_pose, "mut_ddG", ddg)
        setPoseExtraScore(wt_bound_pose, "mut_info", self.mut_info)
        wt_bound_pose.pdb_info().name('wt_%s_%s.pdb' % (prefix, job_id))  # rename decoy.
        poses_to_silent(wt_bound_pose, 'wt_%s_%s.silent' % (prefix, job_id))

        if self.debug:
            mut_bound_pose.dump_pdb(mut_bound_pose.pdb_info().name())
            wt_bound_pose.dump_pdb(wt_bound_pose.pdb_info().name())
            
    def parse_resfile(self, resfile_name):
        xml = rosetta_scripts.XmlObjects.create_from_string("""
        <TASKOPERATIONS>
            <ReadResfile name="read_resfile" filename="%s"/>
        </TASKOPERATIONS>
        <RESIDUE_SELECTORS>
            <Task name="mutable_selector" fixed="0" packable="0" designable="1" task_operations="read_resfile"/>
        </RESIDUE_SELECTORS>
        """ % resfile_name)
        designable_selector = xml.get_residue_selector('mutable_selector')
        return designable_selector

    def minimization(self, pose):
        self.addcst.apply(pose)
        self.min.apply(pose)
        self.clearcst.apply(pose)

    def setup(self):
        # define nbr_shell with +-1 adjacent residue in sequence.
        nbr_shell = NeighborhoodResidueSelector(self.designable_selector, 10.0, True)
        # nbr_only = NeighborhoodResidueSelector(self.designable_selector, 10.0, False)
        not_designable = NotResidueSelector(self.designable_selector)
        self.nbr_adjacent_shell = PrimarySequenceNeighborhoodSelector(1, 1, nbr_shell)
        self.nbr_only_adjacent_shell = AndResidueSelector(self.nbr_adjacent_shell, not_designable)
        context_selector = NotResidueSelector(self.nbr_adjacent_shell)

        # repack only:
        repackonly = OperateOnResidueSubset(RestrictToRepackingRLT(), self.nbr_adjacent_shell)
        repack_nb_only = OperateOnResidueSubset(RestrictToRepackingRLT(), self.nbr_only_adjacent_shell)

        # prevent packing:
        norepack = OperateOnResidueSubset(PreventRepackingRLT(), context_selector)

        # as
        multicool = UseMultiCoolAnnealer(6)

        # set task op
        self.repack_tf = TaskFactory()
        self.repack_tf.push_back(InitializeFromCommandline())
        self.repack_tf.push_back(repackonly)
        self.repack_tf.push_back(norepack)
        self.repack_tf.push_back(multicool)
        self.repack_tf.set_packer_palette(self.pp)

        # set design taskop
        self.mutate_tf = TaskFactory()
        self.mutate_tf.push_back(InitializeFromCommandline())
        self.mutate_tf.push_back(repack_nb_only)
        self.mutate_tf.push_back(norepack)
        self.mutate_tf.push_back(multicool)
        self.mutate_tf.set_packer_palette(self.pp)

        # Resfile Mode;
        if self.resfile:
            self.mutate_tf.push_back(ReadResfile(self.resfile))

        # Using PSSM Profile;
        if self.pssm_resfile:
            self.mutate_tf.push_back(ReadResfile(self.pssm_resfile))

        task = self.mutate_tf.create_task_and_apply_taskoperations(self.pose)
        # task = self.mutate_tf.create_task_and_apply_taskoperations(self.pose)
        print(task)

        # packer
        self.repacker = PackRotamersMover()
        self.repacker.task_factory(self.repack_tf)
        self.repacker.score_function(self.scorefnx)

        # Designer
        self.designer = PackRotamersMover()
        self.designer.task_factory(self.mutate_tf)
        self.designer.score_function(self.scorefnx)

        # clear cst
        self.clearcst = ClearConstraintsMover()

        # MinMover
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)

        self.min = MinMover()
        self.min.set_movemap(movemap)
        self.min.score_function(self.scorefnx_cst)
        self.min.min_type('lbfgs_armijo_nonmonotone')
        # self.min.tolerance(0.000001)
        self.min.max_iter(200)
        # self.min.abs_score_convergence_threshold(1.0)  # 会导致ncaa不收敛一直算min的问题。

        # set up Backrup Protocol & constraint :
        pyrosetta.rosetta.utility.vector1_unsigned_long()
        xml = rosetta_scripts.XmlObjects.create_from_string("""
        <MOVERS>
            <BackrubProtocol name="backrub" mc_kt="1.2" ntrials="%s" recover_low="0"/>
            <AddConstraintsToCurrentConformationMover name="addcst" use_distance_cst="1"
                coord_dev="0.5" min_seq_sep="0" max_distance="10" CA_only="1"
                bound_width="0.0" cst_weight="0.0"/>
        </MOVERS>
        """ % self.ntrials)
        self.backrub_protocol = xml.get_mover('backrub')
        self.backrub_protocol.set_taskfactory(self.repack_tf)
        self.backrub_protocol.set_scorefunction(self.scorefnx)

        if self.cyclization_type:
            # # for skip cyclic-peptide backrup，these will error in backrup;
            peptide_selector = ChainSelector(self.pose.num_chains())

            # # trick to let cyclic-peptide move
            self.nbr_adjacent_shell = AndResidueSelector(self.nbr_adjacent_shell, NotResidueSelector(peptide_selector))

        # backrup on the receptor envs;
        self.backrub_protocol.set_pivots_from_residue_subset(self.nbr_adjacent_shell.apply(self.pose))

        self.addcst = xml.get_mover('addcst')

        # # perturb the peptide structure on the new pose:
        # if self.cyclization_type:
        #     # perturb peptide using genkic:
        #     genkic_mover = self.setup_genkic_mover(self.cyclization_type, self.pose)
        #
        #     # try to perturb 10 times:
        #     print('RUNNING: KIC on ')
        #     for i in range(10):
        #         genkic_mover.apply(self.pose)
        #         if genkic_mover.get_last_move_status() == MoverStatus.MS_SUCCESS:
        #             break
        #         else:
        #             continue

    def setup_genkic_mover(self, cyclization_type, pose):
        '''
        用途: Perturb the cyclic-peptide in the last chain;
        '''
        # NC termini;
        start_res = pose.chain_begin(pose.num_chains())
        end_res = pose.chain_end(pose.num_chains())

        # used for nc closure
        if cyclization_type == 'n_to_c_amide_bond':
            remove_lower_terminus_type_from_pose_residue(pose, start_res)
            remove_upper_terminus_type_from_pose_residue(pose, end_res)
            # find res atom name for N-terminus, C-terminus;
            firstatom = pose.residue(start_res).atom_name(pose.residue(start_res).lower_connect_atom())
            lastatom = pose.residue(end_res).atom_name(pose.residue(end_res).upper_connect_atom())

        # declare termini:
        declare_bonds = DeclareBond()
        if cyclization_type == 'terminal_disulfide':
            # res1: int, atom1: str, res2: int, atom2: str, add_termini: bool, run KIC: bool?
            declare_bonds.set(end_res, 'SG', start_res, 'SG', True, False)
        else:
            declare_bonds.set(end_res, lastatom, start_res, firstatom, False, False)
        declare_bonds.apply(pose)

        # define genKIC
        genkic_mover = GeneralizedKIC()
        genkic_mover.set_selector_scorefunction(self.scorefnx)
        genkic_mover.set_closure_attempts(100)  # small cycles, good for sampling more initialize motif rama;
        genkic_mover.add_filter(filter_type.loop_bump_check)
        genkic_mover.set_min_solution_count(1)

        # # preselect mover;
        # TryDisulfPermutations for disulfide-brigde peptide;
        try_ss = TryDisulfPermutations()
        try_ss.set_consider_already_bonded(False)
        try_ss.set_selector(NotResidueSelector(ResidueIndexSelector(f'{start_res},{end_res}'))) # not NC bond;

        # set preselector:
        filter_steps = ParsedProtocol()
        filter_steps.add_step(try_ss, 'TryDisulfPermutations', None)
        genkic_mover.set_preselection_mover(filter_steps)

        # random select one designable residue within 5 residue distance;
        designable_resid = [index+1 for index, i in enumerate(self.designable_selector.apply(self.pose)) if i == 1]
        anchor_res = random.choice(designable_resid)

        # 智能判断cutpoint:
        cutpoint = pose.chain_end(pose.num_chains())

        # move if the anchor res in the terminus:
        if anchor_res == start_res:
            anchor_res += 1
        elif anchor_res == end_res:
            anchor_res -= 2

        # Set Pivots:
        if cyclization_type == 'terminal_disulfide':
            # Set Pivots: (rsd1: int, at1: str, rsd2: int, at2: str, rsd3: int, at3: str)
            genkic_mover.set_pivot_atoms(anchor_res-1, 'CA', cutpoint, 'SG', anchor_res+1, 'CA')
        elif cyclization_type == 'n_to_c_amide_bond':
            # NC-closure
            # Set Pivots: (rsd1: int, at1: str, rsd2: int, at2: str, rsd3: int, at3: str)
            genkic_mover.set_pivot_atoms(anchor_res-1, 'CA', cutpoint, 'CA', anchor_res+1, 'CA')
        else:
            raise Exception('un-support cyclization_type!')

        # add loop residue;
        backwardN = list(reversed(list(range(start_res, anchor_res))))
        backwardC = list(reversed(list(range(anchor_res+1, end_res+1))))
        loop_indexs = backwardN + backwardC

        # report Loop pdb number;
        print([pose.pdb_info().pose2pdb(num) for num in loop_indexs])
        print(loop_indexs)
        for i in loop_indexs:
            genkic_mover.add_loop_residue(i)

        # add small_perturber & loop residues;
        from pyrosetta.rosetta.utility import vector1_core_id_NamedAtomID
        from pyrosetta.rosetta.core.id import NamedAtomID
        genkic_mover.add_perturber('perturb_dihedral')  #
        for i in loop_indexs:
            if cyclization_type == 'terminal_disulfide' and i in [start_res, end_res]:
                continue
            genkic_mover.add_residue_to_perturber_residue_list(i)
            named_atomID1 = vector1_core_id_NamedAtomID()
            named_atomID1.append(NamedAtomID("N", i))
            named_atomID1.append(NamedAtomID("CA", i))
            named_atomID2 = vector1_core_id_NamedAtomID()
            named_atomID2.append(NamedAtomID("CA", i))
            named_atomID2.append(NamedAtomID("C", i))
            genkic_mover.add_atomset_to_perturber_atomset_list(named_atomID1)
            genkic_mover.add_atomset_to_perturber_atomset_list(named_atomID2)
            genkic_mover.add_value_to_perturber_value_list(0.75)  # 扰动0.01°;
        genkic_mover.set_selector_type(selector_type.lowest_rmsd_selector)  # rmsd;

        # close bond;
        # rsd1: int, at1: str, rsd2: int, at2: str, rsd1_before: int, at1_before: str, rsd2_after: int, at2_after: str, bondlength: float, bondangle1: float, bondangle2: float, torsion: float, randomize_this_torsion: bool, randomize_flanking_torsions: bool
        if cyclization_type == 'terminal_disulfide':
            genkic_mover.close_bond(end_res, 'SG', start_res, 'SG', end_res, 'CB', start_res, 'CB', 2.05, 103, 103, 90, True, True)
            genkic_mover.set_correct_polymer_dependent_atoms(True)
        elif cyclization_type == 'n_to_c_amide_bond':
            # http://cptweb.cpt.wayne.edu/DbD2/help.php#Technical
            genkic_mover.close_bond(end_res, lastatom, start_res, firstatom, end_res, "CA", start_res, "CA", 1.328685, 116.199993, 121.699997, 180, False, False)
            genkic_mover.set_correct_polymer_dependent_atoms(True)
        else:
            raise Exception('un-support cyclization_type!')

        return genkic_mover

    def interface_analysis(self, partners, pose):
        # interface analyzing
        self.scorefnx(pose)

        # analysis
        iam = InterfaceAnalyzerMover(partners, True, self.scorefnx, True, True, True, True)  # no packstat
        iam.apply(pose)

        # ddg Filter jump
        from pyrosetta.rosetta.protocols.simple_ddg import DdgFilter
        jump = 1
        ddg_filter = DdgFilter(0, self.scorefnx, jump, 5)  # 此处应该设置为需要的cutoff，与原有的野生型对照
        ddg_filter.repack(True)
        ddG = ddg_filter.score(pose)

        return ddG


resfile_name = [file for file in os.listdir() if file.endswith('.resfile')][0]
design = FlexDdgDesignMover('1GF1R_R10.pdb', 'C', resfile=resfile_name, ntrials=35000)
ddg = design.apply(prefix='my_output', job_id=1)
    