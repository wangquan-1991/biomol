import os
from pyrosetta import *
from pyrosetta.rosetta.core.simple_metrics.metrics import SelectedResiduesMetric
from pyrosetta.rosetta.core.select.residue_selector import InterGroupInterfaceByVectorSelector, AndResidueSelector, NotResidueSelector, ChainSelector
from pyrosetta.rosetta.core.simple_metrics import get_sm_data
from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
# biopython
from warnings import simplefilter
from Bio import PDB
from Bio import BiopythonWarning
from Bio.PDB.Polypeptide import *
simplefilter('ignore', BiopythonWarning)
pdb_io = PDB.PDBIO()
pdb = PDB.PDBParser()

pyrosetta.init()

pose = pose_from_pdb("1aj7_clean_0001.pdb") #输入你的初始PDB文件，最好是经过relax优化之后的文件

input_res = "45H,46L,47H,48L" #s输入残基编号以及链号，多个残基使用逗号分割

#根据输入来选取残基，这个函数非常重要,此处是在pose中选择了输出残基，之后用于设计
range_selector = ResidueIndexSelector(input_res)
selected = range_selector.apply(pose)
selected_index = SelectedResiduesMetric(range_selector)
selected_index.set_output_in_rosetta_num(True) # True: poseID, False: PDBID

selected_index = SelectedResiduesMetric(range_selector)
selected_index.set_output_in_rosetta_num(True) # True: poseID, False: PDBID
prefix = 'design_'
selected_index.apply(pose,prefix)

sm_data = get_sm_data(pose)
string_metric = sm_data.get_string_metric_data()
string_metric['design_selection']
design_pose_id_list = sm_data.get_string_metric_data()['design_selection']

#在设计残基的10A范围内选择残基,后续会用与repack设置

nbr_repack_selector = NeighborhoodResidueSelector(range_selector, 10.0, False)  # False 代表不包括核心的氨基酸。
nbr_repack_selector.apply(pose)

#选择剩余其他残基，这些残基在设计过程中不进行repack和design
from pyrosetta.rosetta.core.select.residue_selector import OrResidueSelector
core_selector = OrResidueSelector(nbr_repack_selector, range_selector)
core_selector.apply(pose)
static_selector = NotResidueSelector(core_selector)
static_selector.apply(pose)

# The task factory accepts all the task operations
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

# These are pretty standard
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())


# repack only
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), nbr_repack_selector))

#not repack,not design
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), static_selector))

# Enable design 
aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
aa_to_design.aas_to_keep("ADEFGHIKLMNPQRSTVWY") #设置你要做的突变类型
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
    aa_to_design, range_selector))
packer_task = tf.create_task_and_apply_taskoperations(pose)
print(packer_task)
print('Great!')
#set mover map
mm = pyrosetta.rosetta.core.kinematics.MoveMap()
mm.set_bb(True)
mm.set_chi(True)
mm.set_jump(True)

#init()
#init_cmd = '-corrections::beta_nov16 -corrections::beta_nov16_cart' #开启beta打分函数
#init(init_cmd)
#定义desig函数


#Set Fastdesign:

beta_nov16 = create_score_function('beta_nov16') #设置打分函数
final_design = FastDesign(beta_nov16, 3)  #选择打分函数，以及重复次数
final_design.set_task_factory(tf)
final_design.set_movemap(mm)
#final_design.apply(pose)
#final_design.cartesian(True)

circle = 5  #正式运行前，请先预计算单次循环所需时间与资源

# # design for 5 times: code for design pose.
for i in range(circle):      #具体循环次数可以进行修改，建议先进行一次循环确定计算时间，再进行进一步设置
    init()
    init_cmd = '-corrections::beta_nov16 -corrections::beta_nov16_cart' #开启beta打分函数
    init(init_cmd)
    design_pose = Pose() #创建一个名为design_pose的蛋白质结构对象，用于存储每次循环的设计结果
    design_pose.assign(pose)  # 将初始结构starting_pose复制给pose，以便在每次循环中都从相同的初始结构开始设计。
    final_design.apply(design_pose)  ## apply design
    #pose_pdbinfo = pose.pdb_info()
   # output to silent file;
    #output_pdb = "output_{i}.pdb".format(_)
    design_pose.dump_pdb( f'result_{i}.pdb') #每个设计保存为一个PDB文件，当输出过多时，可关闭此选项,f'c18_result_{i}.pdb' 是一个格式化字符串，用于生成不同循环迭代中的文件名，其中 {i} 是循环变量的占位符，将被实际循环迭代的值替换。
    poses_to_silent(design_pose, f'result_{i}.silent') #为每一个结果保存一个silent文件，便于后续能量分析，建议开启

'''
#将所有silent汇总到一个final_result.silent文件便于分析处理
final_silent = 'final_result.silent'
with open(final_silent, 'w') as f_out:
    for i in range(circle): #该数字与之前循环数字一致，无需修改
        input_silent = f'c24_result_{i}.silent'
        with open(input_silent, 'r') as f_in:
            f_out.write(f_in.read())

  # cat silents to savedpath:
    
    os.system('find . -name "*.silent" | xargs cat  > result')
    os.system('find ./ -name "*.silent" | xargs rm -rf {}')
    outname = str(random.randint(0, 99999999999999999999999)) + str(int(time.time()))
    outpath = os.path.join(savedpath, f'{outname}.silent')
    os.system(f'mv result {outpath}')


# to csv: 将关键参数放入csv文件中，便于pandas分析
#os.chdir(savedpath)
import rstoolbox as rs
rules = {'scores_ignore': [''], 'sequence': ""}
df = rs.io.parse_rosetta_file('final_result.silent', rules)
df.to_csv('final_result.csv')

# Remove individual silent files
for i in range(circle): #与循环数量相同
    input_silent = f'design_result_{i}.silent'
    os.remove(input_silent)
'''
