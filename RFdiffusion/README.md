(https://github.com/wangquan-1991/biomol/assets/103472012/b10b22c1-4c27-40c9-a56e-791bec0f759e)# RFDiffusion  

github：https://github.com/RosettaCommons/RFdiffusion  

克隆git 
`git clone https://github.com/RosettaCommons/RFdiffusion.git`  
下载权重文件：  
`cd RFdiffusion`  
`mkdir models && cd models`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt`  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt`  

可选:  
`wget http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt`  

original structure prediction weights  
wget http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt

安装：  
**RFdiffusion**的默认安装存在问题，如果使用docker会正常运行(见dockerfile)，但是使用conda会导致无法正确调用cuda，conda修改安装流程可参考：  
使用修改后的SE3nv.yml文件,手动安装cuda相应模块：  
`conda env create -f env/SE3nv.yml`  
`conda activate SE3nv`  
`cd env/SE3Transformer`  
`pip3 install --force-reinstall torch torchvision torchaudio`  
`pip install --no-cache-dir -r requirements.txt`  
`python setup.py install`  
`cd ../..` # change into the root directory of the repository  
`pip install -e . # install the rfdiffusion module from the root of the repository`  
解压测试数据：  
`tar -xvf examples/ppi_scaffolds_subset.tar.gz -C examples/`  
测试：  
`./scripts/run_inference.py 'contigmap.contigs=[150-150]' inference.output_prefix=test_outputs/test inference.num_designs=10`

## 使用方法及脚本  

`nohup /RFdiffusion/scripts/run_inference.py inference.output_prefix=/your/output/path/ inference.input_pdb=/your/input/pdb_path/H1__wis67_stem_0001.pdb 'contigmap.contigs=[A1-34/3/A35-51/3-5/A52-110/12/A112-207]' inference.num_designs=10 &` #A链1-34不设计，插入3个AA，A链35-51不设计，插入3-5个AA，A链52-110不设计，插入12个AA，A链112-207不设计，共设计10个骨架  
尝试理解一下下列指令：  
`nohup /mnt/wq/RFdiffusion/RFdiffusion-main/scripts/run_inference.py inference.output_prefix=/mnt/wq/flu/RFdiffusion/H1/ inference.input_pdb=/mnt/wq/flu/RFdiffusion/H1/H1__wis67_stem_0001.pdb 'contigmap.contigs=[A1-34/3/A35-51/3-5/A52-107/15-50/A111-207]' inference.num_designs=500 &`  

`nohup /mnt/wq/RFdiffusion/RFdiffusion-main/scripts/run_inference.py inference.output_prefix=/mnt/wq/flu/RFdiffusion/H1/3_4_15/ inference.input_pdb=/mnt/wq/flu/RFdiffusion/H1/3_4_15/H1__wis67_stem_0001.pdb 'contigmap.contigs=[A1-34/3/A35-51/4/A52-107/15/A111-207]' inference.num_designs=500 &`  

`nohup /mnt/wq/RFdiffusion/RFdiffusion-main/scripts/run_inference.py inference.output_prefix=/mnt/wq/flu/RFdiffusion/H1/3_4_25/ inference.input_pdb=/mnt/wq/flu/RFdiffusion/H1/3_4_25/H1__wis67_stem_0001.pdb 'contigmap.contigs=[A1-34/3/A35-51/4/A52-107/25/A111-207]' inference.num_designs=500 &`  


**更多案例以及参数可以参照作者的example ：https://github.com/RosettaCommons/RFdiffusion/tree/main/examples**  
**请注意，生成的骨架结构务必仔细检查（视觉），如果存在明显不合理或者不想要的的结构，人工删除，或者修改相应参数重新生成**  


# proteinMPNN  
RFdiffusion只能够生成骨架，并不会添加侧链，你可以使用rosetta或者proteinMPNN来添加侧链，作者推荐proteinMPNN（注意使用最新权重）  
链接：https://github.com/dauparas/ProteinMPNN   ***可以直接在RFdiffusion环境下调用***  

demo：
`#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --gres=gpu:rtx2080:1
#SBATCH -c 3
#SBATCH --output=example_4_non_fixed.out

#source activate mlfold

folder_with_pdbs="../inputs/PDB_complexes/pdbs/" #输入pdb

output_dir="../outputs/example_4_non_fixed_outputs" #输出路径
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="A C"
#The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
design_only_positions="1 2 3 4 5 6 7 8 9 10, 3 4 5 6 7 8" #design only these residues; use flag --specify_non_fixed  #指定设计残基，也可以指定不设计残基

python ../helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains 

python ../helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python ../helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

python ../protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 2 \  
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
`


# 使用RFdiffusion+ProteinMPNN+Fastrelax进行蛋白设计  
***实际的药物设计需要多个流程模块结合到一起***  
参考：https://github.com/nrbennet/dl_binder_design  

参考流程：  
结晶蛋白/建模 → 前处理，relax等能量最小化 → 确认设计目的以及区域，hotspot等 → 设计（长度，二级结构偏好） → score，视觉检查，筛选 → ProteinMPNN设计序列（偏好性） → AF2预测结构，PLDDT以及pAE计算（例如plddt>90, PAE<10），或者ESMfold等 → 视觉检查，筛选  →  ...  

参考链接：https://github.com/nrbennet/dl_binder_design  

## Install  
***使用AF2预测结构效果会更好，但是计算成本会更高，所以使用阉割版的AF2initial，或者ESMfold，但是ESMfold需要先用WT测试一下***  
1. 第三方库，例如AF2及ESMfold等
2. pyrosetta, 按照相应流程安装
3. 安装模块
`git clone https://github.com/nrbennet/dl_binder_design.git`  
`cd include/`  
`conda env create -f proteinmpnn_fastrelax.yml`





























