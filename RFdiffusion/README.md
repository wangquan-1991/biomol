# RFDiffusion  

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













