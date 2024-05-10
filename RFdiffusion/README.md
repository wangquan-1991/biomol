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

# original structure prediction weights
wget http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt

安装：  
**RFdiffusion**的默认安装存在问题，如果使用docker会正常运行(见dockerfile)，但是使用conda会导致无法正确调用cuda，conda修改安装流程可参考：  
使用修改后的SE3nv.yml文件,手动安装cuda相应模块：  
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

