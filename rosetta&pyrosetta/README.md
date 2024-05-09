rosetta+pyrosetta
用处：
蛋白处理，数据处理，结构优化，线肽、环肽设计，蛋白设计，抗体亲和力成熟，纳米颗粒设计等，还会作为某些AI模型(尤其是baker组开发的)的底层支持

安装：
下载：https://www.pyrosetta.org/downloads   使用学术邮箱注册下载(根据你的系统以及python脚本选择正确的下载文件) #商业化使用请购买，不贵
相应安装包可以见10.10.113.62服务器  /HDD/disk1/wq/soft
下载对应版本的whl文件，安装不易报错
#pyrosetta 安装
sudo apt install python3 python3-pip
pip3 install pyrosetta-2024.19+release.a34b73c-cp312-cp312-linux_x86_64.whl


安装rosetta：下载最新版rosetta，否则你可能需要python2 并且指定python is python2
#安装支持库
sudo apt install build-essential python2 python-is-python2 openmpi-bin  openmpi-common openmpi-doc libopenmpi3 libopenmpi-dev  sudo -H pip3 install SCons
tar -xzvf rosetta_*.tgz  
cd rosetta_*  
cd main/source 
./scons.py -j mode=release bin  
./scons.py bin mode=release extras=mpi  
./scons.py -j mode=release cat=test 
运行测试：
python test/run.py

RFdiffusion/RFdiffusionAA/RFdiffusion antibody
https://github.com/baker-laboratory/rf_diffusion_all_atom

https://github.com/baker-laboratory/RoseTTAFold-All-Atom



