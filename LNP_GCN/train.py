import os
os.environ['OMP_NUM_THREADS'] = '1'  # 设置OpenMP线程数，确保多线程不会干扰计算
import numpy as np
import torch
import torch.nn as nn
from torch.nn import functional as F
import random
import dgl
from dgl.nn.pytorch import GraphConv  # 图卷积层
from dgl.data import DGLDataset  # DGL数据集基础类
from dgl.dataloading import GraphDataLoader  # 图数据加载器
from sklearn.metrics import r2_score  # R2评分
from scipy.stats import spearmanr  # Spearman相关系数
from dgllife.utils import CanonicalAtomFeaturizer  # 分子图的原子特征提取器
from rdkit import Chem  # RDKit化学分子处理库
from dgllife.utils import mol_to_bigraph  # RDKit分子转图的工具函数
import pandas as pd  # 数据处理库
import optuna  # 添加Optuna库
import json  # 添加JSON库

def setup_seed(seed):
    dgl.random.seed(seed)  # 设置DGL的随机种子
    torch.manual_seed(seed)  # 设置PyTorch的随机种子
    torch.cuda.manual_seed_all(seed)  # 设置CUDA的随机种子
    np.random.seed(seed)  # 设置NumPy的随机种子
    random.seed(seed)  # 设置Python随机模块的种子
    torch.backends.cudnn.deterministic = True  # 确保每次计算的结果一致
    torch.backends.cudnn.benchmark = True  # 允许cudnn优化加速
setup_seed(42)  # 设置随机种子为42

def process(one_smiles, num_virtual_nodes=0):
    m = Chem.MolFromSmiles(one_smiles)  # 从SMILES字符串创建分子对象
    node_enc = CanonicalAtomFeaturizer()  # 原子特征提取器
    g = mol_to_bigraph(m, True, node_enc, None, False)  # 转换分子为无向图
    if not num_virtual_nodes:  # 如果不需要虚拟节点，返回该图
        return g
    else:
        num_real_nodes = g.num_nodes()  # 获取图的节点数
        n, d = g.ndata['h'].shape  # 获取节点特征的形状
        real_nodes = list(range(num_real_nodes))
        g.add_nodes(1, {'h': torch.ones(1, d) * 0.02})  # 添加一个虚拟节点
        virtual_src = []
        virtual_dst = []
        for count in range(num_virtual_nodes):  # 添加虚拟节点之间的边
            virtual_node = num_real_nodes + count
            virtual_node_copy = [virtual_node] * num_real_nodes
            virtual_src.extend(real_nodes)
            virtual_src.extend(virtual_node_copy)
            virtual_dst.extend(virtual_node_copy)
            virtual_dst.extend(real_nodes)
        g.add_edges(virtual_src, virtual_dst)
        for ek, ev in g.edata.items():  # 添加虚拟节点的边权重
            ev = torch.cat([ev, torch.zeros(g.num_edges(), 1)], dim=1)
            ev[-num_virtual_nodes * num_real_nodes * 2:, -1] = 1
            g.edata[ek] = ev
        return g

class CMCDB(DGLDataset):
    def __init__(self, cmc_path, phase='train', num_virtual_nodes=0):
        super().__init__(name='CMCDB', )
        self.cmc_data = pd.read_csv(cmc_path).values.tolist()  # 读取数据
        cmc_data_new = []
        for i in self.cmc_data:
            smiles, cmc = i[0], i[1]  # 提取SMILES和目标值
            cmc_data_new.append([smiles, cmc])
        self.cmc_data = cmc_data_new
        smiles = [x[0] for x in self.cmc_data]
        self.mol_grap = {x: process(x, num_virtual_nodes=num_virtual_nodes) for x in smiles if isinstance(x, str)}  # 创建分子图
        self.cmc_data.sort(key=lambda x: x[0])
        random.seed(10)
        random.shuffle(self.cmc_data)

        if phase == 'train':
            self.cmc_data = self.cmc_data[:len(self.cmc_data) * 8 // 10]  # 训练集
        elif phase == 'test':
            self.cmc_data = self.cmc_data[len(self.cmc_data) * 8 // 10 :]  # 测试集

    def __getitem__(self, i):
        one_cmc = self.cmc_data[i]
        mol, cmc = one_cmc
        cmc = cmc / 5.  # 归一化目标值
        mol = self.mol_grap[mol]  # 获取对应的分子图
        return np.array([cmc]), mol

    def __len__(self):
        return len(self.cmc_data)

def rank_descending(input_list):
 sorted_indices = sorted(range(len(input_list)), key=lambda i: input_list[i], reverse=True)
 ranks = [0] * len(input_list)
 for rank, index in enumerate(sorted_indices):
     ranks[index] = rank
 return ranks

class ToyModel(nn.Module):
    def __init__(self, in_dim, hidden_dim, n_classes, ):
        super().__init__()
        self.conv1 = GraphConv(in_dim, hidden_dim)  # 第一层图卷积
        self.conv2 = GraphConv(hidden_dim, hidden_dim)  # 第二层图卷积
        self.classify1 = nn.Linear(hidden_dim, hidden_dim)  # 分类层
        self.classify2 = nn.Linear(hidden_dim, hidden_dim)  # 分类层
        self.classify3 = nn.Linear(hidden_dim, n_classes)  # 输出层
        self.ln1 = nn.LayerNorm(hidden_dim)  # 层归一化
        self.ln2 = nn.LayerNorm(hidden_dim)  # 层归一化
        self.ln3 = nn.LayerNorm(hidden_dim)  # 层归一化
        self.ln4 = nn.LayerNorm(hidden_dim)  # 层归一化

    def forward(self, g):
        h = g.ndata['h'].float()  # 获取节点特征
        h.requires_grad = True
        h1 = F.gelu(self.ln1(self.conv1(g, h)))  # 图卷积和激活函数
        h1 = F.gelu(self.ln2(self.conv2(g, h1)))  # 图卷积和激活函数
        g.ndata['h'] = h1
        hg = dgl.mean_nodes(g, 'h')  # 聚合节点特征
        output = F.relu(self.ln3(self.classify1(hg)))  # 分类层
        output = self.classify3(output)  # 输出结果
        return output

def objective(trial):
    # 定义超参数搜索空间
    hidden_dim = trial.suggest_int('hidden_dim', 64, 512)
    lr = trial.suggest_float('lr', 1e-5, 1e-2, log=True)
    batch_size = trial.suggest_int('batch_size', 16, 128)

    # 创建数据集和数据加载器
    train_db = CMCDB(cmc_path='111.csv', phase='train')
    val_db = CMCDB(cmc_path='111.csv', phase='test')
    train_loader = GraphDataLoader(train_db, batch_size=batch_size, shuffle=True)
    val_loader = GraphDataLoader(val_db, batch_size=batch_size, shuffle=False)

    # 创建模型
    model = ToyModel(in_dim=74, hidden_dim=hidden_dim, n_classes=1)
    loss_fn = nn.MSELoss()
    optim = torch.optim.Adam(model.parameters(), lr)

    # 训练模型
    for epoch in range(50):
        model.train()
        for label, g in train_loader:
            optim.zero_grad()
            res = model(g)
            loss = torch.sqrt(loss_fn(res, label.float()))
            loss.backward()
            optim.step()

    # 验证模型
    model.eval()
    val_loss = 0
    with torch.no_grad():
        for label, g in val_loader:
            res = model(g)
            loss = torch.sqrt(loss_fn(res, label.float()))
            val_loss += loss.item()
    val_loss /= len(val_loader)
    return val_loss

if __name__ == '__main__':
    study = optuna.create_study(direction='minimize')
    study.optimize(objective, n_trials=50)
    print(f"Best trial: {study.best_trial.params}")

    # 保存最佳超参数到JSON文件
    with open('/mnt/data/share/wq/soft/CMC_GCN/0220/code/best_params.json', 'w') as f:
        json.dump(study.best_trial.params, f)

    # 读取最佳超参数
    with open('/mnt/data/share/wq/soft/CMC_GCN/0220/code/best_params.json', 'r') as f:
        best_params = json.load(f)

    hidden_dim = best_params['hidden_dim']
    lr = best_params['lr']
    batch_size = best_params['batch_size']

    train_db = CMCDB(cmc_path='111.csv', phase='train')
    test_db = CMCDB(cmc_path='111.csv', phase='test')
    train_loader = GraphDataLoader(train_db, batch_size=batch_size, shuffle=True)
    val_loader = GraphDataLoader(test_db, batch_size=batch_size, shuffle=False)

    model = ToyModel(in_dim=74, hidden_dim=hidden_dim, n_classes=1)
    loss_fn = nn.MSELoss()
    optim = torch.optim.Adam(model.parameters(), lr)
    epoch = 300
    best_loss = 10
    best_test = 10
    scheduler = torch.optim.lr_scheduler.StepLR(optim, step_size=len(train_loader) * 300, gamma=0.1)

    for e in range(epoch):
        loss_lst = []
        train_res_lst = []
        train_label_lst = []
        for no, i in enumerate(train_loader):
            model.zero_grad()
            label, g = i
            res = model(g)
            loss = torch.sqrt(loss_fn(res, label.float()))
            loss.backward()
            optim.step()
            scheduler.step()
            optim.zero_grad()
            loss_lst.append(loss.item())
            train_res_lst.extend(res.reshape(-1).detach().numpy().tolist())
            train_label_lst.extend(label.reshape(-1).detach().numpy().tolist())

        model.eval()
        loss_lst = []
        val_res_lst = []
        val_label_lst = []
        for no, i in enumerate(val_loader):
            label, g = i
            res = model(g)
            loss = torch.sqrt(F.mse_loss(res, label))
            loss_lst.append(loss.item())
            val_res_lst.extend(res.reshape(-1).detach().numpy().tolist())
            val_label_lst.extend(label.reshape(-1).detach().numpy().tolist())

        sr_test = spearmanr(val_label_lst, val_res_lst)[0]
        a1, a2 = rank_descending(val_res_lst), rank_descending(val_label_lst)
        rank_best30 = len([[a1, a2] for a1, a2 in zip(a1, a2) if a1 <= 30 and a2 <= 30]) / 30
        rank_best10 = len([[a1, a2] for a1, a2 in zip(a1, a2) if a1 <= 10 and a2 <= 10]) / 10
        loss_mean = np.array(loss_lst).mean().round(4)
        r2 = r2_score(val_label_lst, val_res_lst)

        print('test loss', str(e).zfill(3), best_test, loss_mean, 'r2:', r2, 'spearman:', sr_test, 'val rank10:', rank_best10, 'val rank30:', rank_best30)

        torch.save(model.state_dict(), f'lnp_model_simple_{rank_best10}_{rank_best30}.pt')
        model.train()

