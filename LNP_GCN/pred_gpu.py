import torch
import torch.nn as nn
import torch.nn.functional as F
import dgl
from dgl.nn.pytorch import GraphConv
from dgllife.utils import CanonicalAtomFeaturizer, mol_to_bigraph
from rdkit import Chem
import pandas as pd

# 定义ToyModel类（确保与训练时的模型结构一致）
class ToyModel(nn.Module):
    def __init__(self, in_dim, hidden_dim, n_classes):
        super().__init__()
        self.conv1 = GraphConv(in_dim, hidden_dim)
        self.conv2 = GraphConv(hidden_dim, hidden_dim)
        self.classify1 = nn.Linear(hidden_dim, hidden_dim)
        self.classify2 = nn.Linear(hidden_dim, hidden_dim)
        self.classify3 = nn.Linear(hidden_dim, n_classes)
        self.ln1 = nn.LayerNorm(hidden_dim)
        self.ln2 = nn.LayerNorm(hidden_dim)
        self.ln3 = nn.LayerNorm(hidden_dim)
        self.ln4 = nn.LayerNorm(hidden_dim)

    def forward(self, g):
        h = g.ndata['h'].float()
        h.requires_grad = True
        h1 = F.gelu(self.ln1(self.conv1(g, h)))
        h1 = F.gelu(self.ln2(self.conv2(g, h1)))
        g.ndata['h'] = h1
        hg = dgl.mean_nodes(g, 'h')
        output = F.relu(self.ln3(self.classify1(hg)))
        output = self.classify3(output)
        return output

# 加载模型权重
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = ToyModel(in_dim=74, hidden_dim=290, n_classes=1).to(device)
model.load_state_dict(torch.load('lnp_model_simple_0.5_0.5.pt'))
model.eval()

# 处理新的SMILES数据
def process_smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    node_enc = CanonicalAtomFeaturizer()
    g = mol_to_bigraph(m, True, node_enc, None, False)
    return g

# 预测函数
def predict(smiles):
    g = process_smiles(smiles)
    g = g.to(device)
    with torch.no_grad():
        res = model(g)
    return res.item()

# 读取CSV文件并进行预测
def predict_from_csv(csv_file):
    df = pd.read_csv(csv_file)
    predictions = []
    for smiles in df.iloc[:, 0]:  # 假设第一列是SMILES
        prediction = predict(smiles)
        predictions.append(prediction)
    df['prediction'] = predictions
    return df

# 示例：预测test.csv中的SMILES数据
result_df = predict_from_csv('test.csv')
result_df.to_csv('predictions.csv', index=False)
print("Predictions saved to predictions.csv")
