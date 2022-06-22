import torch
import torch.utils.data
import torch.nn.functional as F
from torch_geometric.nn import  GATConv
import torch.nn as nn
from NN.Attention_Layer import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

nHiddenUnits = 8
nHeads = 8 
dropout = 0.4 
alpha = 0.2


class GAT_transformer(torch.nn.Module):
    """
    GAT model
    param d_input: number of input features
    param d_output: dim of output vector
    param s_max: dim for attention layer
    
    return: classification result
    """
    def __init__(self,d_input:int,d_output:int,s_max:int):
        super(GAT_transformer, self).__init__()
        self.d_input = d_input
        self.d_output = d_output
        self.s_max = s_max
        self.gat1 = GATConv(self.d_input, out_channels=nHiddenUnits,
                            heads=nHeads, concat=True, negative_slope=alpha,
                            dropout=dropout, bias=True)
        self.gat2 = GATConv(70, out_channels=d_output,
                            heads=nHeads, concat=False, negative_slope=alpha,
                            dropout=dropout, bias=True)
        self.transformer = SelfAttention(self.s_max)
        self.norm1 = nn.LayerNorm(64)
        self.norm = nn.LayerNorm(6)
        self.dropout = torch.nn.Dropout()

    def forward(self, x, edge_index):        
        x,gat_attn1 = self.gat1(x, edge_index,return_attention_weights =True)
        x = self.norm1(x)
        x = F.elu(x)        
        x_t = self.transformer(x)
        x_t = self.norm(x_t)
        x,gat_att2 = self.gat2(torch.cat((x,x_t),dim=1), edge_index,return_attention_weights =True)
        return F.log_softmax(x, dim=1)
