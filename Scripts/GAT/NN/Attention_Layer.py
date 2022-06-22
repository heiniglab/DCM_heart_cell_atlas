import torch.utils.data
import torch.nn as nn
import torch.nn.functional as F
import math

class Encoder(nn.Module):
    """
    Attentiom layer
    """
    def __init__(self, dim_Q, dim_K, dim_V, s_max, d_model, num_heads, ln=False, skip=True):
        super(Encoder, self).__init__()
        self.dim_V = dim_V
        self.num_heads = num_heads
        self.skip = skip
        self.s_max = s_max
        self.d_model = d_model
        self.fc_q = nn.Linear(dim_Q, d_model)
        self.fc_k = nn.Linear(dim_K, d_model)
        self.fc_v = nn.Linear(dim_K, d_model)
        if ln:
            self.ln0 = nn.LayerNorm(d_model)
            self.ln1 = nn.LayerNorm(d_model)
        self.ff = nn.Sequential(
        nn.Linear(d_model, 4 * d_model),
        nn.ReLU(),
        nn.Linear(4 * d_model, d_model))
        self.fc_rep = nn.Linear(s_max, d_model)
        
    def forward(self, Q, K):
        Q = self.fc_q(Q)
        K, V = self.fc_k(K), self.fc_v(K)
        dim_split = self.d_model // self.num_heads
        Q_ = torch.cat(Q.split(dim_split, -2), 0)
        K_ = torch.cat(K.split(dim_split, -2), 0)
        V_ = torch.cat(V.split(dim_split, -2), 0)
        Q_ = Q_.reshape(1,Q_.shape[0],Q_.shape[1])
        K_ = K_.transpose(-2,-1)
        K_ = K_.reshape(1,K_.shape[0],K_.shape[1])
        A = torch.softmax(Q_.bmm(K_)/math.sqrt(self.d_model), dim=-1)
        V_ = V_.reshape(1,V_.shape[0],V_.shape[1])
        A_1 = A.bmm(V_)
        O = torch.cat((A_1).split(Q.size(0), 0), 2)
        O = torch.cat((Q_ + A_1).split(Q.size(0), 0), 2) if getattr(self, 'skip', True) else \
             torch.cat((A_1).split(Q.size(0), 0), 2)
        O = O if getattr(self, 'ln0', None) is None else self.ln0(O)
        O = O + self.ff(O)
        O = O if getattr(self, 'ln1', None) is None else self.ln1(O)
        O = F.pad(O, (0, self.s_max- O.shape[-1]), 'constant', 0)
        O = self.fc_rep(O)
        O = O.squeeze()
        return O

class SelfAttention(nn.Module):
    def __init__(self, s_max, dim_in=64, dim_out=6, num_heads=2, ln=True, skip=True):
        super(SelfAttention, self).__init__()
        self.Encoder = Encoder(dim_in, dim_in, dim_in, s_max, dim_out, num_heads, ln=ln, skip=skip)
    def forward(self, X):
        return self.Encoder(X, X)