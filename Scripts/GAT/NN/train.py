import numpy as np
from tqdm import tqdm
import torch.nn.functional as F
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import additional_function.modules as modules


def train_batches(model,optimizer,dataloader,device,lf='nll'):
    """
    train 1 epoch
    
    param model: Neural Network
    param optimizer: optimizer
    param datalodaer: dataloader
    param lf: loss function nll or cross entropy
    
    return: Tuple(output, mean epoch loss, mean epoch accuracy) 
    """
    epoch_loss = []
    epoch_acc = []
    for batch in dataloader:
        batch = batch.to(device)
        x, edge_index = batch.x, batch.edge_index
        optimizer.zero_grad()
        output = model(x,edge_index)
        if lf == 'nll':
            loss = F.nll_loss(output, batch.y)
        else:
            loss = F.cross_entropy(output, batch.y)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch_loss.append(loss.item())
        epoch_acc.append(modules.accuracy(output, batch.y).item())
    
    return output,np.mean(epoch_loss),np.mean(epoch_acc)

def val_batches(model,dataloader,device,lf = 'nll'):
    """
    train 1 epoch
    
    param model: Neural Network
    param optimizer: optimizer
    param datalodaer: dataloader
    param lf: loss function nll or cross entropy
    
    return: Tuple(output, mean epoch loss, mean epoch accuracy) 
    """
    epoch_loss = []
    epoch_acc = []
    outputs = np.array([])
    for batch in dataloader:
        batch = batch.to(device)
        x, edge_index = batch.x, batch.edge_index
        output = model(x,edge_index)
        if lf == 'nll':
            loss = F.nll_loss(output, batch.y)
        else:
            loss = F.cross_entropy(output, batch.y)
        epoch_loss.append(loss.item())
        epoch_acc.append(modules.accuracy(output, batch.y).item())
        outputs = np.append(outputs,output.max(1)[1].detach().cpu().numpy())
    return outputs,np.mean(epoch_loss),np.mean(epoch_acc)

def predict_batches(model,dataloader,permutation_index:np.array,true_index:np.array,device:str)->np.array:
    """
    prediction
    
    param model: Neural Network
    param datalodaer: dataloader
    param permutation_index: indexes after permutation
    param true_index: true_index 
    device: cuda or cpu
    
    return: classification results
    """
    outputs = np.array([])
    patient_index = [int(np.where(permutation_index==i)[0]) for i in true_index]
    for batch in tqdm(dataloader):
        batch = batch.to(device)
        x, edge_index = batch.x, batch.edge_index
        output = model(edge_index,x)
        outputs = np.append(outputs,output.max(1)[1].detach().cpu().numpy())
    outputs_valid = outputs[patient_index] 
    return outputs_valid


def train(model,optimizer,
          num_epoch,train_loader,val_loader,
          device:str,lf:str = 'entropy',test_loader = None, verbose = None):
    """
    training procedure
    
    param model: Neural Network
    param optimizer: optimizer
    param train_loader: train dataloader
    param val_loader: val dataloader
    param device: cuda or cpu
    param lf: loss function nll or cross entropy
    param test_loader: test dataloader
    verbose: verbose
    
    return: classification result, total accuracy train,total accuracy test,total loss train,total loss val
    """
    
    total_acc_train = []
    total_loss_train = []
    total_acc_test = []
    total_loss_test = []
    min_loss = 100
    for epoch in range(num_epoch):
        model.train()
        out,epoch_loss,epoch_acc = train_batches(model,optimizer,
                                                 train_loader,device,lf)
        model.eval()
        out,epoch_loss_val,epoch_acc_val = val_batches(model,val_loader,device,lf)
        total_acc_train.append(epoch_acc),total_loss_train.append(epoch_loss)
        total_acc_test.append(epoch_acc_val),total_loss_test.append(epoch_loss_val)
        if verbose:
            print('epoch: ',epoch,'train acc: ',epoch_acc,'train loss: ', epoch_loss,
                'val acc: ',epoch_acc_val,'val loss: ', epoch_loss_val)
        if (epoch%10  == 0) &(test_loader!=None):
            model.eval()
            out,targ,test_loss,test_acc = val_batches(model,val_loader,device,lf)
            print('epoch: ',epoch,'train acc: ',epoch_acc,'train loss: ', epoch_loss,
                'val acc: ',epoch_acc_val,'val loss: ', epoch_loss_val,
                 'test acc: ',test_acc,'test loss: ', test_loss) 
        if epoch_loss_val < min_loss:
            min_loss = epoch_loss_val
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve == 50:
                print('Early stopping!' )
                return model,total_acc_train,total_acc_test,total_loss_train,total_loss_test
    return model,total_acc_train,total_acc_test,total_loss_train,total_loss_test