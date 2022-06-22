from sklearn.metrics import accuracy_score
import NN.train as training
from sklearn.model_selection import train_test_split
from abc import ABC, abstractmethod
from torch_geometric.data import Data
from tqdm import tqdm
from preprocessing.prep_func import *
from anndata import AnnData
from additional_function import modules
from preprocessing import prep_func
import numpy as np
import torch
import json
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


class ModelTrainer(ABC):
    """
    Abstract base class for  trainer.
    This class is an operator which does not cache any states inside.
    """
    def __init__(self, model, args=None):
        self.model = model
        self.id = 0
        self.args = args

    def set_id(self, trainer_id):
        self.id = trainer_id

    @abstractmethod
    def get_model_params(self):
        pass

    @abstractmethod
    def set_model_params(self, model_parameters):
        pass

    @abstractmethod
    def train(self, train_data, device, args=None):
        pass

    @abstractmethod
    def test(self, test_data, device, args=None):
        pass



def save_json(data:dict,name:str):
    """
    save json file
    """
    with open(name, 'w') as fp:
        json.dump(data, fp)
                
def loader_creation(data:dict,device:str,args:dict,shuffle=True):
    """
    data loader creation
    """
    node_features = torch.from_numpy(data['X']).float()
    target = args['target']
    try:
        data[target]
    except:
        target = False
        print('no target')
    else:
        print('target detected')
        labels = data[target]
        labels = torch.LongTensor(labels)
    edge_index,attr = modules.scipysparse2torchsparse(data['adj'])    
    if target:
        dat = Data(x=node_features, edge_index=edge_index, y=labels)
    else:
        dat = Data(x=node_features, edge_index=edge_index)
    dat = dat.to(device)
    cd = modules.ClusterData(dat,num_parts = args['NumParts'])
    loader = modules.ClusterLoader(cd,batch_size = args['BatchSize'],shuffle=shuffle)
    return loader,cd.perm.detach().cpu().numpy()



class GcnTrainer(ModelTrainer):
    """
    main training and prediction module
    """
    def get_model_params(self):
        return self.model.cpu().state_dict()

    def set_model_params(self, model_parameters):
        print("set_model_params")
        self.model.load_state_dict(model_parameters)

    def train(self, train_data:AnnData,  args:dict):
        """
        training procedure
        
        param training_data: AnnData object with training dataset
        param args: dictionary with params of model
        
        return: accuracy on train and val, weights saved in dictionary from args
        """
        model = self.model
        device = args['device']
        device = torch.device(device)
        model.to(device)
        model.train()
        # test_data = None
        # try:
        #     test_data = self.test_data
        # except:
        #     pass
        if args['client_optimizer'] == "sgd":
            optimizer = torch.optim.SGD(model.parameters(), lr=args['lr'])
        else:
            optimizer = torch.optim.Adam(model.parameters(), lr=args['lr'])
        epoch_accs_train = []
        epoch_accs_val = []
        idx_train, idx_test = train_test_split(train_data.obs.index,
                                               train_size=args['training_size'])
        tdata = train_data[(train_data.obs.index.isin(idx_train))]
        val = train_data[(train_data.obs.index.isin(idx_test))]
        tdata = prep_func.graph_pp(tdata,num_comps=100,num_neighbors=30, bbknn=False)
        val = prep_func.graph_pp(val,num_comps=100,num_neighbors=30, bbknn=False)
        gdata_train = prep_func.dictthat(tdata)
        gdata_val = prep_func.dictthat(val)
        cl_train, _ = loader_creation(gdata_train,device,args,True)
        cl_val, _ = loader_creation(gdata_val,device,args,True)
        min_loss = 100
        epochs_no_improve = 0
        for epoch in tqdm(range(args['epochs'])):
            model.train()
            out,epoch_loss,epoch_acc = training.train_batches(model,optimizer,
                                                 cl_train,device,args['lf'])
            model.eval()
            out,epoch_loss_val,epoch_acc_val = training.val_batches(model,cl_val,device,args['lf'])
            if epoch%10 ==0:
                print('train epoch loss: ',epoch_loss,'train epoch accuracy:', epoch_acc)
                print('val epoch loss: ',epoch_loss_val,'val epoch accuracy:', epoch_acc_val)
            epoch_accs_train.append(epoch_acc)
            epoch_accs_val.append(epoch_acc_val)
            if epoch_loss_val < min_loss:
                min_loss = epoch_loss_val
                epochs_no_improve = 0
            else:
                epochs_no_improve += 1
                if epochs_no_improve == 30:
                    print('Early stopping!' )
                    return epoch_accs_train,epoch_accs_val
        if args['save']:
            model_weights = model.cpu().state_dict()
            torch.save(model_weights, args['PATH'])
            data_size = len(tdata)
            optimizers_params = optimizer.state_dict()
            params = dict(zip(['set_size','args'],[data_size,args]))
            print('params',params)
            save_json(params,args['params_name'])
        return epoch_accs_train,epoch_accs_val

    def test(self, test_data:AnnData, args:dict):
        """
        prediction on test data
        
        param test_data: test dataset
        param args: args
        
        return: model, prediction
        """
        print("----------test--------")
        device = 'cpu'#args['device']
        device = torch.device(device)
        model = self.model
        model.eval()
        model.to(device)
        test = prep_func.graph_pp(test_data,num_comps=100,num_neighbors=30, bbknn=False)        
        gdata_train = prep_func.dictthat(test)
        test_reset = test_data.obs.reset_index()
        test_reset_indx = np.array(test_reset.index)
        cl_test, perm_test = loader_creation(gdata_train,device,args,False)
        prediction = training.predict_batches(model,cl_test,perm_test,test_reset_indx,device)
        if args['target'] in test_reset.columns:
            labels_val = test.obs[args['target']]
            labels_ind = labels_val[test_reset_indx]
            print('test accuracy',accuracy_score(prediction,labels_ind))
        return model,prediction,labels_ind        
            
def transform_list_to_tensor(model_params_list:list)->torch.Tensor:
    '''
    create torch.Tensor from list of lists
    '''
    for k in model_params_list.keys():
        model_params_list[k] = torch.from_numpy(np.asarray(model_params_list[k])).float()
    return model_params_list


def transform_tensor_to_list(model_params:torch.Tensor)->list:
    '''
    create list from torch.Tensor
    '''
    for k in model_params.keys():
        model_params[k] = model_params[k].detach().numpy().tolist()
    return model_params

