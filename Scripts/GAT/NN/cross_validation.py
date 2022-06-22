from NN.GAT import GAT_transformer
import NN.train as training
from preprocessing import prep_func
from anndata import AnnData
from NN.FML import loader_creation
from additional_function import modules
import torch
import numpy as np
from sklearn.model_selection import train_test_split
from torch_geometric.data import Data
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

NumParts = 4000
BatchSize = 256
NumParts = 4000 # num sub-graphs
Device = 'cuda' # if no gpu, `Device='cpu'`
WeightDecay=5e-4
fastmode = False # if `fastmode=False`, report validation
nHiddenUnits = 8
nHeads = 8 # number of attention heads
nEpochs = 5000
dropout = 0.4 # applied to all GAT layers
alpha = 0.2 # alpha for leaky_relu
patience = 100 # epochs to beat
clip = None # set `clip=1` to turn on gradient clipping
s_max = 41
rs=829197 # random_seed
Device = 'cuda'


def cross_validation(adata:AnnData,args:dict,verbose=False):
    """
    Leave one out cross validation across all patients

    param adata: anndata object
    param args: params for GAT model
    param verbose: verbose intermediate results of training procedure
    """
    
    device = torch.device(args['device'])
    ctype_encoder = args['ctype_encoder']
    cv_by_patients = []
    cv_patients = []
    cv_proba = []
    accuracy_by_patient = {}
    patients = set(adata.obs[args['Patient']])
    # iterating by patients
    for patient in list(patients):
        print(patient)
        training_patients = set(adata.obs[args['Patient']]).difference(set([patient]))
        indx = adata.obs[adata.obs[args['Patient']] == patient].index
        train_length = 0
        val_length = 0
        adata.obs['target'] = adata.obs[args['Gene']].map(ctype_encoder)
        num_unique_classes = len(set(adata.obs['target']))
        # check that all classes in train, test, val with splitting according patients
        while  (train_length<args['num_classes']) or (val_length<args['num_classes']):
            idx_train, idx_test = train_test_split(adata[adata.obs[args['Patient']].isin(training_patients)].obs.index, train_size=args['training_size'])
            tdata = adata[(adata.obs.index.isin(idx_train))]
            test = adata[(adata.obs.index.isin(idx_test))]
            train_length = len((np.unique(tdata.obs['target'],return_counts = True)[0]))
            val_length = len((np.unique(test.obs['target'],return_counts = True)[0]))
        #generating list with final indexes for current patient and validation sample for graph creation
        # then only test samples considering for accuracy measures    
        final_test_ind = list(idx_test)+list(indx)
        #Graph creation
        tdata = prep_func.graph_pp(tdata,num_comps=100,num_neighbors=30, bbknn=False)
        test = prep_func.graph_pp(test,num_comps=100,num_neighbors=30, bbknn=False)
        tdata.obs['target'] = tdata.obs[args['Gene']].map(ctype_encoder)
        test.obs['target'] = test.obs[args['Gene']].map(ctype_encoder)
        gdata_train = prep_func.dictthat(tdata)
        gdata_test  = prep_func.dictthat(test)
        #graph dataloader 
        dataloader_train, _ = loader_creation(gdata_train,device,args,True)
        dataloader_val, _ = loader_creation(gdata_test,device,args,True)
        model = GAT_transformer(args['input_dime'],args['num_classes'],s_max).to(device)
        optimizer = torch.optim.Adam(model.parameters(),
                                    lr=args['lr'],weight_decay=WeightDecay)
        #trainig procedure                                    
        epoch_accs_train = []
        epoch_accs_val = []
        min_loss = 100
        epochs_no_improve = 0
        assert len(set(tdata.obs[args['Patient']]).intersection(set([patient]))) == 0, 'intersection between train patients and current patient'
        assert len(set(test.obs[args['Patient']]).intersection(set([patient]))) == 0, 'intersection between train patients and current patient'
        for epoch in range(args['epochs']):
            model.train()
            out,epoch_loss,epoch_acc = training.train_batches(model,optimizer,
                                                dataloader_train,device,args['lf'])
            model.eval()
            out,epoch_loss_val,epoch_acc_val = training.val_batches(model,dataloader_val,device,args['lf'])
            if epoch%10 == 0 and verbose == True:
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
        print("best train accuracy is: ",np.max(epoch_accs_train))
        print("best test accuracy is: ",np.max(epoch_accs_val))
        # test graph generation
        test = adata[(adata.obs.index.isin(final_test_ind))]
        test = prep_func.graph_pp(test,num_comps=100,num_neighbors=30, bbknn=False)
        test.obs['target'] = test.obs[args['Gene']].map(ctype_encoder)
        gdata_test  = prep_func.dictthat(test)
        node_features_val = torch.from_numpy(gdata_test['X']).float()
        labels_val = gdata_test['target']
        labels_val = torch.LongTensor(labels_val)
        edge_index_val,attr_val = modules.scipysparse2torchsparse(gdata_test['adj'])
        d_test = Data(x=node_features_val,edge_index=edge_index_val,y=labels_val)
        d_test = d_test.to(args['device'])
        x, edge_index = d_test.x, d_test.edge_index
        #test prediction
        predictions = model(x, edge_index) 
        ind = (np.where(np.array(gdata_test[args['Patient']]) == patient)[0])
        labels_ind = labels_val[ind]
        # subsetting only prediction from current patient
        predictions = predictions[ind]
        #accuracy measurments
        accuracy_patient = modules.accuracy(predictions,labels_ind).item()
        print('-----'*20)
        print('accuracy for patient :',patient,' is ',accuracy_patient)
        print('-----'*20)
        predictions_int = predictions.detach().cpu().max(1)[1]
        cv_by_patients.append(predictions_int)
        cv_proba.append(torch.exp(predictions.detach().cpu().max(1)[0]))
        cv_patients.append(patient)
        accuracy_by_patient[patient] = accuracy_patient
    return accuracy_by_patient,cv_patients,cv_proba,cv_by_patients





    
    
    