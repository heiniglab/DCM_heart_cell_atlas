label_encoder = {'TTN' : 0, 'LMNA': 1,'PKP2': 2,'RBM20': 3,'mutation negative' : 4}
print(label_encoder)

args = {
    'ctype_encoder':label_encoder,
    'Gene' : 'Gene',
    'device': 'cuda',
    'input_dime':3000,
    'num_classes':5,
    's_max':41,
    'Patient':'Patient',
    'training_size':0.8,
    'lf':'nll',
    'client_optimizer':'adam',
    'epochs':600,
    'lr':1e-5,
    'target':'target',
    'NumParts':4000,
    'BatchSize':256,
    'save':True,
    'PATH':'model_weights/model2_weights.ckpt',
    'params_name':'model_weights/model2_params.json'}


