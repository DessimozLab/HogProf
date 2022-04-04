import torch
from torch_geometric.nn import HeteroConv , Linear , TransformerConv , MFConv    
import torch_geometric.transforms as T
import torch.nn.functional as F
import torch_geometric.transforms as T
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.data import HeteroData ,InMemoryDataset
from torch.optim.lr_scheduler import ReduceLROnPlateau
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

import pickle
import numpy as np

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

with open('trainingset_nosectors.pkl' , 'rb')as trainout:
    trainingdata = pickle.loads(trainout.read())
with open('testgset_nosectors' , 'rb')as trainout:
    testingdata = pickle.loads(trainout.read())

print('training samples ' ,len(trainingdata))
print('testing samples' , len(testingdata))

data = trainingdata[0]
class HeteroGCN(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels, num_layers):
        super().__init__()
        
        self.convs = torch.nn.ModuleList()
        self.lins =  torch.nn.ModuleList()
        self.lins2 =  torch.nn.ModuleList()
        
        for i in range(num_layers):
            conv = HeteroConv({
                
                ('phylonodes_up', 'phylolink_up', 'phylonodes_up'):MFConv((-1,-1),  int( hidden_channels/(i+1)) ),
                ('phylonodes_down', 'phylolink_down', 'phylonodes_down'):MFConv((-1,-1),  int( hidden_channels/(i+1)) ),
                ('phylonodes_down', 'phylolink_down_up', 'phylonodes_up'):MFConv((-1,-1),  int( hidden_channels/(i+1)) ),
                ('phylonodes_up', 'phylolink_up_down', 'phylonodes_down'):MFConv((-1,-1),  int( hidden_channels/(i+1)) ),
                ('phylonodes_down', 'informs', 'godnode'):TransformerConv((-1,-1),  int( hidden_channels/(i+1)) ),
                ('phylonodes_up', 'informs', 'godnode'):TransformerConv((-1,-1),  int( hidden_channels/(i+1)) ),
                #('godnode', 'informs', 'phylonodes_down'):TransformerConv((-1,-1),  int( hidden_channels/(i+1)) ),
                #('godnode', 'informs', 'phylonodes_up'):TransformerConv((-1,-1),  int( hidden_channels/(i+1)) ),
            
            } , aggr='sum')
            
            self.convs.append(conv)

            for vectype in  ['phylonodes_up', 'phylonodes_down' , 'sectornode' , 'godnode' ]:
                lin1 = Linear(-1 , int( hidden_channels/(i+1)))
                self.lins.append( lin1 )
                
            print( 'hidden units' , int( hidden_channels/(i+1)) )
            print( 'layer' , i )

        for vectype in ['phylonodes_up', 'phylonodes_down' , 'sectornode' , 'godnode' ]:
            lin2 = Linear(-1 , out_channels)
            self.lins2.append( lin2 ) 
    def forward(self, x_dict, edge_index_dict):
        lins = iter(self.lins)
        for i,conv in enumerate(self.convs):
            x_dict = conv(x_dict , edge_index_dict)
            #x_dict = {key: F.dropout(x , p = .10 , training = self.training ) for key, x in x_dict.items()}
            for key, x in x_dict.items():
                x_dict[key] = next(lins)(x)
        lins2 = iter(self.lins2)
        for key, x in x_dict.items():
            x_dict[key] =  next(lins2)(x)        
        return {key: F.tanh(x) for key, x in x_dict.items()}
        #return x_dict


model = HeteroGCN(hidden_channels=50, out_channels=2, num_layers=5)
model = model.double()

with torch.no_grad():  # Initialize lazy modules.
    out = model(data.x_dict , data.edge_index_dict)

optimizer = torch.optim.Adadelta(model.parameters(), lr=0.1, weight_decay=0)

trainloader = DataLoader( trainingdata , batch_size = 10 , shuffle=True)
testloader = DataLoader(testingdata , batch_size = 10 )

print('start training')
import warnings
lastauc = 0
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = model.to(device)
    for epoch in range(10000):
        model.train()
        losses1=[]
        losses2 =[]
        losses3 =[]

        for i,data in enumerate(trainloader):
            data = data.to(device)
            optimizer.zero_grad()
            out = model(data.x_dict ,data.edge_index_dict)
            loss1 =  F.smooth_l1_loss(out['phylonodes_up'].double(), data['phylonodes_up'].y.double())
            loss2 =  F.smooth_l1_loss(out['phylonodes_down'].double(), data['phylonodes_down'].y.double())
            loss3 =  F.smooth_l1_loss(out['godnode'].double(), data['godnode'].y.double())
            loss = (.5*loss1 + 2.5*loss3 +  .5* loss2 )
            loss.backward()
            optimizer.step()
            losses3.append(float(loss3.to('cpu')))
            losses2.append(float(loss2.to('cpu')))
            losses1.append(float(loss1.to('cpu')))
        print('losses', np.mean(losses1), np.mean(losses2), np.mean(losses3) )
        model.eval()
        aucs = []
        aucsn = []
        for i,testdata in enumerate(testloader):
            testdata = testdata.to(device)
            pred = model(testdata.x_dict ,testdata.edge_index_dict)
            pred = {idx:x.to('cpu') for idx,x in pred.items()}
            truth = testdata['godnode']['y'][:,0].to('cpu').detach().numpy()
            predy =  pred['godnode'][:,0].to('cpu').detach().numpy()
            fpr, tpr, _ = roc_curve(  truth  ,predy )
            aucs.append(auc(fpr, tpr))
            truth = testdata['phylonodes_down']['y'][:,0].to('cpu').detach().numpy()
            predy =  pred['phylonodes_down'][:,0].to('cpu').detach().numpy()
            fpr, tpr, _ = roc_curve(  truth  ,predy )
            aucsn.append(auc(fpr, tpr))
        meanauc = np.mean(aucs)
        print('auc',meanauc)
        if meanauc > lastauc:
            lastauc = meanauc
            print('saving')
            torch.save(model.state_dict(), './phylographnet_script_final_l1.torch')
            print('done')
        meanauc = np.mean(aucsn)
        print('auc_n',meanauc)