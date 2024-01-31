#This file contain the functions related to the Graph Convolutional Network



from rdkit import Chem
from rdkit.Chem import AllChem

import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.data import Data
import torch.nn.functional as F

# The model

class GCN(torch.nn.Module):
    """ Graph convolutionnal network using torch_geometric.nn.GCNconv() wich is 
    the implementation of the work done in this paper : https://arxiv.org/abs/1609.02907
    and with the help of Vaitea
    """

    def __init__(self,nin,nhid1,nhid2,nhid3):

        super().__init__()

        self.conv1 = GCNConv(nin,nhid1)
        self.conv2 = GCNConv(nhid1,nhid2)
        self.conv3 = GCNConv(nhid2,nhid3)
        self.fc = nn.Linear(nhid3, 1)     #output layer, since we predict continuous values we use Linear layer

    def forward(self, data):
        """ data : an object of the class DataBatch that inherits from 'torch_geometric.data.Dataset',
            obtained from a DataLoader object forch_geometric.loader 
            data : contain x, edge_index, batch
            x          : Tensor of the values of every graph of the batch, shape (nAtomsBatch, nFeatures)
            edge_index : Tensor of every edge,  shape (2,nBonds*2)
            batch      : Tensor of The molecule from wich the atom comes from in the batch, shape (nAtomsBatch)
        """
        
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = F.relu(self.conv3(x, edge_index))
        x = global_mean_pool(x, batch)        #pool the data accross every node to reduce output dimention
        x = self.fc(x)

        return x
    
# Graph manipulation
def atom_features(atom):

    
    """ Create an atom features vector. 
    """

    return torch.tensor([
        atom.GetAtomicNum(),    # Atomic number
        atom.GetDegree(),       # Degree
        atom.GetFormalCharge(), # Formal charge
        atom.GetIsAromatic(),   # Aromaticity

        # Features that can be added:
        #atom.GetTotalNumHs(),      # Total number of hydrogens
        #atom.GetImplicitValence(), # Implicit valence
        #atom.GetExplicitValence(), # Explicit valence
        #atom.GetFormalCharge(),    # Formal charge
        #atom.GetMass(),            # Atomic mass
        #atom.GetHybridization(),   # Hybridization state

    ], dtype=torch.float)
    
class Graph:
    """ Create a graph of a molecule from its smile
        NB : target can be None if unknown

    """

    def __init__(self, smile, target,addHs = False) :
        """ smile  : the molecule's smile
            target : the knwn affinity value to target
            addHs  : if set to True : add the Hydrogen on the molecule
        """

        self.smiles = smile
        self.target = target
        
        #create the molecule
        self.smiles_to_mol(addHs)
        
        #create the graph
        self.mol_to_graph(self.target)
        
    
    def smiles_to_mol(self,addHs) :

        # Use MolFromSmiles from RDKit to get a molecule object
        mol = Chem.MolFromSmiles(self.smiles)
        
        # If a valid mol is not returned, set mol as None and exit
        if mol is None:
            self.mol = None
            return

        if addHs :
            self.mol = Chem.AddHs(mol) 
        else : 
            self.mol=mol

    def mol_to_graph(self,target) :
        
        # Get list of atoms in molecule
        self.atoms = self.mol.GetAtoms() #List of objects of class Atom
        
        node_features = torch.stack([atom_features(atom) for atom in self.atoms])
        
        edges = []
        for bond in self.mol.GetBonds():
        
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            edges.append((i, j))
            edges.append((j, i))  # Add both directions

        edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
        
        self.data = Data(x=node_features, edge_index=edge_index, target=target) #edge index = index lié aux noeuds

def df_to_graphs(dataframe, KnownScores = True, addHs = False):
    """ take a dataframe with columns "Smiles" and "Affinity" and return a list of graphs of the molecules
        addHs : if set to true, the molecules will be added his hydrogens
    """

    graphs = [] 

    if KnownScores : 
        for i in range(len(dataframe)) : 
            smile = dataframe['Smile'][i]
            score = float(dataframe['Affinity'][i])
    
            graphs.append(Graph(smile,score,addHs).data) 

    else : 
        for i in range(len(dataframe)) : 
            smile = dataframe['Smile'][i]
            score = None
    
            graphs.append(Graph(smile,score,addHs).data) 

    return graphs