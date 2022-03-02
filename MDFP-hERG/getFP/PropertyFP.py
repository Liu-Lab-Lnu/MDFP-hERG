from descriptastorus import raw
import numpy
from descriptastorus import MolFileIndex
from descriptastorus.descriptors import rdNormalizedDescriptors
from rdkit import Chem
import logging
import pandas as pd

index = MolFileIndex.MakeSmilesIndex("./herg.smiles", "herg", hasHeader=False, smilesColumn=0, nameColumn=1) # 5 ge canshu ketiao

# make the normalized descriptor generator
generator = rdNormalizedDescriptors.RDKit2DNormalized()
columns = generator.columns # list of tuples:  (descriptor_name, numpytype) ...
descriptor_name = [x[0] for x in columns] # list name ketiao

# example for converting a smiles string into the values
def rdkit_2d_normalized_features(smiles: str):
    # n.b. the first element is true/false if the descriptors were properly computed
    results = generator.process(smiles)
    processed, features = results[0], results[1:]
    if processed is None:
        logging.warning("Unable to process smiles %s", smiles)
    # if processed is None, the features are are default values for the type
    return features # list zhi

r = raw.MakeStore(columns, 939, "herg_store") # (columns, shujuliang+2 , kuming)

for i in range(0,938): # Data size
    if i == 703:
        continue
    index_mol = index.getMol(i)
    index_name = index.getName(i)
    features = rdkit_2d_normalized_features(index_mol)
    r.putRow(i, features)
    print(index_name, " ", i)

r = raw.RawStore("herg_store")
head = r.colnames
datain = []
index_name = []

for i in range(0,938): # Data size
    if i == 703:
        continue
    putout = r.get(i)
    datain.append(putout)
    index_name.append(index.getName(i))

herg_fp = pd.DataFrame(columns = head,index = index_name,data=datain)
herg_fp.to_csv("./PropertyFP.csv",encoding="gbk")
