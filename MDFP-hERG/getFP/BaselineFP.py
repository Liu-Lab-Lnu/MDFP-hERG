import numpy,os
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.linear_model import Lasso
from sklearn.ensemble import GradientBoostingRegressor
import pandas as pd

# helper function for generating MDFP+
def _getStatsData(data):
  ave = numpy.mean(data)
  stdv = numpy.std(data)
  med = numpy.median(data)
  return [ave, stdv, med]

# function to generate the MDFP+
def getMDFPplus(m):
  fp = []
  fp.append(m.GetNumHeavyAtoms())
  fp.append(AllChem.CalcNumRotatableBonds(m))
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#7]')))) # nitrogens
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#8]')))) # oxygens
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#9]')))) # fluorines
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#15]')))) # phosphorous
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#16]')))) # sulfurs
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#17]')))) # chlorines
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#35]')))) # bromines
  fp.append(len(m.GetSubstructMatches(Chem.MolFromSmarts('[#53]')))) # iodines
  return fp

# paths
cwd = os.getcwd()
path = cwd+'/'
smi_file = "herg.smiles"
smi = []
h_name = []

for line in open(path+smi_file,'r'):
    line = line.rstrip().split()
    smi.append(line[0])
    h_name.append(line[1])

# Here, the times series extracted from the MD trajectories are read in.
# For this example, the arrays are filled with random numbers.
tot_ene = [1, 2, 3, 4]
tot_lj = [1, 2, 3, 4]
tot_coulomb = [1, 2, 3, 4]
intra_ene = [1, 2, 3, 4]
intra_lj = [1, 2, 3, 4]
intra_coulomb = [1, 2, 3, 4]
rgyr = [1, 2, 3, 4]
sasa = [1, 2, 3, 4]

# generate the MDFP+
this_data = [tot_ene, tot_lj, tot_coulomb, intra_ene, intra_lj, intra_coulomb, rgyr, sasa]
head = ['HeavyAtoms','RotatableBonds','nitrogens','oxygens','fluorines','phosphorous','sulfurs','chlorines','bromines','iodines']
datain = []
findd = []
cnt = 1
for m in smi:
    try:
        mol = Chem.MolFromSmiles(m)
        mdfp = getMDFPplus(mol)
        datain.append(mdfp)
        print(m)
        cnt = cnt+1
    except:
        findd.append(cnt)
        continue
print(findd)
print(cnt)
herg_fp = pd.DataFrame(columns = head,data=datain)
herg_fp.to_csv("./otbfp.csv",encoding="gbk")

