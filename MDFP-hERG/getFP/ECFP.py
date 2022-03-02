#! /usr/bin/python
# coding: utf-8
from rdkit import Chem
from rdkit.Chem import AllChem
naismiles = open('herg.smiles', 'r', encoding='utf-8')
cnt = 0
def add_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for i in range(atoms):
        mol.GetAtomWithIdx(i).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(i).GetIdx()))
    return mol
cwd = os.getcwd()
path = cwd+'/'
ot = open("path+otecfp.csv", 'w+')

for naismile in naismiles:
    try:
        smile = naismile.split(' ')[0]
        COMPID = naismile.split(' ')[1][:-1]
        m1 = Chem.MolFromSmiles(smile)
        add_atom_index(m1)
        bi = {}
        fp = AllChem.GetMorganFingerprint(m1, 2, bitInfo=bi)
        cnt += 1
        print(cnt)
        print(COMPID, ',', fp.GetNonzeroElements(), '\n', file=ot)
    except:
        print('error', '\n', file=ot)
        continue
