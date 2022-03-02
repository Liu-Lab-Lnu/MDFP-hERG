from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Descriptors3D, rdMolDescriptors
import rdkit


def fix_valence_charge(mol):
    for atom in mol.GetAtoms():
        explicitValence = 0
        for bond in atom.GetBonds():
            explicitValence = explicitValence + bond.GetBondTypeAsDouble()
        atom.SetFormalCharge(0)
        if atom.GetSymbol() == 'N' and explicitValence > 3:
            atom.SetFormalCharge(int(explicitValence - 3))
        elif atom.GetSymbol() == 'O' and explicitValence > 2:
            atom.SetFormalCharge(int(explicitValence - 2))
        elif atom.GetSymbol() == 'C' and explicitValence > 4:
            atom.SetFormalCharge(int(explicitValence - 4))
    return Chem.SanitizeMol(mol)


suppl = Chem.SDMolSupplier('1-all.sdf', sanitize=False)
mols = [x for x in suppl]
[fix_valence_charge(x) for x in mols]


ot1 = open('ot1.csv', 'a+')
ot2 = open('ot2.csv', 'a+')
ot3 = open('ot3.csv', 'a+')
ot4 = open('ot4.csv', 'a+')
ot5 = open('ot5.csv', 'a+')
ot6 = open('ot6.csv', 'a+')
ot7 = open('ot7.csv', 'a+')
ot8 = open('ot8.csv', 'a+')
ot9 = open('ot9.csv', 'a+')
ot10 = open('ot10.csv', 'a+')
ot11 = open('ot11.csv', 'a+')
ot12 = open('ot12.csv', 'a+')
ot13 = open('ot13.csv', 'a+')
ot14 = open('ot14.csv', 'a+')
ot15 = open('ot15.csv', 'a+')


cnt = 0


cnt = 0
print("PBF_fps", file=ot1)
for mol in mols:
    PBF_fps = rdMolDescriptors.CalcPBF(mol)
    cnt += 1
    print(cnt, PBF_fps, file=ot2, sep=",")


cnt = 0
print("PMI1_fps", file=ot2)
for mol in mols:
    PMI1_fps = rdMolDescriptors.CalcPMI1(mol)
    cnt += 1
    print(cnt, PMI1_fps, file=ot3, sep=",")


cnt = 0
print("PMI2_fps", file=ot3)
for mol in mols:
    PMI2_fps = rdMolDescriptors.CalcPMI2(mol)
    cnt += 1
    print(cnt, PMI2_fps, file=ot4, sep=",")


cnt = 0
print("PMI3_fps", file=ot4)
for mol in mols:
    PMI3_fps = rdMolDescriptors.CalcPMI3(mol)
    cnt += 1
    print(cnt, PMI3_fps, file=ot5, sep=",")


cnt = 0
print("NPR2_fps", file=ot5)
for mol in mols:
    NPR2_fps = rdMolDescriptors.CalcNPR2(mol)
    cnt += 1
    print(cnt, NPR2_fps, file=ot7, sep=",")


cnt = 0
print("RG_fps", file=ot6)
for mol in mols:
    RG_fps = rdMolDescriptors.CalcRadiusOfGyration(mol)
    cnt += 1
    print(cnt, RG_fps, file=ot8, sep=",")


cnt = 0
print("SF_fps", file=ot7)
for mol in mols:
    SF_fps = rdMolDescriptors.CalcInertialShapeFactor(mol)
    cnt += 1
    print(cnt, SF_fps, file=ot9, sep=",")


cnt = 0
print("Aspher_fps", file=ot8)
for mol in mols:
    Aspher_fps = rdMolDescriptors.CalcAsphericity(mol)
    cnt += 1
    print(cnt, Aspher_fps, file=ot11, sep=",")


cnt = 0
print("Spher_fps", file=ot9)
for mol in mols:
    Spher_fps = rdMolDescriptors.CalcSpherocityIndex(mol)
    cnt += 1
    print(cnt, Spher_fps, file=ot12, sep=",")


ot1.close()
ot2.close()
ot3.close()
ot4.close()
ot5.close()
ot6.close()
ot7.close()
ot8.close()
ot9.close()
