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


ot = open('ot.csv', 'a+')
print("mol", end = "", file=ot)


for mol in mols:
    PBF_fps = rdMolDescriptors.CalcPBF(mol)
    print(",", PBF_fps, file=ot, end=",", sep="")


for mol in mols:
    PMI1_fps = rdMolDescriptors.CalcPMI1(mol)
    print(PMI1_fps, file=ot, end=",")


for mol in mols:
    PMI2_fps = rdMolDescriptors.CalcPMI2(mol)
    print(PMI2_fps, file=ot, end=",")


for mol in mols:
    PMI3_fps = rdMolDescriptors.CalcPMI3(mol)
    print(PMI3_fps, file=ot, end=",")


for mol in mols:
    NPR2_fps = rdMolDescriptors.CalcNPR2(mol)
    print(NPR2_fps, file=ot, end=",")


for mol in mols:
    RG_fps = rdMolDescriptors.CalcRadiusOfGyration(mol)
    print(RG_fps, file=ot, end=",")


for mol in mols:
    SF_fps = rdMolDescriptors.CalcInertialShapeFactor(mol)
    print(SF_fps, file=ot, end=",")


for mol in mols:
    Aspher_fps = rdMolDescriptors.CalcAsphericity(mol)
    print(Aspher_fps, file=ot, end=",")


for mol in mols:
    Spher_fps = rdMolDescriptors.CalcSpherocityIndex(mol)
    print(Spher_fps, file=ot, end="")


print("\n", end = "", file=ot)

ot.close()
