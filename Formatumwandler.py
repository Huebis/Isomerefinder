from rdkit import Chem
from rdkit.Chem import AddHs


#Lesen eines files im SMILE FORMAT in ein 1D Array
def Readfile():
    with open('C7H6O2.smi', 'r', encoding='utf-8') as file:
        zeilen = file.readlines()  # Reads the entire file


    isomerarray = []
    for eintrag in zeilen:
        isomerarray.append(eintrag[:-2])
    return isomerarray


def AddWasserstoff(isomerearray):

    for t in range(len(isomerearray)):
        molekulobject = Chem.MolFromSmiles(isomerearray[t])
        molekulobject = AddHs(molekulobject)
        isomerearray[t] = Chem.MolToSmiles(molekulobject)
    return isomerearray