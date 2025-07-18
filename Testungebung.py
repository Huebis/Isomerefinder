from rdkit import Chem
from rdkit.Chem import Draw
def test():
    # SMILES-String eingeben
    smiles = 'C(=C1(C(#C(C(#C)))))(C(#C(C(#C(C2)))))(C1(C(=C(C=2))))'  # Beispiel: Benzol mit expliziten Doppelbindungen

    # Molekül erzeugen
    mol = Chem.MolFromSmiles(smiles)
    print(Chem.MolToSmiles(mol))

    mol = Chem.AddHs(mol)
    # Bild erzeugen
    img = Draw.MolToImage(mol)

    # Bild anzeigen
    img.show()


#[[0, [3, 1], [1, 2]], [0, [3, 0], [1, 3]], [0, [1, 0], [3, 3]], [0, [3, 2], [1, 1]]]


if (True and False) or True:
    print("jdjdjdd")