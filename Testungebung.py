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


if False:
    raise ValueError("Ungültiger Wert!")


test = [[0, [3, 1], [1, 11]], [0, [3, 0], [1, 11]], [1, [1, 10], [1, 7], [1, 3]], [1, [1, 7], [1, 5], [1, 2]], [1, [1, 9], [1, 6], [1, 14]], [1, [1, 13], [1, 10], [1, 3]], [1, [1, 12], [1, 9], [1, 4]], [1, [1, 3], [1, 8], [1, 2]], [2, [1, 15], [1, 7]], [2, [1, 4], [1, 6]], [3, [1, 2], [1, 5]], [3, [1, 0], [1, 1]], [5, [1, 6]], [5, [1, 5]], [6, [1, 4]], [6, [1, 8]], [7, [1, 17]], [7, [1, 16]]]

print(test[11])