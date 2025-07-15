from rdkit import Chem
from rdkit.Chem import Draw

# SMILES-String eingeben
smiles = 'C(#C(O(O1)))(C(#C1))'  # Beispiel: Benzol mit expliziten Doppelbindungen

# Molekül erzeugen
mol = Chem.MolFromSmiles(smiles)
print(Chem.MolToSmiles(mol))

mol = Chem.AddHs(mol)
# Bild erzeugen
img = Draw.MolToImage(mol)

# Bild anzeigen
img.show()


#[[0, [3, 1], [1, 2]], [0, [3, 0], [1, 3]], [0, [1, 0], [3, 3]], [0, [3, 2], [1, 1]]]
