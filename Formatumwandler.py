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

def AddWasserstoff2():
    isomere = 'O1C2CC3(OC3C#C1)C2'

    def AnzahlAtombindungen(char):
        if char == 'O':
            return 2
        if char == 'C':
            return 4
        if char == 'c':
            return 3
    def Bindungen(char):
        if char == '=':
            return 2
        if char == '#':
            return 3
        return 1  #Wenn nicht ein = oder ein # kommen muss vorhin eine Zahl, ein Atom oder eine Klammer kommen, dies impliziert dann direkt, dass es automatisch ein einfach Bindung ist


    for a in range (len(isomere)):
        if isomere[a] == 'C' or isomere[a] == 'c' or isomere[a] == 'O':
            anzahlBindungen = AnzahlAtombindungen(isomere[a])
            print("Hallo")
            print("a:" + str(a))

            bindungenohneWasserstoff = 0
            if a > 0:
                bindungenohneWasserstoff += Bindungen(isomere[a-1])



            for b in range(1,len(isomere)-a,1):
                if isomere[a+b].isdigit():
                    bindungenohneWasserstoff += 1


                if not isomere[a+b].isdigit():
                    c = 0
                    if isomere[a+b] == '(':
                        bindungenohneWasserstoff += Bindungen(isomere[a + b+1])
                        c += 1
                        anzahlklammern = 1 # nach rechts offene Klammer ist +1, nach links offene -1
                        while(True):
                            c += 1
                            if isomere[a+b+c] == '(':
                                anzahlklammern += 1
                            if isomere[a+b+c] == ')':
                                anzahlklammern -= 1
                                if anzahlklammern == 0:
                                    if not isomere[a+b+c+1] == '(':
                                        c += 1
                                        break
                                    else:
                                        bindungenohneWasserstoff += Bindungen(isomere[a + b + c + 2])
                                        c += 3
                                        anzahlklammern = 1  # nach rechts offene Klammer ist +1, nach links offene -1



                    if isomere[a+b+c] == ')':
                        break
                    else:
                        bindungenohneWasserstoff += Bindungen(isomere[a+b+c])
                        break

            print(anzahlBindungen -bindungenohneWasserstoff)





