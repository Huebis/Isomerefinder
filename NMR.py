import Klassen
import copy

def Spektralanalyse_verschwundeD2Opeaks():
    tempnmr = copy.deepcopy(Klassen.Molekuelinfo.nmrdaten)

    for a in Klassen.Molekuelinfo.d20nmrdaten:
        tempnmr.remove(a)

    print(tempnmr)

    mindestanzahlOH = 0
    for a in range(len(tempnmr)-1):
        if tempnmr[a][0] - tempnmr[a+1][0] < 25:
            mindestanzahlOH += 1

    Klassen.Molekuelinfo.oxygeniumsubstitution[5] = mindestanzahlOH
    return



def MaximalanzahlOH():
    if Klassen.Molekuelinfo.nmrdaten == None or Klassen.Molekuelinfo.d20nmrdaten == None:
        return Klassen.Molekuelinfo.isomere[2]
    elif Klassen.Molekuelinfo.oxygeniumsubstitution[5] == 0:
        return 0
    else:
        return Klassen.Molekuelinfo.isomere[2]


def MaximalanzahlCOOH():
    if Klassen.Molekuelinfo.nmrdaten == None or Klassen.Molekuelinfo.d20nmrdaten == None:
        return round(Klassen.Molekuelinfo.isomere[2]/2)
    elif Klassen.Molekuelinfo.oxygeniumsubstitution[5] == 0:
        return 0
    else:
        return round(Klassen.Molekuelinfo.isomere[2]/2)



