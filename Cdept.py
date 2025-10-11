# from os import MFD_ALLOW_SEALING

import Klassen

def Cmindestanzahl():
    if Klassen.Molekuelinfo.cNRMdaten != None:
        Klassen.Molekuelinfo.AddCarbonsubstitutionsgrad(5,len(Klassen.Molekuelinfo.cNRMdaten))
        return
    elif Klassen.Molekuelinfo.cdeptdaten != None:
        Klassen.Molekuelinfo.AddCarbonsubstitutionsgrad(5, len(Klassen.Molekuelinfo.cdeptdaten))
        return

    return







def plusminus(messwert, wert , abweichung,):
    if messwert < wert - abweichung:
        return False
    if messwert > wert + abweichung:
        return False
    return True



# Bei weitem noch nicht fertig. Ich möchte, dass das Programm erst nach der MS analyse laufen kann, damit ich mithilfe der genauen Anzahl an C atomen nie CO und C1 und C2 besser bestimmen kann und hoffentlich auch aus dem Gleichungssystem (mit natürlichen Lösugen) alle Möglichkeiten entnehmen kann.
def CdeptundSymetrieIDanalyse():

    for a in range(len(Klassen.Molekuelinfo.Carbonsubstitutionsgrad)):
        Klassen.Molekuelinfo.Carbonsubstitutionsgrad[a] = 0



    messdaten = None
    if Klassen.Molekuelinfo.cNRMdaten != None:
        messdaten = Klassen.Molekuelinfo.cNRMdaten

        if len(Klassen.Molekuelinfo.cNRMdaten) == Klassen.Molekuelinfo.isomere[0]:
            Klassen.Molekuelinfo.cSymetrie = False
            Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] = len(Klassen.Molekuelinfo.cNRMdaten) - len(Klassen.Molekuelinfo.cdeptdaten)
        else:
            Klassen.Molekuelinfo.cSymetrie = True
            Klassen.Molekuelinfo.anzahlcSymetrieelemente = Klassen.Molekuelinfo.isomere[0] - len(Klassen.Molekuelinfo.cNRMdaten)



    if Klassen.Molekuelinfo.cdeptdaten != None:
        messdaten = Klassen.Molekuelinfo.cdeptdaten

        for messwert in messdaten:
            if messwert[1] > 0:
                Klassen.Molekuelinfo.Carbonsubstitutionsgrad[4] += 1
            if messwert[1] < 0:
                Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2] += 1

        # Hier weiterarbeiten das nächste Mal alle CH1 CH3 zuordnen ( wie auch CH2) nachher dann Gleichungssystem lösen (sollte hoffentlich kein Problem sein. Aber mal schaue.
    return

## Anzahl Alkohole / Anzahl Aldehyde / Anzahl Ketone / Anzahl Carbonsäuren / Anzahl Ester + Carbonsäuren / Anzahl Carbonsäuren + Alkohole
    oxygeniumsubstitution = [0,0,0,0,0,0]
    #Anzahl C / Anzahl CH / Anzahl CH2 / Anzahl CH3 / Anzahl an CH + CH3 welche noch nicht klar sind / Anzahl Benzol weitere CHX gruppen, welche man aber noch nicht zugeordnet hat
    Carbonsubstitutionsgrad= [0,0,0,0,0,0]

def MindestanzahlAldehyde():
    if Klassen.Molekuelinfo.cdeptdaten == None:
        return 0
    anzahlAldehyde = 0
    for a in Klassen.Molekuelinfo.cdeptdaten:
        if a[0] >= 190:
            if a[0] <= 250:
                anzahlAldehyde += 1
    return anzahlAldehyde




def MaximalAldehyde():
    if Klassen.Molekuelinfo.cdeptdaten == None:
        return round(Klassen.Molekuelinfo.isomere[2]/2)

    anzahlAldehyde = 0
    for a in Klassen.Molekuelinfo.cdeptdaten:
        if a[0] >= 190:
            if a[0] <= 250:
                anzahlAldehyde += 1
    if anzahlAldehyde == 0 and Klassen.Molekuelinfo.cdeptdaten[0][0] < 160:
        return 0
    else:
        return round(Klassen.Molekuelinfo.isomere[2] / 2) + 1
"""
    if Klassen.Molekuelinfo.cSymetrie == False:
        return anzahlAldehyde
    else:
       
"""

def MindestanzahlKetone():
    if Klassen.Molekuelinfo.cdeptdaten == None:
        return 0
    anzahlKetone = 0
    tempnmr = Klassen.Molekuelinfo.cNRMdaten.copy()

    for a in Klassen.Molekuelinfo.cdeptdaten:
        tempnmr.remove(a[0])


    for a in tempnmr:
        if a >= 190:
            if a <= 250:
                anzahlKetone += 1

    return anzahlKetone

def MaximalKetone():
    if Klassen.Molekuelinfo.cdeptdaten == None:
        return round(Klassen.Molekuelinfo.isomere[2] / 2) + 1

    tempnmr = Klassen.Molekuelinfo.cNRMdaten.copy()

    for a in Klassen.Molekuelinfo.cdeptdaten:
        tempnmr.remove(a[0])

    anzahlKetone = 0
    for a in tempnmr:
        if a >= 190:
            if a <= 250:
                anzahlKetone += 1

    if anzahlKetone == 0 and tempnmr[0] < 160:
        return 0
    else:
        return round(Klassen.Molekuelinfo.isomere[2] / 2) + 1
"""

    if Klassen.Molekuelinfo.cSymetrie == False:
        return anzahlKetone
    else:

"""
def ElementeimBereich(cNMR,anfang,ende):
    sum = 0

    for a in cNMR:
        if a >= anfang:
            if a <= ende:
                sum += 1
    return sum


def CNMRgruppenkonfigurationsplausibilitätskontrolle_beikeinerCSymetrie(gruppenkonfiguration):

    print("Ich überprüfe die CNMR gruppenkonfigurationsplausibiltät (bei CSymetrie")
    tempCNMR = Klassen.Molekuelinfo.cNRMdaten




# Noch nicht fertig

def CNMRgruppenkonfigurationsplausibilitätskontrolle_beiCSymetrie(gruppenkonfiguration):
    tempCNMR = Klassen.Molekuelinfo.cNRMdaten.copy()

    gruppenkonfiguration.pop(4) #Alkohole sind nicht relevant
    #gruppenkonfiguration.pop(8) # Ether sind nicht relevant (aber muss ich behalten für die
    # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 /  4 = Aldheyd / 5 = Keton / 6 = Carbonsäure / Ether


