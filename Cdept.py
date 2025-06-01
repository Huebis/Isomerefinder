from os import MFD_ALLOW_SEALING

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
def CdeptIDanalyse():

    for a in Klassen.Molekuelinfo.Carbonsubstitutionsgrad:
        Klassen.Molekuelinfo.Carbonsubstitutionsgrad[a] = 0

    for a in Klassen.Molekuelinfo.oxygeniumsubstitution:
        Klassen.Molekuelinfo.oxygeniumsubstitution[a] = 0


    messdaten = None
    if Klassen.Molekuelinfo.cNRMdaten != None:
        messdaten = Klassen.Molekuelinfo.cNRMdaten

    elif Klassen.Molekuelinfo.cdeptdaten != None:
        messdaten = Klassen.Molekuelinfo.cdeptdaten

    if messdaten != None
        for messwert in messdaten:

            if plusminus(messwert, 175,5):
                Klassen.Molekuelinfo.oxygeniumsubstitution[3] += 1
                # Zuerst muss ausgeschlossen werden, dass die Carbonsäure eigenlich ein COH=O-H ist.
                if Klassen.Molekuelinfo.isomere[0] != 1 and Klassen.Molekuelinfo.isomere[1] != 2 and Klassen.Molekuelinfo.isomere[2] != 2:
                    Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] += 1

            if plusminus(messwert, 200,5):
                Klassen.Molekuelinfo.oxygeniumsubstitution[1] += 1
                if Klassen.Molekuelinfo.isomere[0] != 1 and Klassen.Molekuelinfo.isomere[1] != 1 and Klassen.Molekuelinfo.isomere[2] != 2: # gleiches hier, nur beim Aldehyd
                    Klassen.Molekuelinfo.Carbonsubstitutionsgrad[1] += 1

            if plusminus(messwert, 210,5):
                Klassen.Molekuelinfo.oxygeniumsubstitution[2] += 1
                Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] += 1

            #if plusminus(messwert, 128,15):
                #print("test")

    if Klassen.Molekuelinfo.cdeptdaten != None:
        # Hier weiterarbeiten das nächste Mal alle CH1 CH3 zuordnen ( wie auch CH2) nachher dann Gleichungssystem lösen (sollte hoffentlich kein Problem sein. Aber mal schaue.
    return

## Anzahl Alkohole / Anzahl Aldehyde / Anzahl Ketone / Anzahl Carbonsäuren / Anzahl Ester + Carbonsäuren / Anzahl Carbonsäuren + Alkohole
    oxygeniumsubstitution = [0,0,0,0,0,0]
    #Anzahl C / Anzahl CH / Anzahl CH2 / Anzahl CH3 / Anzahl an CH + CH3 welche noch nicht klar sind / Anzahl Benzol weitere CHX gruppen, welche man aber noch nicht zugeordnet hat
    Carbonsubstitutionsgrad= [0,0,0,0,0,0]



