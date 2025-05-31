import Klassen

def Cmindestanzahl():
    if Klassen.Molekuelinfo.cNRMdaten != None:
        Klassen.Molekuelinfo.AddCarbonsubstitutionsgrad(5,len(Klassen.Molekuelinfo.cNRMdaten))
    else:
        Klassen.Molekuelinfo.AddCarbonsubstitutionsgrad(5, len(Klassen.Molekuelinfo.cdeptdaten))
    return





# Bei weitem noch nicht fertig. Ich möchte, dass das Programm erst nach der MS analyse laufen kann, damit ich mithilfe der genauen Anzahl an C atomen nie CO und C1 und C2 besser bestimmen kann und hoffentlich auch aus dem Gleichungssystem (mit natürlichen Lösugen) alle Möglichkeiten entnehmen kann.
def CdeptIDanalyse():
    cdeptdaten = Klassen.Molekuelinfo.cdeptdaten
    for a in cdeptdaten:
        if a[1] < 0:
            Klassen.Molekuelinfo.AddCarbonsubstitutionsgrad(4)
        if a[1] > 0:
            Klassen.Molekuelinfo.AddCarbonsubstitutionsgrad(2)
    if Klassen.Molekuelinfo.cNRMdaten != None:
        cNMRlen = len(Klassen.Molekuelinfo.cNRMdaten)




def cdeptIDanalyse(cdeptdaten, lenofcNMR):
    cdeptID = [0,0,0]
    ceptID[0] = lenofcNMR - len(cdeptdaten)
    for a in cdeptdaten:
        if a[1] < 0:
            cdeptID[1] += 1
        else:
            cdeptID[2] += 1
    return cdeptID
