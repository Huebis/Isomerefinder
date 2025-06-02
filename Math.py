import Klassen

def FunktionelleGruppensubstitution():
    if Klassen.Molekuelinfo.cSymetrie == False:

        möglichkeiten = []
        gruppenkonfiguration = [0,0,0,0,0,0,0,0]

        # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 Keton / Aldehyd / Carbonsäure

        übrigeC = Klassen.Molekuelinfo.isomere[0] - Klassen.Molekuelinfo.Carbonsubstitutionsgrad[1] -  Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2] - Klassen.Molekuelinfo.Carbonsubstitutionsgrad[3]
        übrigeO = Klassen.Molekuelinfo.isomere[1]
        übrigeH = Klassen.Molekuelinfo.isomere[2]
