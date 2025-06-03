import Klassen

def FunktionelleGruppensubstitution():
    if Klassen.Molekuelinfo.cSymetrie == False:

        if Klassen.Molekuelinfo.cdeptdaten != None and Klassen.Molekuelinfo.cNRMdaten != None:

            möglichkeiten = []
            gruppenkonfiguration = [0,0,0,0,0,0,0,0]

            # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 Keton + Aldehyd / Carbonsäure / Ether

            Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] = len(Klassen.Molekuelinfo.cNRMdaten) - len(Klassen.Molekuelinfo.cdeptdaten)

            übrigeC = Klassen.Molekuelinfo.isomere[0] - Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] -  Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]
            übrigeO = Klassen.Molekuelinfo.isomere[1] - Klassen.Molekuelinfo.oxygeniumsubstitution[1] - Klassen.Molekuelinfo.oxygeniumsubstitution[2]
            übrigeH = Klassen.Molekuelinfo.isomere[2] - 2*Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]

            #definierung der schon bekannten gruppen

            gruppenkonfiguration[0] = Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0]
            gruppenkonfiguration[2] = Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]
            gruppenkonfiguration[5] = Klassen.Molekuelinfo.oxygeniumsubstitution[1] + Klassen.Molekuelinfo.oxygeniumsubstitution[2]
            gruppenkonfiguration[6] = übrigeC - Klassen.Molekuelinfo.Carbonsubstitutionsgrad[4]

            übrigeC -= gruppenkonfiguration[6]
            übrigeO -= 2*gruppenkonfiguration[6]
            übrigeH -= gruppenkonfiguration[6]

            if übrigeC < 0 or übrigeO < 0 or übrigeH < 0:
                exit("Fehler bei Math")

            wasserstoffverbrauch = 0

            möglichkeiten = []
            zusatuzgruppen = [0,0,0]
            obergrenze = [übrigeC,übrigeO, übrigeO]
            # 0 = CH1 / CH3 werden nicht benögitgt, da diese impliziert daraus klar sind (C0 und CH2 sind auchs schon klar
            # 1 = Anzahl Alkohole / 2 = Ester /
            while zusatuzgruppen[-1] > obergrenze[-1]:
                zusatuzgruppen[0] += 1
                for a in range(len(zusatuzgruppen)-1):
                    if zusatuzgruppen[a] > obergrenze[a]:
                        zusatuzgruppen[a+1] += 1
                        zusatuzgruppen[a] = 0

                # überprüfung ob diese Kombination überhaupt Möglich ist.
                if zusatuzgruppen[1] + zusatuzgruppen[2] == übrigeO:
                    if zusatuzgruppen[0] +(übrigeC-zusatuzgruppen[0])*3 + zusatuzgruppen[1] == übrigeH:
                        möglichkeiten.append(zusatuzgruppen.copy())

            print(gruppenkonfiguration)
            print(möglichkeiten)

            return


## Anzahl Alkohole / Anzahl Aldehyde / Anzahl Ketone / Anzahl Carbonsäuren / Anzahl Ester + Carbonsäuren / Anzahl Carbonsäuren + Alkohole
    oxygeniumsubstitution = [0,0,0,0,0,0]
    #Anzahl C / Anzahl CH / Anzahl CH2 / Anzahl CH3 / Anzahl an CH + CH3 welche noch nicht klar sind / weitere CHX gruppen, welche man aber noch nicht zugeordnet hat
    Carbonsubstitutionsgrad= [0,0,0,0,0,0]