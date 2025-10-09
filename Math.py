import Klassen
import Cdept
import NMR
def FunktionelleGruppensubstitutionVERSUCH1GESCHEITERTEINFACHEREIDEE():
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
            # 1 = Anzahl Alkohole / 2 = Ether /
            while zusatuzgruppen[-1] > obergrenze[-1]:
                zusatuzgruppen[0] += 1
                for a in range(len(zusatuzgruppen)-1):
                    if zusatuzgruppen[a] > obergrenze[a]:
                        zusatuzgruppen[a+1] += 1
                        zusatuzgruppen[a] = 0

                # überprüfung ob diese Kombination überhaupt Möglich ist.
                if zusatuzgruppen[1] + zusatuzgruppen[2] == übrigeO:
                    if zusatuzgruppen[0] +(übrigeC-zusatuzgruppen[0])*3 + zusatuzgruppen[1] == übrigeH:
                        if zusatuzgruppen[0] >= Klassen.Molekuelinfo.Carbonsubstitutionsgrad[1] and übrigeC-zusatuzgruppen[0] >= Klassen.Molekuelinfo.Carbonsubstitutionsgrad[3]:
                            möglichkeiten.append(zusatuzgruppen.copy())

            print(gruppenkonfiguration)
            print(möglichkeiten)

            return

    else:

        möglichkeiten = []
        gruppenkonfiguration = [0, 0, 0, 0, 0, 0, 0]

        # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 Keton + Aldehyd / Carbonsäure / Ether

        Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] = len(Klassen.Molekuelinfo.cNRMdaten) - len(
            Klassen.Molekuelinfo.cdeptdaten)

        übrigeC = Klassen.Molekuelinfo.isomere[0] - Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0] - \
                  Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]
        übrigeO = Klassen.Molekuelinfo.isomere[1] - Klassen.Molekuelinfo.oxygeniumsubstitution[1] - \
                  Klassen.Molekuelinfo.oxygeniumsubstitution[2]
        übrigeH = Klassen.Molekuelinfo.isomere[2] - 2 * Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]

        # definierung der schon bekannten gruppen

        gruppenkonfiguration[0] = Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0]
        gruppenkonfiguration[2] = Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]



## Anzahl Alkohole / Anzahl Aldehyde / Anzahl Ketone / Anzahl Carbonsäuren / Anzahl Ester + Carbonsäuren / Anzahl Carbonsäuren + Alkohole
    oxygeniumsubstitution = [0,0,0,0,0,0]
    #Anzahl C / Anzahl CH / Anzahl CH2 / Anzahl CH3 / Anzahl an CH + CH3 welche noch nicht klar sind / weitere CHX gruppen, welche man aber noch nicht zugeordnet hat
    Carbonsubstitutionsgrad= [0,0,0,0,0,0]


def SatzderCsymetrie(gruppenkonfiguration):
    # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 = Aldheyd / 6 = Keton / 7 = Carbonsäure / 8 = Ether / 9 = Ester

    #sum = 0

    #if gruppenkonfiguration[3] > 1:
        #sum += gruppenkonfiguration[3] -1
    #if gruppenkonfiguration[5] > 1:
        #sum += gruppenkonfiguration[5] -1
    #if gruppenkonfiguration[7] > 1:
        #sum += gruppenkonfiguration[7] -1


    #if sum >= Klassen.Molekuelinfo.anzahlcSymetrieelemente:
        #return True
    #else:
        #return False













    # Veralteter Ansatz, da mann nur schauen muss wie viele einstämmige es hat
    sum = 0 #die Summe aller potenzieller Symmetrischen Funktioneller Gruppen
    for a in gruppenkonfiguration:
            if a >= 2:
                sum += a-1

    if gruppenkonfiguration[4] >= 2:
        sum -= a-1
    if gruppenkonfiguration[8] >= 2:
        sum -= a-1


    if sum >= Klassen.Molekuelinfo.anzahlcSymetrieelemente:
        return True
    else:
        return False



def Plausibilitaetskontrolle(C,O,H,gruppenkonfiguration):
    if sum(gruppenkonfiguration) - gruppenkonfiguration[4] - gruppenkonfiguration[8] != C:
        return False
    if gruppenkonfiguration[4] + gruppenkonfiguration[5] + gruppenkonfiguration[6] + gruppenkonfiguration[7] * 2 +gruppenkonfiguration[8] != O:
        return False
    if gruppenkonfiguration[1] + gruppenkonfiguration[2] * 2 + 3 * gruppenkonfiguration[3] + gruppenkonfiguration[4] +gruppenkonfiguration[5] + 2 * gruppenkonfiguration[7] != H:
        return False

    #Das Molekül muss eine ganzzahlige DBÄ haben, ansonsten ist es nicht möglich (Formel : https://de.wikipedia.org/wiki/Doppelbindungs%C3%A4quivalent).
    DBÄ = (gruppenkonfiguration[0]*2 - gruppenkonfiguration[3] - gruppenkonfiguration[4] - gruppenkonfiguration[5] - gruppenkonfiguration[7] + gruppenkonfiguration[1])/2
    if DBÄ < 0:
        return False
    if DBÄ % 1 != 0:
        return False
    # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 = Aldheyd / 6 = Keton / 7 = Carbonsäure / 8 = Ether / 9 = Ester

    # Überprüfung ob die Mindestanzahl an OH bzw. Carboxylgruppe erfüllt ist

    if gruppenkonfiguration[4] + gruppenkonfiguration[7] < Klassen.Molekuelinfo.oxygeniumsubstitution[5]:
        return False





    if Klassen.Molekuelinfo.cSymetrie == False:
        if Klassen.Molekuelinfo.cdeptdaten != None:
            if gruppenkonfiguration[2] != Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]:
                return False

            if gruppenkonfiguration[1] + gruppenkonfiguration[3] +gruppenkonfiguration[5] != Klassen.Molekuelinfo.Carbonsubstitutionsgrad[4]:
                return False
            if gruppenkonfiguration[0] + gruppenkonfiguration[6] + gruppenkonfiguration[7] != Klassen.Molekuelinfo.Carbonsubstitutionsgrad[0]:
                return False




    elif Klassen.Molekuelinfo.cSymetrie == True:
        if Klassen.Molekuelinfo.cdeptdaten != None:
            if gruppenkonfiguration[2] < Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]:
                return False
            if gruppenkonfiguration[1] + gruppenkonfiguration[3] < Klassen.Molekuelinfo.Carbonsubstitutionsgrad[4]:
                return False
        if Klassen.Molekuelinfo.cNRMdaten != None:
            if not SatzderCsymetrie(gruppenkonfiguration):
                return False

            #if not Cdept.CNMRgruppenkonfigurationsplausibilitätskontrolle_beiCSymetrie(gruppenkonfiguration):
                #return False

    else:
        if Klassen.Molekuelinfo.cdeptdaten != None:
            if gruppenkonfiguration[2] < Klassen.Molekuelinfo.Carbonsubstitutionsgrad[2]:
                return False
            if gruppenkonfiguration[1] + gruppenkonfiguration[3] < Klassen.Molekuelinfo.Carbonsubstitutionsgrad[4]:
                return False
        if Klassen.Molekuelinfo.cNRMdaten != None:
            if gruppenkonfiguration[5] + gruppenkonfiguration[6] < Klassen.Molekuelinfo.oxygeniumsubstitution[6]:
                return False
            if Klassen.Molekuelinfo.oxygeniumsubstitution[6] == 0 and gruppenkonfiguration[5] + gruppenkonfiguration[6] != 0:
                return False


    return True


def FunktionelleGruppensubstitution():

    C = Klassen.Molekuelinfo.isomere[0]
    O = Klassen.Molekuelinfo.isomere[2]
    H = Klassen.Molekuelinfo.isomere[1]

    möglichkeiten = []

    # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 = Aldheyd / 6 = Keton / 7 = Carbonsäure / 8 = Ether
    gruppenkonfiguration = [0]*9
    obergrenze = [C,C,round(H/2), round(H/3)+1, NMR.MaximalanzahlOH(), Cdept.MaximalAldehyde(),Cdept.MaximalKetone(),NMR.MaximalanzahlCOOH(),O]
    untergrenze = [0,0,0,0,0,Cdept.MindestanzahlAldehyde(),Cdept.MindestanzahlKetone(),0,0] # Muss modifiziert werden

    while gruppenkonfiguration[-1] <= obergrenze[-1]:
        for a in range(len(gruppenkonfiguration)-1):
            if gruppenkonfiguration[a] <= obergrenze[a]:
                break
            else:
                gruppenkonfiguration[a] = untergrenze[a]
                gruppenkonfiguration[a+1] += 1


        # Plausibilitätskontrolle
        if Plausibilitaetskontrolle(C,O,H,gruppenkonfiguration):
            print("ich habe es geschaft")
            möglichkeiten.append(gruppenkonfiguration.copy())

        gruppenkonfiguration[0] += 1



    print(len(möglichkeiten))

    isomergruppen = [Klassen.Molekuelinfo() for i in range(len(möglichkeiten))]

    for a in range(len(möglichkeiten)):
        isomergruppen[a].gruppenkonfiguration = möglichkeiten[a]
    return isomergruppen


