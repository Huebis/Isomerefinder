import Klassen

#Wie in der Maturaarbeit beschrieben, ist die Funktion aufgrund komischer MS-Spektren nicht gebrauch

#ganzes File wird momentan nicht gebraucht



def IsDoucorrekt(C,O,H):
    dou = (2*C - H)
    if dou < 0:
        return False
    if dou % 2 == 1:
        return False
    return True


def Summenformelerkennung():
    molaremasse = Klassen.Molekuelinfo.msmainpeak

    #Zeile 1 Anzahl C, danach Anzahl O und zum Schluss Anzahl H
    möglicheSF = []

    C = sum(Klassen.Molekuelinfo.Carbonsubstitutionsgrad)
    O = (molaremasse - C*12) // 16
    H = molaremasse - 16*O


    while(True):

        if IsDoucorrekt(C,O,H):
            möglicheSF.append([C,O,H])

        if O > 0:
            O -= 1
            H += 16
        else:

            C += 1
            H -= 12
            if H < 0:
                break
            O = H//16
            H -= O*16


    Summenformelranking(möglicheSF)

    return


def Summenformelranking(möglichesummenformeln):
    molaremasse = Klassen.Molekuelinfo.msmainpeak
    msspektrum = Klassen.Molekuelinfo.msdata
    oneaddpeak = 0
    twoaddpeak = 0
    for a in range(len(msspektrum)-1,-1,-1):
        if msspektrum[a][0] == molaremasse:
            peak = msspektrum[a][1]
            if a < len(msspektrum) - 1:
                oneaddpeak = msspektrum[a+1][1]
                if a < len(msspektrum) - 2:
                    twoaddpeak = msspektrum[a + 2][1]

            break

    accuracy0 = 4
    Summenformel = [0,0,0]



    for a in range(len(möglichesummenformeln)):
        prozentualoneaddpeak = oneaddpeak*100/peak
        prozentualtwoaddpeak = twoaddpeak * 100 / peak
        C13wahrscheinlichkeit = möglichesummenformeln[a][0]*1.1
        H2wahrscheinlichkeit = möglichesummenformeln[a][2]*0.016
        O17wahrscheinlichkeit = möglichesummenformeln[a][1]*0.04
        O18wahrscheinlichkeit = möglichesummenformeln[a][1]*0.2
        accuracy = abs(1- (C13wahrscheinlichkeit+H2wahrscheinlichkeit+O17wahrscheinlichkeit)/prozentualoneaddpeak)

        if prozentualtwoaddpeak != 0:
            accuracy += abs(1-((C13wahrscheinlichkeit*(C13wahrscheinlichkeit-1.1) + H2wahrscheinlichkeit*(H2wahrscheinlichkeit-0.016) + O18wahrscheinlichkeit + O17wahrscheinlichkeit*(O17wahrscheinlichkeit-0.04)+C13wahrscheinlichkeit*H2wahrscheinlichkeit+C13wahrscheinlichkeit*O17wahrscheinlichkeit + H2wahrscheinlichkeit*O17wahrscheinlichkeit)/prozentualtwoaddpeak))

        if accuracy0 > accuracy:
            accuracy0 = accuracy
            Summenformel[0] = möglichesummenformeln[a][0]
            Summenformel[1] = möglichesummenformeln[a][1]
            Summenformel[2] = möglichesummenformeln[a][2]

    Klassen.Molekuelinfo.Replaceisomere(Summenformel)
    return Summenformel



