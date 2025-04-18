

def IsDoucorrekt(C,O,H):
    dou = (2*C - H)
    if dou < 0:
        return False
    if dou % 2 == 1:
        return False
    return True


def Summenformelerkennung(ms_spektrum,minC):



    molaremasse = 122

    #Zeile 1 Anzahl C, danach Anzahl O und zum Schluss Anzahl H
    möglicheSF = []

    C = minC
    O = (molaremasse - minC*12) // 16
    H = molaremasse - 16*O


    while(True):

        if IsDoucorrekt(C,O,H):
            möglicheSF.append([[C],[O],[H]])

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


    return möglicheSF
