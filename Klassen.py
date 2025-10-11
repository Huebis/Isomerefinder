import subprocess
import random
import copy
from logging import raiseExceptions
import math
import itertools

import Testcase

import rdkit.Chem.AllChem
from rdkit import Chem
from rdkit.Chem import Draw

from rdkit.Chem import FindAllSubgraphsOfLengthN


class Molekuelinfo():

    # Anzahl C / Anzahl 0/ Anzahl H
    isomere = [0,0,0]
    #0 Anzahl Alkohole /1 Anzahl Aldehyde /2 Anzahl Ketone /3 Anzahl Carbonsäuren /4 Anzahl Ester + Carbonsäuren /5 Anzahl Carbonsäuren + Alkohole /6 Anzahl Aldehyde + Ketone
    oxygeniumsubstitution = [0]*7
    #Anzahl C / Anzahl CH / Anzahl CH2 / Anzahl CH3 / Anzahl an CH + CH3 welche noch nicht klar sind / weitere CHX gruppen, welche man aber noch nicht zugeordnet hat
    Carbonsubstitutionsgrad= [0,0,0,0,0,0]

    msmainpeak = None
    msdata = None
    cNRMdaten = None
    cdeptdaten = None
    nmrdaten = None
    cSymetrie = None
    anzahlcSymetrieelemente = None
    d20nmrdaten = None
    funktionelleGruppenkonfiguration = None

    @classmethod
    def Replaceisomere(cls, list):
        cls.isomere = list

    @classmethod
    def Addoxygeniumsubstitution(cls, number, count=1):
        cls.oxygeniumsubstitution[number] += count

    @classmethod
    def AddCarbonsubstitutionsgrad(cls, number, count=1):
        cls.Carbonsubstitutionsgrad[number] += count
        #cls.Carbonsubstitutionsgrad[5] -= count


    @classmethod
    def Printinfo(cls):
        print("isomere:", cls.isomere)
        print("oxygenium:", cls.oxygeniumsubstitution)
        print("Carbonsubstitutionsgrad:", cls.Carbonsubstitutionsgrad)
        print("Csymetry:", cls.cSymetrie)

    @classmethod
    def VorbereitungNMRdatenfürCalcheuristik(cls):
        nmrsignalbundel = []


        daten = copy.deepcopy(cls.d20nmrdaten)

        temp = [daten[0]]
        daten.pop(0)
        while daten != []:
            if abs(temp[-1][0] - daten[0][0]) < 25: #Wert anpassen falls nötig
                temp.append(daten[0])
                daten.pop(0)
            else:
                nmrsignalbundel.append(temp)
                temp = [daten[0]]
                daten.pop(0)
        nmrsignalbundel.append(temp)

        """
        print("NMRBunde")
        print(nmrsignalbundel)
        """

        herzunterschiedeinbundels = []

        for nmrbundel in nmrsignalbundel:
            temp = []
            for a in range(len(nmrbundel)):
                for b in range(a+1,len(nmrbundel)):
                    if abs(nmrbundel[a][0] -nmrbundel[b][0]) < 25:
                        temp.append(abs(nmrbundel[a][0] -nmrbundel[b][0]))

            herzunterschiedeinbundels.append(temp)

        cls.gruppierted20nmrdaten = [[None] for a in range(len(nmrsignalbundel))]

        #zuerst bekommt jede Gruppe die passende Chemsiche Verschiebung

        for a in range(len(nmrsignalbundel)):
            cls.gruppierted20nmrdaten[a][0] = nmrsignalbundel[a][round(len(nmrsignalbundel[a])/2)][1]

        #jetzt wird noch geschaut welche Signale sich mit welchen Koppeln lassen
        for a in range(len(nmrsignalbundel)):
            for b in range(len(nmrsignalbundel)):
                if b != a:
                    if any(abs(x-y) <= 0.11 for x in herzunterschiedeinbundels[a]  for y in herzunterschiedeinbundels[b]): #wert noch anpassen falls nötig
                        cls.gruppierted20nmrdaten[a].append(b)









    def __init__(self):
        # Ein mutables Objekt wie eine Liste wird in allen Instanzen geteilt, wenn es im Konstruktor nicht explizit initialisiert wird
        # self.CarbonID = [0] * 5
        self.gruppenkonfiguration = None
        self.gruppierted20nmrdaten = None

    def EntwicklungIsomerelist(self):

        # 0 = C / 1 = CH / 2 = CH2 / 3 = CH3 / 4 = OH / 5 = Aldheyd / 6 = Keton / Carbonsäure / Ether (gruppenkonfiguration)

        summenformelstring = ""
        if self.gruppenkonfiguration[0] != 0:
            summenformelstring += "C" + str(self.gruppenkonfiguration[0])
        if self.gruppenkonfiguration[1] != 0:
            summenformelstring += "N" + str(self.gruppenkonfiguration[1])
        if self.gruppenkonfiguration[2] != 0:
            summenformelstring += "O" + str(self.gruppenkonfiguration[2])
        if self.gruppenkonfiguration[3] != 0:
            summenformelstring += "S(val=1)" + str(self.gruppenkonfiguration[3])
        if self.gruppenkonfiguration[4] != 0:
            summenformelstring += "P(val=1)" + str(self.gruppenkonfiguration[4])
        if self.gruppenkonfiguration[5] != 0:
            summenformelstring += "F" + str(self.gruppenkonfiguration[5])
        if self.gruppenkonfiguration[6] != 0:
            summenformelstring += "I(val=2)" + str(self.gruppenkonfiguration[6])
        if self.gruppenkonfiguration[7] != 0:
            summenformelstring += "Cl" + str(self.gruppenkonfiguration[7])
        if self.gruppenkonfiguration[8] != 0:
            summenformelstring += "Br(val=2)" + str(self.gruppenkonfiguration[8])
        if self.gruppenkonfiguration[9] != 0:
            summenformelstring += "H(val=2)" + str(self.gruppenkonfiguration[9])

        result = subprocess.run('java -jar /home/eliahh/Workspace/Matura/test/MAYGEN/target/MAYGEN-1.8.jar -f "'+summenformelstring +'" -setElements -v -t -smi -m -o /home/eliahh/Workspace/Matura/test/Isomergruppen', shell=True, capture_output=True, text=True)

        print(result)

        def test(self):
            print("isomere:", self.isomere)
            print("oxygenium:", self.oxygeniumsubstitution)
            print("Carbonsubstitutionsgrad:", self.Carbonsubstitutionsgrad)
            print("gruppenkonfiguration", self.gruppenkonfiguration)




class individuum():

    def Anzahloffeneverbindungen(self,Atom):
        anzahlmöglicheverbindungen = self.Anzahlverbindungen(Atom[0])

        anzahlverbrauchteverbindungen = 0
        for a in range(1, len(Atom), 1):
            anzahlverbrauchteverbindungen += Atom[a][0]

        return anzahlmöglicheverbindungen - anzahlverbrauchteverbindungen

    def AnzahlverbindungenzwischenzweiAtomen(self,Atom1, posititionAtom2):
        for a in range(1, len(Atom1), 1):
            if Atom1[a][1] == posititionAtom2:
                return Atom1[a][0]
        return 0

    def Deletverbindung(self, positionatom1, positionatom2):
        for a in range(1, len(self.molekularstruktur[positionatom1])):
            if self.molekularstruktur[positionatom1][a][1] == positionatom2:
                self.molekularstruktur[positionatom1].pop(a)
                break
        for a in range(1, len(self.molekularstruktur[positionatom2])):
            if self.molekularstruktur[positionatom2][a][1] == positionatom1:
                self.molekularstruktur[positionatom2].pop(a)
                break
        return

    def Creatverbindung(self, positionatom1, positionatom2, wertigkeitderVerbindung):
        if self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstruktur[positionatom1],positionatom2) == 0:
            self.molekularstruktur[positionatom1].append([wertigkeitderVerbindung, positionatom2])
            self.molekularstruktur[positionatom2].append([wertigkeitderVerbindung, positionatom1])
        else:
            wertigkeitderaltenVerbindung = self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstruktur[positionatom1],positionatom2)
            self.Deletverbindung(positionatom1,positionatom2)
            self.molekularstruktur[positionatom1].append([wertigkeitderVerbindung+wertigkeitderaltenVerbindung, positionatom2])
            self.molekularstruktur[positionatom2].append([wertigkeitderVerbindung+wertigkeitderaltenVerbindung, positionatom1])

    def Isindirektverbunden(self,positionatom1, positionatom2):
        #überprüft ob zwei Atome direkt oder indirket via anderer Atome miteinander verbunden ist

        schondurchlofeneelementepositionen = [positionatom1]
        nochzudurchlaufendeelementpositionen = [positionatom1]

        while nochzudurchlaufendeelementpositionen != []:
            atom = self.molekularstruktur[nochzudurchlaufendeelementpositionen[0]]
            nochzudurchlaufendeelementpositionen.pop(0)
            for a in range(1, len(atom)):
                if atom[a][1] == positionatom2:
                    return True
                if not atom[a][1] in schondurchlofeneelementepositionen:
                    schondurchlofeneelementepositionen.append(atom[a][1])
                    nochzudurchlaufendeelementpositionen.append(atom[a][1])

        return False


    def isligit(self, anzahlverbindungenüberspringen = False, anazahlAtomeüberspringen = False, topologischzusammenhängendüberspringen = False):

        def AlleAtomehabenalleverbindungenverbraucht():
            for atom in self.molekularstruktur:
                if self.Anzahloffeneverbindungen(atom) != 0:
                    #print("Error nicht alle Atome haben die richtige Anzahl an Verbindungen" + str(atom))
                    return False
            return True

        def AlleAtomekommeninderrichtigneAnzahlvor():
            elemente = [0 for a in range(len(self.elemente))]
            for atom in self.molekularstruktur:
                elemente[atom[0]] += 1
            if elemente == self.elemente:
                return True

            #print("Error es kommen nicht alle Atome in der richtigen Anzahl vor")
            return False



        def Topologischzusammenhängend():
            überprüfteatome = [False for a in range(len(self.molekularstruktur))]
            verbindungenderAtomeüberprüfen = []

            überprüfteatome[0] = True
            verbindungenderAtomeüberprüfen.append(0)

            while verbindungenderAtomeüberprüfen != []:
                momentanprüfend = verbindungenderAtomeüberprüfen[0]
                verbindungenderAtomeüberprüfen.pop(0)
                for a in range(1,len(self.molekularstruktur[momentanprüfend])):
                    if not überprüfteatome[self.molekularstruktur[momentanprüfend][a][1]]:
                        überprüfteatome[self.molekularstruktur[momentanprüfend][a][1]] = True
                        verbindungenderAtomeüberprüfen.append(self.molekularstruktur[momentanprüfend][a][1])

            for a in überprüfteatome:
                if not a:
                    #print("Error das Atom ist nicht zusammenhängend")
                    return False

            return True





        if not AlleAtomehabenalleverbindungenverbraucht() and not anzahlverbindungenüberspringen :
            return False
        if not AlleAtomekommeninderrichtigneAnzahlvor() and not anazahlAtomeüberspringen:
            return False
        if not Topologischzusammenhängend() and not topologischzusammenhängendüberspringen:
            return False
        return True

    def Iszyklisch(self, elementnummer):
        #return True wenn in einem zyklischen Verbindung ist und gibt False wenn nicht
        if self.molekularstruktur[elementnummer][0] >= self.elementgruppengrenzen[3]: #Aufgrund der einfachbindung ist es nicht möglich
            return False

        for a in range(1,len(self.molekularstruktur[elementnummer])-1,1):
            wertikeitderverbindung = self.molekularstruktur[elementnummer][a][0]
            atom2 = self.molekularstruktur[elementnummer][a][1]

            self.Deletverbindung(elementnummer, atom2)
            bool = self.isligit(True, True, False)

            self.Creatverbindung(elementnummer, atom2, wertikeitderverbindung)

            if bool:
                return True

        return False

    def Cyclogrösse(self, position): #Zuerst muss überprüft werden ob es überhaupt Cyclisch ist
        #Die Funktion gibt zurück wie gross die grösste Cyclo ist (int)


        entscheidungen = [[1,position]]
        cyclos = []
        while entscheidungen != []:
            #print("Entscheidungen")
            #print(entscheidungen)
            #print(len(entscheidungen))
            if len(self.molekularstruktur[entscheidungen[-1][1]]) > entscheidungen[-1][0]:
                positionneuesatom = self.molekularstruktur[entscheidungen[-1][1]][entscheidungen[-1][0]][1]
                neuesatom = self.molekularstruktur[positionneuesatom]

                if position == positionneuesatom:
                    cyclos.append(copy.deepcopy(entscheidungen))
                    entscheidungen[-1][0] += 1
                elif len(neuesatom) == 2: # bzw. hat nur Verbindung zu einem Atom
                    entscheidungen[-1][0] += 1
                elif not any(entscheidung[1] == positionneuesatom for entscheidung in entscheidungen):
                    entscheidungen.append([1,positionneuesatom])
                else:
                    entscheidungen[-1][0] += 1

            else:
                entscheidungen.pop(-1)
                if entscheidungen != []:
                    entscheidungen[-1][0] += 1


        if cyclos == []:
            raise Exception("Fehler bei Cyclogrösse, ist gar kein Cyclo")
        # Da bei entscheidungen alles abgelaufen wird, werden pro Cyclo auch zwei Einträge gemacht (einmal von der einen Seite her und einmal von der anderen)
        # es kann sein, dass ein Atom in mehr als einem Cyclo ist, dann muss der wichtigste genommen werden.
        #da ein Aromat das wichtigste ist, wird dies zuerst mal herausgefiltert
        else:
            for möglichkeiten in cyclos:
                if len(möglichkeiten) == 6:
                    #abklären ob es ein Aromat ist
                    bool = True
                    nächsteverbindung = 2
                    for tuple in möglichkeiten:
                        if not self.molekularstruktur[tuple[1]][tuple[0]][0] == nächsteverbindung:
                            bool = False
                            break
                        if nächsteverbindung == 2:
                            nächsteverbindung -= 1
                        else:
                            nächsteverbindung += 1
                    if bool:
                        return 0

            # Falls es kein Aromat ist, wird einfach die höchte Cyclo verbindung genommen für eine möglichst genaue aproximation

            max = 0

            for a in cyclos:
                if len(a) > max:
                    max = len(a)
            return max


    def PositionendernachbarnwelcheeineWasserstoffkopplungeingehen(self,position):
        # Die Funktion gibt wieder welche Positionen bzw. Atome eine H-Kopplung haben müssen mit dem inputAtom und dies immer als Tuple. Der zweite Wert gibt an, wie oft die Kopplung im NMR mindestens vorkommen muss

        kopplungen = []


        """
        #zuerst muss überprüft werden ob es sich bei der Position um ein Atom in einem Aromat handelt, da der Weg dann ein bisschen spezieller ist
        if self.Iszyklisch(position):
            if self.Cyclogrösse(position) == 0: # Bedeutung es ist ein Aromat

        """
        # 0 = C / 1 = CH / 2 = CH2 / 3 = Keton / 4 = Ether / 5 = CH3 / 6 = OH / 7 = Aldehyd / 8 = Carbonsäure
        # self.elementgruppengrenzen = [0,1,2,5]

        #Normallfall falls es kein Aromat ist
        atom = self.molekularstruktur[position]

        for a in range(1,len(atom)):
            if self.molekularstruktur[atom[a][1]][0] == self.elementgruppengrenzen[3]: # ist ein CH3
                kopplungen.append([atom[a][1],3])
            elif self.molekularstruktur[atom[a][1]][0] == self.elementgruppengrenzen[2]: # ist ein CH2
                kopplungen.append([atom[a][1],1]) # Aufgrund von Stereoisomer
            elif self.molekularstruktur[atom[a][1]][0] == self.elementgruppengrenzen[1]: # ist ein CH
                kopplungen.append([atom[a][1],1])
            elif self.molekularstruktur[atom[a][1]][0] == self.elementgruppengrenzen[3] + 2: # ist ein C=O(H)
                kopplungen.append([atom[a][1],1])


        if kopplungen != []:
            return kopplungen
        return None




    def CalcHeuristik(self):
        def BestrafungKetonAlkoholverbindungen(): # 0 = C / 1 = CH / 2 = CH2 / 3 = Keton / 4 = Ether / 5 = CH3 / 6 = OH / 7 = Aldehyd / 8 = Carbonsäure
            #Das Program kann ein Keton an ein Alkohol setzen, dies ist aber nicht möglich, da dies eine Carbonsäure geben wird und diese Anzahl schon genau bestimmt ist.
            #Abzug von 1 für jede solche KetonAlkoholverbindung

            anzahlKetonAlkoholverbindungen = 0

            for atom in self.molekularstruktur:
                if atom[0] == self.elementgruppengrenzen[2] +1:
                    for a in range(1,len(atom)):
                        if atom[a][1] == self.elementgruppengrenzen[3] + 1:
                            anzahlKetonAlkoholverbindungen += 1
            return anzahlKetonAlkoholverbindungen


        def BestrafungEthertransformationzuKeton_Aldehyd():
            anzahlfalscherEther = 0 #für Jeden Ether mit einer Doppelbindung gibt es einen Punkt (darf eigenlich nicht sein
            for atom in self.molekularstruktur:
                if atom[0] == self.elementgruppengrenzen[2]+2:
                    if len(atom) == 2:
                        anzahlfalscherEther += 1

            return anzahlfalscherEther


        def MSpeaküberprüfung():

            def MassezählermitStartpunkt(startpunkt, atomebeidenenegestoptwerdensoll):
                nochzubearbeitetenelemente = [startpunkt]
                schonbearbeiteteelemente = atomebeidenenegestoptwerdensoll

                masse = 0

                while nochzubearbeitetenelemente != []:
                    atom = self.molekularstruktur[nochzubearbeitetenelemente[0]]
                    schonbearbeiteteelemente.append(nochzubearbeitetenelemente[0])
                    nochzubearbeitetenelemente.pop(0)
                    for a in range(1,len(atom),1):
                        if not (atom[a][1] in schonbearbeiteteelemente):
                            nochzubearbeitetenelemente.append(atom[a][1])
                    masse += self.elementmassen[atom[0]]

                return masse





            # zuerst wird der alpha cleavage untersucht bei Ketonen für jede passende Alpha cleavage gibt es einen Minus Punkt und für jede Klivage, die nicht existiert gibt es 1 Punkt
            mstreffer = 0

            for position,atom in enumerate(self.molekularstruktur):
                if atom[0] == self.elementgruppengrenzen[2] + 1: # ist ein Keton
                    if not self.Iszyklisch(position) and len(atom) > 2:
                        masse1 = MassezählermitStartpunkt(atom[1][1], [position]) + 28 # 28 für das hinzufügen des Ketons
                        masse2 = MassezählermitStartpunkt(atom[2][1],[position]) + 28  # 28 für das hinzufügen des Ketons

                        for mswert in Molekuelinfo.msdata:
                            if mswert[0] == masse1:
                                mstreffer -= 5
                                masse1 = 0
                            if mswert[0] == masse2:
                                mstreffer -= 5
                                masse2 = 0
                        if masse1 != 0:
                            mstreffer += 5
                        if masse2 != 0:
                            mstreffer += 5


            return mstreffer




            #zuerst muss überprüft werden ob es sich wirklich um ein Keton oder um ein Keten handelt und dann muss noch ausgeschlossen werden, dass es nicht in einer Zyklischen Verbindung ist

        def NMRspektrumAnalyse():

            def BewertungMethylgruppen():
                methylgruppennummer = self.elementgruppengrenzen[3]
                CH2gruppennummer = self.elementgruppengrenzen[2]
                CH1gruppennummer = self.elementgruppengrenzen[1]
                
                if self.elemente[methylgruppennummer] == 0: # Überprüfung ob es überhaupt Methylgruppen hat
                    return 0


                # 0 = Singlet / 1 = Doublet / 2 = Triplet
                methylkopplungen = [0,0,0]


                for atom in self.molekularstruktur:
                    if atom[0] == methylgruppennummer:
                        if atom[1][1] == CH2gruppennummer:
                            methylkopplungen[2] += 1
                        elif atom[1][1] == CH1gruppennummer:
                            methylkopplungen[1] += 1
                        else:
                            methylkopplungen[0] += 1

                nmrmethylgruppen = []

                cvletzerpeakfrequenz = -50
                for wert in Molekuelinfo.d20nmrdaten:
                    if wert[1] < 1.2: #Ausschluss, damit es kein CH2 usw. sein kann. Ist nicht wirklich eine fixer Wert, aber eine gute Approximation
                        if abs(wert[0] - cvletzerpeakfrequenz) >= 25:
                            nmrmethylgruppen.append(1)
                        else:
                            nmrmethylgruppen[-1] += 1


                fehlerkorrektur = 0 # Um das Fehlinterpretieren vorzubeigen bei überschneidungen, wird dieser Wert erhöht, wenn eine Wert vorkommt, welcher nicht alleine stehen kann und daher eine Kombination vorliegen muss
                for nmrgruppe in nmrmethylgruppen:
                    if nmrgruppe <= 3:
                        methylkopplungen[nmrgruppe-1] -= 1
                    elif nmrgruppe == 4 or nmrgruppe == 5:
                        fehlerkorrektur += 2
                    else:
                        fehlerkorrektur += 3

                abweichung = 0
                for wert in methylkopplungen:
                    abweichung += abs(wert)
                abweichung -= fehlerkorrektur

                if abweichung <= 0:
                    return 0
                else:
                    return abweichung

            def BewertungvonCH2undCH1gruppenn():
                approximation = [None for a in range(len(self.molekularstruktur))]

                for position, atom in enumerate(self.molekularstruktur):
                    if atom[0] == self.elementgruppengrenzen[3]:  # ist ein CH3
                        approximation[position] = [1]
                        # spezifikationen
                        if atom[1][1] == self.elementgruppengrenzen[2] + 2:
                            approximation[position][0] += 2

                    elif atom[0] == self.elementgruppengrenzen[2]:  # ist ein CH2
                        if len(atom) == 2:  # es ist ein Olephin
                            approximation[position] = [5.5]
                        else:
                            approximation[position] = [2]


                    elif atom[0] == self.elementgruppengrenzen[1]:  # ist ein CH
                        if self.Iszyklisch(position):
                            if self.Cyclogrösse(
                                    position) == 0:  # Spezialfall wenn es ein Aromat ist, ansonsten hat es ca den wert 2.3
                                approximation[position] = [7.3]

                            else:
                                approximation[position] = [2.3]

                        else:
                            approximation[position] = [2.3]


                    elif atom[0] == self.elementgruppengrenzen[3] + 2:  # ist ein C=O(H)
                        approximation[position] = [9.5]

                for a in range(len(approximation)):
                    if approximation[a] != None:
                        kopplungen = self.PositionendernachbarnwelcheeineWasserstoffkopplungeingehen(a)
                        if kopplungen != None:
                            for kopplung in kopplungen:
                                approximation[a].append(kopplung)

                # Nun hat es in der Liste noch gewisse Stellen mit None. Diese Müssen jetzt gelöscht werden und alle Positionen müssen angepasst werden

                # zuerst müssen noch lehre [] gelöscht werden, welche wege self.PositionendernachbarnwelcheeineWasserstoffkopplungeingehen

                for a in range(len(approximation) - 1, -1, -1):
                    if approximation[a] == None:
                        for b in range(len(approximation)):
                            if approximation[b] != None:
                                for c in range(1, len(approximation[b])):
                                    if approximation[b][c][0] > a:
                                        approximation[b][c][0] -= 1
                        approximation.pop(a)
                """
                print("approximation")
                print(approximation)
                print(self.molekularstruktur)
                """


                # nun werden die NMR daten Analysiert
                nmrwerte = copy.deepcopy(Molekuelinfo.gruppierted20nmrdaten)

                # jeder approximationswert wird mit jedem nmrwert ausprobiert und verglichen
                def nmrwertemitapproximationvergleich(nmrwert, aproximierterwert):
                    unterschied = (nmrwert[0] - aproximierterwert[0]) ** 2
                    unterschied += abs(len(nmrwert) - len(aproximierterwert)) * 2
                    # Idee ob es ohne besser funktionert
                    """
                    if len(nmrwert) > len(aproximierterwert):
                        unterschied -= len(aproximierterwert) * 4
                    else:
                        unterschied -= len(nmrwert) * 4
                    """
                    return unterschied

                vergleiche = []

                for a in range(len(approximation)):
                    vergleiche.append([])
                    for b in range(len(nmrwerte)):
                        vergleiche[a].append(nmrwertemitapproximationvergleich(nmrwerte[b], approximation[a]))

                transformation = [None for a in
                                  range(len(approximation))]  # jeder approximationswert bekommt ein NMRwert zugewisen

                # Es darf auch mehrfachbelegungen geben (aber in einer genau bestimmen Anzahl)
                anzahlmehrfachbelegungen = len(approximation) - len(nmrwerte)

                for temp in range(len(approximation)):
                    # da der Wert möglichst klein sein soll , wird jetzt das Atom ausgewält, bei welchem das beste, das schlechteste ist
                    besterwertjederspalte = [None for a in range(len(approximation))]
                    for a in range(len(approximation)):
                        spalenwerte = []
                        for b in range(len(nmrwerte)):
                            if vergleiche[a][b] != None:
                                spalenwerte.append(vergleiche[a][b])
                        if spalenwerte != []:
                            besterwertjederspalte[a] = min(spalenwerte)
                    """
                    print("besterwertjederspalte")
                    print(besterwertjederspalte)
                    print("vergleiche")
                    print(vergleiche)
                    print("Transformation")
                    print(transformation)
                    print("Übrige mehrfachbelegungen")
                    print(anzahlmehrfachbelegungen)
                    """
                    schlechtesterWert = max([wert for wert in besterwertjederspalte if wert is not None])

                    breakbool = False
                    for a, wert in enumerate(besterwertjederspalte):

                        if wert == schlechtesterWert:
                            for b in range(len(nmrwerte)):
                                if vergleiche[a][b] == wert:
                                    if b in transformation:
                                        anzahlmehrfachbelegungen -= 1

                                    transformation[a] = b
                                    for c in range(len(nmrwerte)):
                                        vergleiche[a][c] = None

                                    if anzahlmehrfachbelegungen <= 0:
                                        for spalte in vergleiche:
                                            spalte[b] = None

                                    breakbool = True
                                    break
                        if breakbool:
                            break
                """
                print("nmrwerte")
                print(nmrwerte)
                print("approximationen")
                print(approximation)
                """
                # Noch nicht final aber jetzt kann mal alles getestet werden
                summe = 0
                for a in range(len(approximation)):
                    summe += nmrwertemitapproximationvergleich(nmrwerte[transformation[a]], approximation[a])

                return summe

            #print("Methylgruppen")
            """
            print("Methylgruppen")
            print(BewertungMethylgruppen())
            print("NMR ALLES")
            print(BewertungvonCH2undCH1gruppenn())
            """

            return  BewertungvonCH2undCH1gruppenn() # + BewertungMethylgruppen() wirkt mehr als Störfaktor, als es hilft














        if False:
            print("AbzugKetonAlkoholverbindungen")
            print(AbzugKetonAlkoholverbindungen())
            print("AbzugEthertransformationzuKeton_Aldehyd")
            print(AbzugEthertransformationzuKeton_Aldehyd())
            print("MSpeaküberprüfung")
            print(MSpeaküberprüfung())
            print("NMRspektrumAnalyse")
            print(NMRspektrumAnalyse())
        self.heuristikwert = BestrafungKetonAlkoholverbindungen()
        self.heuristikwert += BestrafungEthertransformationzuKeton_Aldehyd()
        self.heuristikwert += 3*MSpeaküberprüfung()
        self.heuristikwert += 1.5*NMRspektrumAnalyse()
        #print(heuristikwert)

        return #heuristikwert



    def Muation(self):
        #Alle Mutationsfunktionen returnen True wenn eine Mutation durchgeführt werden konnte und False, wenn keine Mutation durchgeführt wurde. Die Manipulation wird direkt am self.molekularstruktur abgeändert


        def MutationzweiergleichwertigerAtome(spezifischeposition = None):

            if spezifischeposition == None:
                positionatom1 = random.randint(self.elementgruppengrenzen[1], len(self.molekularstruktur) - 1)
            else:
                positionatom1 = spezifischeposition

            wertigkeitatom1 = self.Anzahlverbindungen(self.molekularstruktur[positionatom1][0])
            gleichwertigeatome = []
            for count, atom in enumerate(self.molekularstruktur):
                if atom[0] != self.molekularstruktur[positionatom1][0] and self.Anzahlverbindungen(
                        atom[0]) == wertigkeitatom1:
                    gleichwertigeatome.append(count)

            if gleichwertigeatome != []:
                positionatom2 = gleichwertigeatome[random.randint(0, len(gleichwertigeatome) - 1)]
                atom1 = self.molekularstruktur[positionatom1][0]
                atom2 = self.molekularstruktur[positionatom2][0]
                self.molekularstruktur[positionatom1][0] = atom2
                self.molekularstruktur[positionatom2][0] = atom1
                return True

            return False

        def MutationreduzierungMehrfachbindungzuanderemOrt():

            #überprüfen ob es überhaupt Mehrfachbindungen hat
            listeallerMehrfachbindungen = []
            for atomposition,atom in enumerate(self.molekularstruktur):
                for a in range(1,len(atom)):
                    if atom[a][0] > 1:
                        listeallerMehrfachbindungen.append([atomposition, atom[a][1]])

            if listeallerMehrfachbindungen == []:
                return False

            #print("liste aller möglichkeiten")
            #print(listeallerMehrfachbindungen)
            #Jetzt müssen alle einfach oder Doppelbindungen gesucht werden, an welchem an beiden Atomen davon zwei Einzelbindungen sind, welche mann dann "herübernehmen kann"
            substitutionsmöglichkeiten = [] # Immer Viereinträge einfachbindung mit Atom1, Atom1, Atom2 (Atom1 und Atom2 sind der Ort bei dem die Verbindung erhöht wird), einfachbindung mit Atom2


            def isteineweitereeinfachbindungbeibeidenvorhanden(positionatom1, positionatom2):
                boolatom1 = False
                boolatom2 = False
                positionatomineinfachverbindungmitatom1 = 0
                positionatomineinfachverbindungmitatom2 = 0
                for a in range(1,len(self.molekularstruktur[positionatom1])):
                    if self.molekularstruktur[positionatom1][a][0] == 1 and self.molekularstruktur[positionatom1][a][1] != positionatom2:
                        boolatom1 = True
                        positionatomineinfachverbindungmitatom1 = self.molekularstruktur[positionatom1][a][1]

                for a in range(1,len(self.molekularstruktur[positionatom2])):
                    if self.molekularstruktur[positionatom2][a][0] == 1 and self.molekularstruktur[positionatom2][a][1] != positionatom1:
                        boolatom2 = True
                        positionatomineinfachverbindungmitatom2 = self.molekularstruktur[positionatom2][a][1]

                if boolatom1 and boolatom2:
                        return True, positionatomineinfachverbindungmitatom1, positionatomineinfachverbindungmitatom2
                else:
                    return False,None,None

            for atomposition,atom in enumerate(self.molekularstruktur):
                if self.Anzahlverbindungen(atom[0]) > 1:
                    for a in range(1,len(atom)):
                        if atom[a][0] == 1:
                            antwort = isteineweitereeinfachbindungbeibeidenvorhanden(atomposition,atom[a][1])
                            if antwort[0] and (self.Anzahlverbindungen(atom[0]) >= 3 or self.Anzahlverbindungen(self.molekularstruktur[atom[a][1]][0]) >= 3) :
                                substitutionsmöglichkeiten.append([antwort[1],atomposition,atom[a][1],antwort[2]])
                        elif atom[a][0] == 2:
                            antwort = isteineweitereeinfachbindungbeibeidenvorhanden(atomposition,atom[a][1])
                            if antwort[0] and (self.Anzahlverbindungen(atom[0]) >= 4 or self.Anzahlverbindungen(self.molekularstruktur[atom[a][1]][0]) >= 4) :
                                substitutionsmöglichkeiten.append([antwort[1],atomposition,atom[a][1],antwort[2]])


            # jetzt wird sich für eine Mehrfachverbindung entschieden
            randommehrfachbindung = random.randint(0, len(listeallerMehrfachbindungen) - 1)
            reduzierdesAtom1 = listeallerMehrfachbindungen[randommehrfachbindung][0]
            reduzierdesAtom2 = listeallerMehrfachbindungen[randommehrfachbindung][1]

            #print("Substitutionsmöglichkeiten")
            #print(substitutionsmöglichkeiten)

            # die Substitutionsmöglichkeiten werden jetzt noch gefiltert
            for a in range(len(substitutionsmöglichkeiten)-1,-1,-1):
                if substitutionsmöglichkeiten[a][0] == reduzierdesAtom1 or substitutionsmöglichkeiten[a][0] == reduzierdesAtom2:
                    substitutionsmöglichkeiten.pop(a)

                elif substitutionsmöglichkeiten[a][3] == reduzierdesAtom1 or substitutionsmöglichkeiten[a][3] == reduzierdesAtom2:
                    substitutionsmöglichkeiten.pop(a)

                elif substitutionsmöglichkeiten[a][1] == reduzierdesAtom1 or substitutionsmöglichkeiten[a][1] == reduzierdesAtom2:
                    substitutionsmöglichkeiten.pop(a)

                elif substitutionsmöglichkeiten[a][2] == reduzierdesAtom1 or substitutionsmöglichkeiten[a][2] == reduzierdesAtom2:
                    substitutionsmöglichkeiten.pop(a)


            if substitutionsmöglichkeiten == []:
                return False
            #print(substitutionsmöglichkeiten)



            randomsubstitution = random.randint(0,len(substitutionsmöglichkeiten)-1)



            einfachbindungzusubstitutionsAtom1 = substitutionsmöglichkeiten[randomsubstitution][0]
            einfachbindungzusubstitutionsAtom2 = substitutionsmöglichkeiten[randomsubstitution][3]
            substitutionsAtom1 = substitutionsmöglichkeiten[randomsubstitution][1]
            substitutionsAtom2 = substitutionsmöglichkeiten[randomsubstitution][2]

            self.Creatverbindung(reduzierdesAtom1,reduzierdesAtom2,-1)

            self.Deletverbindung(einfachbindungzusubstitutionsAtom1,substitutionsAtom1)
            self.Creatverbindung(einfachbindungzusubstitutionsAtom1, reduzierdesAtom1, 1)
            self.Deletverbindung(einfachbindungzusubstitutionsAtom2, substitutionsAtom2)
            self.Creatverbindung(einfachbindungzusubstitutionsAtom2,reduzierdesAtom2,1)
            self.Creatverbindung(substitutionsAtom1, substitutionsAtom2, 1)

            return True



        def MutationreduzierungMehrfachbindungzueinemCyclo():
            return False

        def Muationentfernungcyclofüreinemehrfachbindung():
            return False

        def MutationAtomherausreisenundanderswoneueinsetzen():

            positionatom = random.randint(0, len(self.molekularstruktur) - 1)
            wertigkeitatom = self.Anzahlverbindungen(self.molekularstruktur[positionatom][0])

            if wertigkeitatom == 1:
                if MutationzweiergleichwertigerAtome(positionatom):
                    return True
            if wertigkeitatom == 2:
                if len(self.molekularstruktur[positionatom]) == 3:  # bzw es zwei einfachbindungen hat
                    möglichesubstitutionen = []
                    # Suche eine einfachbindung und schmuggle dich rein
                    for position, atom in enumerate(self.molekularstruktur):
                        for a in range(1, len(atom)):
                            if atom[a][0] == 1 and position != positionatom and atom[a][1] != positionatom:
                                möglichesubstitutionen.append([position, atom[a][1]])
                    if möglichesubstitutionen == []:
                        return False

                    substitutionsverbindung = random.randint(0, len(möglichesubstitutionen) - 1)

                    # das zwei Wertigkeatom muss zuerst von seinen zwei verbindungen herausgelöst werden, und die beidenen offnene stellen wieder zusammengeflickt werden
                    positionatom2 = self.molekularstruktur[positionatom][1][1]
                    positionatom3 = self.molekularstruktur[positionatom][2][1]
                    self.Deletverbindung(positionatom, positionatom2)
                    self.Deletverbindung(positionatom, positionatom3)
                    self.Creatverbindung(positionatom2, positionatom3, 1)

                    # das herausgelöste Atom wird nun wieder an bestimmten ort wieder eingesetzt
                    positionsubstitiontsatom1 = möglichesubstitutionen[substitutionsverbindung][0]
                    positionsubstitiontsatom2 = möglichesubstitutionen[substitutionsverbindung][1]
                    self.Deletverbindung(positionsubstitiontsatom1, positionsubstitiontsatom2)

                    # neue verbindungen werden erstellt
                    self.Creatverbindung(positionsubstitiontsatom1, positionatom, 1)
                    self.Creatverbindung(positionsubstitiontsatom2, positionatom, 1)



                    return True
                else:
                    return False

            if wertigkeitatom == 3:
                if len(self.molekularstruktur[positionatom]) == 4:  # bzw es drei einfachbindungen hat. Deshalb werden dann zwei miteinander verknüpft, der dritte bleibt und das Atom kann ziwschen eine einfach Bindung sich einfügen:
                    möglichesubstitutionen = []

                    # Suche eine einfachbindung und schmuggle dich rein
                    for position, atom in enumerate(self.molekularstruktur):
                        for a in range(1, len(atom)):
                            if atom[a][0] == 1 and position != positionatom and atom[a][1] != positionatom: #die beiden letzeren schliessen aus, dass die dritte Verbindung die ist, welche gespalten wird, anonsten kann es passieren, dass das ausgewählte Atom eine Bindung mit sich selbst eingeht, was unbedingt vermieden werden muss
                                möglichesubstitutionen.append([position, atom[a][1]])
                    if möglichesubstitutionen == []:
                        return False

                    substitutionsverbindung = random.randint(0, len(möglichesubstitutionen) - 1)

                    # zwei einfachbindungen werden aufgespalten. Es spielt keine Rolle und da die Reihenfolge sowiso zufällig ist, nehme ich einfach die ersten zwei
                    positionatom2 = self.molekularstruktur[positionatom][1][1]
                    positionatom3 = self.molekularstruktur[positionatom][2][1]
                    positionverbleibendesatom = self.molekularstruktur[positionatom][3][1]

                    if self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstruktur[positionatom2], positionatom3) >= 3: #muss ausgeschlossen werden, da es sonst eine Vierfachbindung geben könnte, welche nicht existieren darf
                        return False

                    if self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstruktur[positionatom2], positionatom3)*2 == self.Anzahlverbindungen(self.molekularstruktur[positionatom2][0]) + self.Anzahlverbindungen(self.molekularstruktur[positionatom3][0]) - 2 : #Zudem muss ausgeschlossen werden, dass es sich um ein Dreieck handelt und durch diese Manipulation zwei Teile vom Rest des Moleküls sich loslösen
                        return False




                    self.Deletverbindung(positionatom, positionatom2)
                    self.Deletverbindung(positionatom, positionatom3)

                    #Rein Topologisch muss es jetzt immer noch alles zusammenhängend sein, wäre dies nicht mehr der Fall muss der Mutationsversuch abgebrochen werden, da sich daraus zwei unterschiedliche Moleküle bilden würden

                    if not self.isligit(True,True,False):# es muss nur die Topologie geprüft werden
                        # Fall ist aufgetreten und nun wird alles wieder rückgäng gemacht werden
                        self.Creatverbindung(positionatom, positionatom2,1)
                        self.Creatverbindung(positionatom, positionatom3,1)
                        return False


                    self.Creatverbindung(positionatom2, positionatom3, 1)



                    # das herausgelöste Atom wird nun wieder an bestimmten ort wieder eingesetzt
                    positionsubstitiontsatom1 = möglichesubstitutionen[substitutionsverbindung][0]
                    positionsubstitiontsatom2 = möglichesubstitutionen[substitutionsverbindung][1]
                    self.Deletverbindung(positionsubstitiontsatom1, positionsubstitiontsatom2)

                    # neue verbindungen werden erstellt
                    self.Creatverbindung(positionsubstitiontsatom1, positionatom, 1)
                    self.Creatverbindung(positionsubstitiontsatom2, positionatom,1)



                    return True


                elif len(self.molekularstruktur[positionatom]) == 3:  # Bedeutet, dass das Atom eine Doppelbindung und eine einfachbindung hat
                    möglichesubstitutionen = []
                    # Suche eine einfachbindung und schmuggle dich rein
                    for position, atom in enumerate(self.molekularstruktur):
                        for a in range(1, len(atom)):
                            if atom[a][0] == 1 and position != positionatom:
                                möglichesubstitutionen.append([position, a])
                    if möglichesubstitutionen == []:
                        return False

                    substitutionsverbindung = random.randint(0, len(möglichesubstitutionen) - 1)

                    # die Einfachbindung wird aufgespalten und die Doppelbindung wird zu einer Einfachbindung. Es ist nicht optimal, um das Molekül perfekt wegzunehmen, aber ansonsten hat man probleme "verschliese" der offenene Verbindungen
                    if self.molekularstruktur[positionatom][1][0] == 1:
                        positionatomeinfachbindung = self.molekularstruktur[positionatom][1][1]
                        positionatomdoppelbindung = self.molekularstruktur[positionatom][2][1]
                    else:
                        positionatomeinfachbindung = self.molekularstruktur[positionatom][2][1]
                        positionatomdoppelbindung = self.molekularstruktur[positionatom][1][1]

                    positionsubstitiontsatom1 = möglichesubstitutionen[substitutionsverbindung][0]
                    positionsubstitiontsatom2 = self.molekularstruktur[positionsubstitiontsatom1][möglichesubstitutionen[substitutionsverbindung][1]][1]

                    self.Deletverbindung(positionatom, positionatomeinfachbindung)
                    self.Deletverbindung(positionatom, positionatomdoppelbindung)
                    self.Creatverbindung(positionatomeinfachbindung, positionatomdoppelbindung, 1)
                    self.Creatverbindung(positionatomdoppelbindung, positionatom, 1)

                    positionverbleibendesatom = self.molekularstruktur[positionatom][1][
                        1]  # es kann sein, dass mit diesem Atom eine neue Bindung eingeganen werden soll, dann muss es aber anstatt einer neuer einfachbindung eine doppelbindung geben

                    # das herausgelöste Atom wird nun wieder an bestimmten ort wieder eingesetzt
                    self.Deletverbindung(positionsubstitiontsatom1, positionsubstitiontsatom2)

                    # neue verbindungen werden erstellt
                    if positionverbleibendesatom != positionsubstitiontsatom1:
                        self.Creatverbindung(positionsubstitiontsatom1, positionatom, 1)
                    else:
                        self.Deletverbindung(positionatom, positionverbleibendesatom)
                        self.Creatverbindung(positionatom, positionverbleibendesatom, 2)

                    if positionverbleibendesatom != positionsubstitiontsatom2:
                        self.Creatverbindung(positionsubstitiontsatom2, positionatom, 1)
                    else:
                        self.Deletverbindung(positionatom, positionverbleibendesatom)
                        self.Creatverbindung(positionatom, positionverbleibendesatom, 2)  #

                    return True






                # es ist mit abstand der schwierigste Fall und somit lasse ich ihn für den Moment aus, sollte keine riesigen auswirkungen haben


                #if wertigkeitatom == 4:
                    #return False










        randommutationszeiger = random.random()

        if randommutationszeiger < 0.4:
            #print("MutationzweiergleichwertigerAtome")
            for a in range(10): # Wird mehrmals gemacht, da es nicht immer funktionert. Wenn es funktioniert, wird abgebrochen
                zwischenspeichermolekularstruktur = copy.deepcopy(self.molekularstruktur)
                if MutationzweiergleichwertigerAtome():
                    if self.isligit():
                        if not self.testfunktiondersymetriebeidenBindungen(self.molekularstruktur):
                            raise Exception("Bindungsymetrie ist nicht gegeben")
                        return True
                    else:
                        self.molekularstruktur = copy.deepcopy(zwischenspeichermolekularstruktur)



        elif randommutationszeiger < 0.7:
            #print("MutationAtomherausreisenundanderswoneueinsetzen ")
            for a in range(10):
                zwischenspeichermolekularstruktur = copy.deepcopy(self.molekularstruktur)
                if MutationAtomherausreisenundanderswoneueinsetzen():
                    if self.isligit():
                        if not self.testfunktiondersymetriebeidenBindungen(self.molekularstruktur):
                            raise Exception("Bindungsymetrie ist nicht gegeben")
                        return True
                    else:
                        self.molekularstruktur = copy.deepcopy(zwischenspeichermolekularstruktur)
        else:
            #print("MutationverschiebungMehrfachbindung")
            for a in range(2):
                zwischenspeichermolekularstruktur = copy.deepcopy(self.molekularstruktur)
                if MutationreduzierungMehrfachbindungzuanderemOrt():
                    # Es gibt leider sehr viele Spezialfälle, bei welchem diese Mutation nicht geht. Diese sind aber meistens so speziel, dass es sich nicht lohnt, diese manuel einzelt herauszufiltern. Daher überprüfe ich das Molekül und wenn es nicht mehr zusammenhängend ist und falls dies der Fall ist, muss ich alles wieder rückgängig machen.
                    if self.isligit():
                        if not self.testfunktiondersymetriebeidenBindungen(self.molekularstruktur):
                            raise Exception("Bindungsymetrie ist nicht gegeben")
                        return True
                    else:
                        self.molekularstruktur = copy.deepcopy(zwischenspeichermolekularstruktur)



        return False # Mutation hat nicht funktioniert.

        
        
        
        

    def SMilestransformator(self):
        def Übersetzung(atomnummer, wertigkeit):
            # 0 = C / 1 = CH / 2 = CH2 / 3 = Keton / 4 = Ether / 5 = CH3 / 6 = OH / 7 = Aldehyd / 8 = Carbonsäure
            if wertigkeit == 1:
                if atomnummer == 0:
                    return "(C)"
                if atomnummer == 1:
                    return "(C)"
                if atomnummer == 2:
                    return "(C)"
                if atomnummer == 3:
                    return "(C(=O))"
                if atomnummer == 4:
                    return "(O)"
                if atomnummer == 5:
                    return "(C)"
                if atomnummer == 6:
                    return "(O)"
                if atomnummer == 7:
                    return "(C(=O))"
                if atomnummer == 8:
                    return "(C(=O)(O))"
            if wertigkeit == 2:
                if atomnummer == 0:
                    return "(=C)"
                if atomnummer == 1:
                    return "(=C)"
                if atomnummer == 2:
                    return "(=C)"
                if atomnummer == 3:
                    return "(=C(=O))"
                if atomnummer == 4:
                    return "(=O)"

            if wertigkeit == 3:
                if atomnummer == 0:
                    return "(#C)"
                if atomnummer == 1:
                    return "(#C)"

        def LängederÜbersetzung(atomnummer):
            if atomnummer == 0:
                return 3
            if atomnummer == 1:
                return 3
            if atomnummer == 2:
                return 3
            if atomnummer == 3:
                return 7
            if atomnummer == 4:
                return 3
            if atomnummer == 5:
                return 3
            if atomnummer == 6:
                return 3
            if atomnummer == 7:
                return 7
            if atomnummer == 8:
                return 10

        def PositionenimStringanpassen(liste, ortderänderung, anzahlzeichen):

            for a in range(len(liste)):
                if liste[a] != None:
                    if liste[a][0] >= ortderänderung:
                        liste[a][0] += anzahlzeichen
                    if liste[a][1] >= ortderänderung:
                        liste[a][1] += anzahlzeichen
            return liste


        gemachteAtomgruppen = [False for a in range(len(self.molekularstruktur))]
        ortdergemachtenAtomgruppen = [None for a in range(len(self.molekularstruktur))]
        verbindungennochüberprüfen= []
        cyclozähler = 0





        #ich fange immer mit dem ersten an

        gemachteAtomgruppen[0] = True
        ortdergemachtenAtomgruppen[0] = [0,2]
        verbindungennochüberprüfen.append(0)

        outputstring = Übersetzung(self.molekularstruktur[0][0],1)

        outputstring = "|" + outputstring[1:-1] + "|"


        molekularstrukturcopy = copy.deepcopy(self.molekularstruktur) # Da im prozess die molekularstruktur verändern bzw. gelöscht wird, wird sie hier in einer anderen Variabel gespeichert und am Schluss wieder zum originalzustand zurückgeführt

        while verbindungennochüberprüfen != []:



            momentanesAtom = verbindungennochüberprüfen[0]

            verbindendesAtom = self.molekularstruktur[momentanesAtom][1][1]

            wertigkeitderverbindung = self.molekularstruktur[momentanesAtom][1][0]
            self.Deletverbindung(momentanesAtom,verbindendesAtom)



            if gemachteAtomgruppen[verbindendesAtom]: # ergibt True, wenn es schon gemacht wurde

                cyclozähler += 1

                if cyclozähler < 10:
                    zusatzstring = str(cyclozähler)
                    längenaddition = 1
                else:
                    zusatzstring = "%" + str(cyclozähler)
                    längenaddition = len(zusatzstring)

                mehrfachbindungsstringaddition = ""
                längenadditiondurchmehrfachbindung = 0
                if wertigkeitderverbindung == 2:
                    mehrfachbindungsstringaddition = "="
                    längenadditiondurchmehrfachbindung = 1
                if wertigkeitderverbindung == 3:
                    mehrfachbindungsstringaddition = "#"
                    längenadditiondurchmehrfachbindung = 1


                ortimstringumanzusetzen1 = ortdergemachtenAtomgruppen[momentanesAtom][0] + 2
                ortimstringumanzusetzen2 = ortdergemachtenAtomgruppen[verbindendesAtom][0] + 2

                # Da es vorkommen kann, dass vor dem Atom C ein Zeichen für eine Doppelbindung oder eine Dreifachbindung ist, muss der Ort jetzt noch angepasst werden, dass dies nicht passieren kann
                #es kann aber maximal um eine stelle verschoben werden

                if outputstring[ortimstringumanzusetzen1-1] == "#" or outputstring[ortimstringumanzusetzen1-1] == "=":
                    ortimstringumanzusetzen1 += 1
                if outputstring[ortimstringumanzusetzen2-1] == "#" or outputstring[ortimstringumanzusetzen2-1] == "=":
                    ortimstringumanzusetzen2 += 1



                outputstring = outputstring[:ortimstringumanzusetzen1] + mehrfachbindungsstringaddition + zusatzstring + outputstring[ortimstringumanzusetzen1:]
                ortdergemachtenAtomgruppen = PositionenimStringanpassen(ortdergemachtenAtomgruppen, ortimstringumanzusetzen1, längenadditiondurchmehrfachbindung + längenaddition)

                if ortimstringumanzusetzen1 < ortimstringumanzusetzen2:
                    ortimstringumanzusetzen2 += längenadditiondurchmehrfachbindung + längenaddition

                outputstring = outputstring[:ortimstringumanzusetzen2] + zusatzstring + outputstring[ortimstringumanzusetzen2:]
                ortdergemachtenAtomgruppen = PositionenimStringanpassen(ortdergemachtenAtomgruppen,ortimstringumanzusetzen2,längenaddition)



            else:
                #wird vor der schliessenenden Klammer am ende Hinzugefügt und in der liste OrtdergemachtenAtomgruppen vermekrt, daraufhin müssen alle Werte der Liste angepasst werden
                gemachteAtomgruppen[verbindendesAtom] = True
                verbindungennochüberprüfen.append(verbindendesAtom)
                ortimstringumanzusetzen = ortdergemachtenAtomgruppen[momentanesAtom][1]  # Minus 1 wegen der Klammer die man überspringen muss

                längenadditiondurchmehrfachbindung = 0
                if wertigkeitderverbindung > 1:
                    längenadditiondurchmehrfachbindung = 1

                if wertigkeitderverbindung >= 4:
                    raise ValueError("Mehr als eine Dreifachbindung, dies ist nicht möglich")


                outputstring = outputstring[:ortimstringumanzusetzen] + Übersetzung(self.molekularstruktur[verbindendesAtom][0],wertigkeitderverbindung) + outputstring[ortimstringumanzusetzen:]

                ortdergemachtenAtomgruppen = PositionenimStringanpassen(ortdergemachtenAtomgruppen,ortimstringumanzusetzen,LängederÜbersetzung(self.molekularstruktur[verbindendesAtom][0])+längenadditiondurchmehrfachbindung)

                ortdergemachtenAtomgruppen[verbindendesAtom] = [ortimstringumanzusetzen,ortimstringumanzusetzen + LängederÜbersetzung(self.molekularstruktur[verbindendesAtom][0]) + längenadditiondurchmehrfachbindung -1]





            for a in range(len(verbindungennochüberprüfen)-1,-1,-1):
                if (len(self.molekularstruktur[verbindungennochüberprüfen[a]]) == 1):
                    verbindungennochüberprüfen.pop(a)


        #Die eckigen Klammern müsseen entfern werden, da es um das erste Atom keine Klammern haben darf. Damit der Algoryhmus aber immer perfekt klappt, muss man sie am Anfang hinzufügen und jetzt wieder rausnehmen

        outputstring = outputstring.replace("|","")


        self.smistring = outputstring    #Muss bei Smi einfach so sein, es darf am Anfang und am Ende keine Klammer haben. Da es die äusserten Klammern sind, spielt es auch keine Rolle
        print(self.smistring)
        # Da self.molekularstruktur jetzt eigentlich gelöscht ist, wird dies wieder neu mit der Copy vom Anfang überspielt
        self.molekularstruktur = copy.deepcopy(molekularstrukturcopy)

        return


    def DarstellungMolekülinSMI(self, mitWasserstoff = False, mitBild = True, printSMLIESstring = False):
        if self.smistring != None:

            mol = Chem.MolFromSmiles(self.smistring)
            if mitWasserstoff:
                mol = Chem.AddHs(mol)
            rdkit.Chem.AllChem.Compute2DCoords(mol)
            if mitBild:
                img = Draw.MolToImage(mol, size=(2000, 2000), kekulize=True, bgcolor=(255, 255, 255))
                img.show()
            self.smistring = Chem.MolToSmiles(mol)
            if printSMLIESstring:
                print(self.smistring)
        return


    def Anzahlverbindungen(self, elementnummer):
        if elementnummer >= self.elementgruppengrenzen[3]:
            return 1
        if elementnummer >= self.elementgruppengrenzen[2]:
            return 2
        if elementnummer >= self.elementgruppengrenzen[1]:
            return 3
        if elementnummer >= self.elementgruppengrenzen[0]:
            return 4

    def Strukturintertialgenerator(self,anzahlmutationen):
        molekularsturktur = [[0] * 1 for _ in range(sum(self.elemente))]




        count = 0
        for a in range(len(self.elemente)):
            for b in range(self.elemente[a]):
                molekularsturktur[count][0] = a
                count += 1


        #jetzt kann das wirkliche zusammenbauen eines provisorischen Moleküles beginnen

        offeneverbindungen = self.Anzahloffeneverbindungen(molekularsturktur[0])
        Doppelbindungsequivalenz = (self.elemente[0]*2 + self.elemente[1] -self.elemente[5] - self.elemente[6] -self.elemente[7] -self.elemente[8])/2 + 1


        for a in range(1,len(molekularsturktur),1):
            maxoffeneverbindungen = 0
            ortmax = -0
            for b in range(a):
                if self.Anzahloffeneverbindungen(molekularsturktur[b]) > maxoffeneverbindungen:
                    maxoffeneverbindungen = self.Anzahloffeneverbindungen(molekularsturktur[b])
                    ortmax = b



            if self.Anzahlverbindungen(molekularsturktur[a][0]) == 4:
                if maxoffeneverbindungen >= 3 and Doppelbindungsequivalenz >= 2 and offeneverbindungen >= 3:
                    molekularsturktur[a].append([3,ortmax])
                    molekularsturktur[ortmax].append([3,a])
                    offeneverbindungen -= 2
                    Doppelbindungsequivalenz -= 2

                elif maxoffeneverbindungen >= 2 and Doppelbindungsequivalenz >= 1:
                    molekularsturktur[a].append([2, ortmax])
                    molekularsturktur[ortmax].append([2, a])
                    Doppelbindungsequivalenz -= 1
                    #Anzahl der offenenverbindungen bleibt gleich
                elif maxoffeneverbindungen >= 1:
                    molekularsturktur[a].append([1, ortmax])
                    molekularsturktur[ortmax].append([1, a])
                    offeneverbindungen += 2

            if self.Anzahlverbindungen(molekularsturktur[a][0]) == 3:
                if maxoffeneverbindungen >= 3 and Doppelbindungsequivalenz >= 2 and offeneverbindungen >= 4:
                    molekularsturktur[a].append([3, ortmax])
                    molekularsturktur[ortmax].append([3, a])
                    offeneverbindungen -= 3
                    Doppelbindungsequivalenz -= 2

                elif maxoffeneverbindungen >= 2 and Doppelbindungsequivalenz >= 1 and offeneverbindungen >= 2:
                    molekularsturktur[a].append([2, ortmax])
                    molekularsturktur[ortmax].append([2, a])
                    Doppelbindungsequivalenz -= 1
                    offeneverbindungen -= 1

                elif maxoffeneverbindungen >= 1:
                    molekularsturktur[a].append([1, ortmax])
                    molekularsturktur[ortmax].append([1, a])
                    offeneverbindungen += 1

            if self.Anzahlverbindungen(molekularsturktur[a][0]) == 2:

                if maxoffeneverbindungen >= 2 and Doppelbindungsequivalenz >= 1 and offeneverbindungen >= 3:
                    molekularsturktur[a].append([2, ortmax])
                    molekularsturktur[ortmax].append([2, a])
                    Doppelbindungsequivalenz -= 1
                    offeneverbindungen -= 2

                elif maxoffeneverbindungen >= 1:
                    molekularsturktur[a].append([1, ortmax])
                    molekularsturktur[ortmax].append([1, a])

            if self.Anzahlverbindungen(molekularsturktur[a][0]) == 1:
                molekularsturktur[a].append([1, ortmax])
                molekularsturktur[ortmax].append([1, a])
                offeneverbindungen -= 1



        #Es kann sein, dass es noch Atome hat, welche keine 2 Verbindungen oder mehr eingehen können, aber noch nicht genug verbindungen haben, dafür muss das Doppeblindungsequivalenz aber grösser als 0 sein. Danach versuche das Program die "Löcher" zu schliessen

        if Doppelbindungsequivalenz > 0:


            restlicheoffenestellen = []
            for atom in molekularsturktur:
                restlicheoffenestellen.append(self.Anzahloffeneverbindungen(atom))

            
            # zuerst dreifachoffene stellen werden gesucht
            positionen = []

            for count, wert in enumerate(restlicheoffenestellen):
                if wert == 3:
                    positionen.append(count)

            if positionen != []:

                a = 0
                while a < len(positionen):
                    for b in range(a+1,len(positionen)):
                        if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],positionen[b]) == 0:

                            molekularsturktur[positionen[a]].append([3, positionen[b]])
                            molekularsturktur[positionen[b]].append([3, positionen[a]])
                            Doppelbindungsequivalenz -= 3
                            positionen.pop(b)
                            break
                    if self.Anzahloffeneverbindungen(molekularsturktur[positionen[a]]) == 0:
                        positionen.pop(a)
                    else:
                        a += 1

            # zweifachoffene Stellen versuchen zu schliessen
            for count, wert in enumerate(restlicheoffenestellen):
                if wert == 2:
                    positionen.append(count)


            if positionen != []:
                a = 0
                while a < len(positionen):
                    for b in range(a+1, len(positionen)):
                        if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],positionen[b]) == 0:
                            molekularsturktur[positionen[a]].append([2, positionen[b]])
                            molekularsturktur[positionen[b]].append([2, positionen[a]])
                            Doppelbindungsequivalenz -= 2
                            positionen.pop(b)
                            break
                    if self.Anzahloffeneverbindungen(molekularsturktur[positionen[a]]) <= 1:
                        positionen.pop(a)
                    else:
                        a += 1


            # die einfachen STellen werden am schluss auch noch versucht zu schliessen
            positionen = []
            for count, wert in enumerate(restlicheoffenestellen):
                if wert >= 1:
                    positionen.append(count)

            if positionen != []:

                a = 0
                while a < len(positionen):
                    for b in range(a+1, len(positionen)):

                        if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],positionen[b]) == 0:
                            molekularsturktur[positionen[a]].append([1, positionen[b]])
                            molekularsturktur[positionen[b]].append([1, positionen[a]])
                            Doppelbindungsequivalenz -= 1
                            if self.Anzahloffeneverbindungen(molekularsturktur[positionen[b]]) == 0:
                                positionen.pop(b)
                            break
                        elif self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],positionen[b]) <= 2:

                            for c in range(1,len(molekularsturktur[positionen[a]])):
                                if molekularsturktur[positionen[a]][c][1] == positionen[b]:
                                    molekularsturktur[positionen[a]][c][0] += 1

                            for c in range(1,len(molekularsturktur[positionen[b]])):
                                if molekularsturktur[positionen[b]][c][1] == positionen[a]:
                                    molekularsturktur[positionen[b]][c][0] += 1

                            if self.Anzahloffeneverbindungen(molekularsturktur[positionen[b]]) == 0:
                                positionen.pop(b)
                            break

                    if self.Anzahloffeneverbindungen(molekularsturktur[positionen[a]]) == 0:
                        positionen.pop(a)
                    else:
                        a += 1
            
            # es kann leider am Schluss noch sein, dass eine Molekül noch zwei freie Verbindungen hat, diese werden jetzt noch mit einer Art der Elektophilen Additon geschlossen
            for count, atom in enumerate(molekularsturktur):
                if (self.Anzahloffeneverbindungen(atom)) == 2:
                    #Nun wird versucht eine Doppelbindung oder dreifachbindung um eins zu reduzieren und dann mit den neuen zwei offnene stellen das Atom damit zu verbinden
                    breakbool = False
                    for a in range(len(molekularsturktur)):
                        for b in range(1,len(molekularsturktur[a])):
                            if molekularsturktur[a][b][0] >= 2 and a != count:
                                breakbool = True
                                molekularsturktur[a][b][0] -= 1
                                for c in range(1,len(molekularsturktur[molekularsturktur[a][b][1]])):
                                    if molekularsturktur[molekularsturktur[a][b][1]][c][1] == a:
                                        molekularsturktur[molekularsturktur[a][b][1]][c][0] -= 1
                                        break
                                Doppelbindungsequivalenz -= 1


                                if self.AnzahlverbindungenzwischenzweiAtomen(atom,a) == 0:
                                    molekularsturktur[count].append([1,a])
                                    molekularsturktur[a].append([1,count])
                                else:
                                    for c in range(1,len(molekularsturktur[a])):
                                        if molekularsturktur[a][c][1] == count:
                                            molekularsturktur[a][c][0] += 1
                                            break
                                if self.AnzahlverbindungenzwischenzweiAtomen(atom,molekularsturktur[a][b][1]) == 0:
                                    molekularsturktur[count].append([1,molekularsturktur[a][b][1]])
                                    molekularsturktur[molekularsturktur[a][b][1]].append([1,count])

                                else:
                                    for c in range(1,len(molekularsturktur[molekularsturktur[a][b][1]])):
                                        if molekularsturktur[a][c][1] == count:
                                            molekularsturktur[a][c][0] += 1
                                            break

                        if breakbool:
                            break





        if Doppelbindungsequivalenz // 1 == Doppelbindungsequivalenz and Doppelbindungsequivalenz >= 0:
            for atom in molekularsturktur:
                if self.Anzahloffeneverbindungen(atom) != 0:
                    print(self.gruppenkonfiguration)

        if False:
            for atom in molekularsturktur:
                print(Anzahloffeneverbindungen(atom))
            print("Doppelbindungsequivalenz:" + str(Doppelbindungsequivalenz))
            print(molekularsturktur)

        self.molekularstruktur = molekularsturktur
        for a in range(anzahlmutationen):
            self.Muation()


        return


    def ChildauszweiParents(self, molekularstrukturMutter, molekularstrukturVater):

        def Anzahlverbundeneelemente(molekularstruktur):
            # Die Funktion gibt wieder wie viele Elemente mit dem Element an der Position 0 zusammenhängen

            schondurchlofeneelementepositionen = [0]
            nochzudurchlaufendeelementpositionen = [0]
            anzahlelemente = 0

            while nochzudurchlaufendeelementpositionen != []:
                atom = molekularstruktur[nochzudurchlaufendeelementpositionen[0]]
                nochzudurchlaufendeelementpositionen.pop(0)
                anzahlelemente += 1
                for a in range(1, len(atom)):
                    if not atom[a][1] in schondurchlofeneelementepositionen:
                        schondurchlofeneelementepositionen.append(atom[a][1])
                        nochzudurchlaufendeelementpositionen.append(atom[a][1])

            return anzahlelemente

        def Molekularstrukturzerschneider(molekularstruktur):
            #die Molekülstruktur wird bis zu 3 mal zerschnitten und sobald es zwei einzelne Stücke gibt es return Molekularstruktur
            #Es wird genau maximal 3 mal derschnitten, dass egal wie die Molekularstruktur aussieht es zerschnitten werden kann (Problem wegen Cycloverbindungen)


            testcopy = copy.deepcopy(molekularstruktur)
            anzahlelemente = sum(self.elemente)

            for a in range(3):
                ausgewähltesatom = random.randint(0, anzahlelemente-1)
                if len(molekularstruktur[ausgewähltesatom]) != 1:
                    zweitesatom = molekularstruktur[ausgewähltesatom][random.randint(1,len(molekularstruktur[ausgewähltesatom])-1)][1]

                    for a in range(1,len(molekularstruktur[ausgewähltesatom])):
                        if molekularstruktur[ausgewähltesatom][a][1] == zweitesatom:
                            molekularstruktur[ausgewähltesatom].pop(a)
                            break

                    for a in range(1,len(molekularstruktur[zweitesatom])):
                        if molekularstruktur[zweitesatom][a][1] == ausgewähltesatom:
                            molekularstruktur[zweitesatom].pop(a)
                            break

                    if anzahlelemente != Anzahlverbundeneelemente(molekularstruktur):
                        #self.testfunktion(testcopy,molekularstruktur)
                        return molekularstruktur,Anzahlverbundeneelemente(molekularstruktur)
            return False

        def Molekularstrukturzusammenfügen(zerschnittenemolekularstrukturMutter, zerschnittenemolekularstrukturVater):
            molekularstruktur1 = [] # kommt alles was mit [0] bei der Molekularstruktur des Vaters verbunden ist und der rest der Mutter
            molekularstruktur1übersetzungsaltepositionen = ["None" for a in range(sum(self.elemente))]
            molekularstruktur2 = [] # kommt alles WAS mit [0] bei der Molekularstruktur der Mutter verbunden ist und der Rest des Vater
            molekularstruktur2übersetzungsaltepositionen = ["None" for a in range(sum(self.elemente))]



            schondurchlofeneelementepositionen = [0]
            nochzudurchlaufendeelementpositionen = [0]




            while nochzudurchlaufendeelementpositionen != []:
                atom = zerschnittenemolekularstrukturVater[nochzudurchlaufendeelementpositionen[0]]
                molekularstruktur1.append(atom)
                molekularstruktur1übersetzungsaltepositionen[nochzudurchlaufendeelementpositionen[0]] = len(molekularstruktur1)-1
                nochzudurchlaufendeelementpositionen.pop(0)
                for a in range(1, len(atom)):
                    if not atom[a][1] in schondurchlofeneelementepositionen:
                        schondurchlofeneelementepositionen.append(atom[a][1])
                        nochzudurchlaufendeelementpositionen.append(atom[a][1])



            #bereinigung der positionen
            # Da bei einer Verbindung immer angegeben wird wie stark sie ist und dann die Position des anderen Elements, stimmmt jetzt alles nicht mehr da die Positionen nun anderst sind. Daher ist in den übersetzungsaltpositionen die alten Positionen gespeichert, damit dies nun korrigiert werden kann.

            for atom in molekularstruktur1:
                for a in range(1,len(atom)):
                    atom[a][1] = molekularstruktur1übersetzungsaltepositionen[atom[a][1]]

            #die restlichen Elemente gehen zur molekularstruktur2
            for position,atom in enumerate(zerschnittenemolekularstrukturVater):
                if not position in schondurchlofeneelementepositionen:
                    molekularstruktur2.append(atom)
                    molekularstruktur2übersetzungsaltepositionen[position] = len(molekularstruktur2)-1


            for atom in molekularstruktur2:
                for a in range(1,len(atom)):
                    atom[a][1] = molekularstruktur2übersetzungsaltepositionen[atom[a][1]]

            längemolekularstruktur1 = len(molekularstruktur1)
            längemolekularstruktur2 = len(molekularstruktur2)




            schondurchlofeneelementepositionen = [0]
            nochzudurchlaufendeelementpositionen = [0]

            while nochzudurchlaufendeelementpositionen != []:
                atom = zerschnittenemolekularstrukturMutter[nochzudurchlaufendeelementpositionen[0]]
                molekularstruktur2.append(atom)
                molekularstruktur2übersetzungsaltepositionen[nochzudurchlaufendeelementpositionen[0]] = len(molekularstruktur2)-1
                nochzudurchlaufendeelementpositionen.pop(0)
                for a in range(1, len(atom)):
                    if not atom[a][1] in schondurchlofeneelementepositionen:
                        schondurchlofeneelementepositionen.append(atom[a][1])
                        nochzudurchlaufendeelementpositionen.append(atom[a][1])

            for a in range(längemolekularstruktur2,len(molekularstruktur2)):
                for b in range(1, len(molekularstruktur2[a])):
                    molekularstruktur2[a][b][1] = molekularstruktur2übersetzungsaltepositionen[molekularstruktur2[a][b][1]]


            # die restlichen Elemente gehen zur molekularstruktur2
            for position, atom in enumerate(zerschnittenemolekularstrukturMutter):
                if not position in schondurchlofeneelementepositionen:
                    molekularstruktur1.append(atom)
                    molekularstruktur1übersetzungsaltepositionen[position] = len(molekularstruktur1)-1


            for a in range(längemolekularstruktur1,len(molekularstruktur1)):
                for b in range(1, len(molekularstruktur1[a])):
                    molekularstruktur1[a][b][1] = molekularstruktur1übersetzungsaltepositionen[molekularstruktur1[a][b][1]]

            return molekularstruktur1,molekularstruktur2




            
            
            









            return molekularstruktur1,molekularstruktur2

        def Molekularstrukturflicken(molekularstruktur):

            def AtomausMolekularstrukturlöschen(molekularstruktur,position):
                molekularstruktur.pop(position)
                for atom in molekularstruktur:
                    for a in range(len(atom)-1,0,-1):
                        if atom[a][1] == position:
                            atom.pop(a)
                        elif atom[a][1] > position:
                            atom[a][1] -= 1
                return molekularstruktur



            #zuerst muss geschaut werden, dass alle Elemente in der richtigen Häufigkeit auftritt.
            elemente = self.elemente.copy()
            if not self.testfunktiondersymetriebeidenBindungen(molekularstruktur):
                raise Exception("Bindungsymetrie ist nicht gegeben")




            for atom in molekularstruktur:
                elemente[atom[0]] -= 1

            #print(elemente)

            anzahlelementewelchemanhinzufügenkann = sum(elemente)
            #zuerst wird versucht 4 wertige elemente zu generieren und zu löschen, damit die Anzahl stimmt
            #danach wird versucht 3 wertige Elemente mit anderen 3 wertigen Elemente zu ersetzen und falls möglich zuletzt noch


            #Muss noch auf Spezialfälle überprüft werden

            #print("ich bin im Loop gefangen Molekularstrukturflicken")
            count = 0
            while count < len(elemente):
                #print(elemente)
                #print(count)
                #print("Elementewelche man hinzufügen kann : " + str(anzahlelementewelchemanhinzufügenkann))
                if elemente[count] == 0:
                    count += 1

                #zuerst muss noch der Spezialfall geklärt werden, wenn count beim hintersten Glied angekommen ist

                elif count == len(elemente) -1:
                    if elemente[count] > 0:
                        molekularstruktur.append([count])
                        elemente[count] -= 1
                        anzahlelementewelchemanhinzufügenkann -= 1
                    elif elemente[count] < 0:
                        for a in range(len(molekularstruktur)):
                            if molekularstruktur[a][0] == count:
                                molekularstruktur = AtomausMolekularstrukturlöschen(molekularstruktur,a)
                                break
                        elemente[count] += 1
                        anzahlelementewelchemanhinzufügenkann += 1

                elif all(a >= 0 for a in elemente): # WEnn man nichts mehr tauschen kann, sondern nur noch hinzufügen
                    molekularstruktur.append([count])
                    elemente[count] -= 1
                    anzahlelementewelchemanhinzufügenkann -= 1

                elif all(a <= 0 for a in elemente): # wenn man nichts mehr tauschen kann, sondern nur noch löschen kann
                    for a in range(len(molekularstruktur)):
                        if molekularstruktur[a][0] == count:
                            molekularstruktur = AtomausMolekularstrukturlöschen(molekularstruktur,a)
                            break
                    elemente[count] += 1
                    anzahlelementewelchemanhinzufügenkann += 1

                #Normalfall
                elif elemente[count] > 0:
                    for position,wert in enumerate(elemente):
                        if wert < 0:
                            if self.Anzahlverbindungen(position) < self.Anzahlverbindungen(count) and anzahlelementewelchemanhinzufügenkann > 0:
                                molekularstruktur.append([count])
                                elemente[count] -= 1
                                anzahlelementewelchemanhinzufügenkann -= 1
                            else:
                                for atom in molekularstruktur:
                                    if atom[0] == position:
                                        atom[0] = count
                                        elemente[count] -= 1
                                        elemente[position] += 1
                                        break
                            break
                elif elemente[count] < 0:
                    for position,wert in enumerate(elemente):
                        if wert > 0:
                            if self.Anzahlverbindungen(position) < self.Anzahlverbindungen(count) and anzahlelementewelchemanhinzufügenkann < 0:
                                for a in range(len(molekularstruktur)):
                                    if molekularstruktur[a][0] == count:
                                        molekularstruktur = AtomausMolekularstrukturlöschen(molekularstruktur,a)
                                        break
                                elemente[count] += 1
                                anzahlelementewelchemanhinzufügenkann += 1
                            for atom in molekularstruktur:
                                if atom[0] == count:
                                    atom[0] = position
                                    elemente[count] += 1
                                    elemente[position] -= 1
                                    break
                            break


            #Nun müssen die offenen Verbindungen wieder geflickt werden

            # um die anderen Funktionen zu nutzen muss die molekularstruktur zuerst in self.molekularstruktur umgewandelt werden


            if not self.testfunktiondersymetriebeidenBindungen(molekularstruktur):
                raise Exception("Bindungsymetrie ist nicht gegeben2")



            self.molekularstruktur = molekularstruktur


            #Zu aller erst muss geaschaut werden, dass die DBÄ nicht zu hoch ist, andernfalls müssen solange Mehrfachbindungen gesplitet werden, bis es aufgeht.
            anzahlübrigerDoppelbindungsen = (self.elemente[0] * 2 + self.elemente[1] - self.elemente[5] - self.elemente[6] - self.elemente[7] - self.elemente[8]) / 2 + 1

            for atom in self.molekularstruktur:
                for a in range(1,len(atom)):
                    anzahlübrigerDoppelbindungsen -= (atom[a][0]-1)/2

            if anzahlübrigerDoppelbindungsen % 1 != 0:
                raise Exception("anzahlübrigeDoppelbindungen ist nicht eine ganze Zahl")

            while anzahlübrigerDoppelbindungsen < 0:
                for position,atom in enumerate(self.molekularstruktur):
                    for a in range(1, len(atom)):
                        if atom[a][0] > 1 and anzahlübrigerDoppelbindungsen < 0:
                            wertigkeit = atom[a][0]
                            position2 = atom[a][1]
                            self.Deletverbindung(position,position2)
                            self.Creatverbindung(position,position2,wertigkeit-1)
                            anzahlübrigerDoppelbindungsen += 1
                            break



            # Es kann nun sein, dass die molekularstruktur in vielen verschiedenen einzelstücke fragmentiert ist, diese Einzelstücke müssen wieder miteinander Verbunden werden

            offenestellen = []
            for position,atom in enumerate(self.molekularstruktur):
                if self.Anzahloffeneverbindungen(atom) != 0:
                    offenestellen.append([self.Anzahloffeneverbindungen(atom),position])

            #print("Molekularstruktur")
            #print(molekularstruktur)
            #print(self.elemente)

            # untersucht ob zwei Teile nicht miteinander verbunden sind und verbindet sie dann
            """
            print("offene Stellen")
            print(offenestellen)
            print("Molekularstruktur")
            print(molekularstruktur)
            print(self.molekularstruktur)
            
            """

            for a in range(len(offenestellen)):
                for b in range(a+1,len(offenestellen)):
                    if offenestellen[a][0] > 0 and offenestellen[b][0] > 0:
                        if not self.Isindirektverbunden(offenestellen[a][1],offenestellen[b][1]):
                            self.Creatverbindung(offenestellen[a][1],offenestellen[b][1],1)
                            offenestellen[a][0] -= 1
                            offenestellen[b][0] -= 1

            # Alle werden gelöscht, welche keine offene Stellen mehr haben

            for a in range(len(offenestellen)-1,-1,-1):
                if offenestellen[a][0] == 0:
                    offenestellen.pop(a)




            """
            # Nun werden versucht die restlichen 3fach Verbindungenen zu füllen
            if anzahlübrigerDoppelbindungsen > 1:
                for a in range(len(offenestellen)):
                    for b in range(a+1,len(offenestellen)):
                        if offenestellen[a][0] > 2 and offenestellen[b][0] > 2 and anzahlübrigerDoppelbindungsen >= 2 and self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstruktur[offenestellen[a][1]],offenestellen[b][1]) == 0:
                            self.Creatverbindung(offenestellen[a][1],offenestellen[b][1],3)
                            offenestellen[a][0] -= 3
                            offenestellen[b][0] -= 3
                            anzahlübrigerDoppelbindungsen -= 2
            
            # Nun werden versucht die doppelbindungen zu füllen
            
            if anzahlübrigerDoppelbindungsen >= 1:
                for a in range(len(offenestellen)):
                    for b in range(a+1,len(offenestellen)):
                        if offenestellen[a][0] > 1 and offenestellen[b][0] > 1 and anzahlübrigerDoppelbindungsen >= 1 and self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstruktur[offenestellen[a][1]],offenestellen[b][1]) < 2:
                            self.Creatverbindung(offenestellen[a][1],offenestellen[b][1],2)
                            offenestellen[a][0] -= 2
                            offenestellen[b][0] -= 2
                            anzahlübrigerDoppelbindungsen -= 1
            
            #  nun werden versucht die restlichen einfachbindungen zu füllen
            for a in range(len(offenestellen)):
                for b in range(a + 1, len(offenestellen)):
                    if offenestellen[a][0] > 1 and offenestellen[b][0] > 1 and anzahlübrigerDoppelbindungsen >= 1 and self.AnzahlverbindungenzwischenzweiAtomen(self.molekularstrukturolekularstruktur[offenestellen[a][1]], offenestellen[b][1]) < 2:
                        self.Creatverbindung(offenestellen[a][1], offenestellen[b][1], 2)
                        offenestellen[a][0] -= 2
                        offenestellen[b][0] -= 2
                        anzahlübrigerDoppelbindungsen -= 1
            
            """

            #Kopie von Strukturinitalisierung

            # Um nicht den ganzen Code neu zu schreiben muss der Name des anzahlübrigenDoppelbindungsen in Doppelbindungsequivalenz umgewandelt werden und self.molekularstruktur zurück zu molekularstruktur

            molekularsturktur = self.molekularstruktur
            Doppelbindungsequivalenz = anzahlübrigerDoppelbindungsen

            if Doppelbindungsequivalenz > 0:

                restlicheoffenestellen = []
                for atom in molekularsturktur:
                    restlicheoffenestellen.append(self.Anzahloffeneverbindungen(atom))

                # zuerst dreifachoffene stellen werden gesucht
                positionen = []

                for count, wert in enumerate(restlicheoffenestellen):
                    if wert == 3:
                        positionen.append(count)

                if positionen != []:

                    a = 0
                    while a < len(positionen):
                        for b in range(a + 1, len(positionen)):
                            if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],
                                                                         positionen[b]) == 0:
                                molekularsturktur[positionen[a]].append([3, positionen[b]])
                                molekularsturktur[positionen[b]].append([3, positionen[a]])
                                Doppelbindungsequivalenz -= 3
                                positionen.pop(b)
                                break
                        if self.Anzahloffeneverbindungen(molekularsturktur[positionen[a]]) == 0:
                            positionen.pop(a)
                        else:
                            a += 1

                # zweifachoffene Stellen versuchen zu schliessen
                for count, wert in enumerate(restlicheoffenestellen):
                    if wert == 2:
                        positionen.append(count)

                if positionen != []:
                    a = 0
                    while a < len(positionen):
                        for b in range(a + 1, len(positionen)):
                            if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],
                                                                         positionen[b]) == 0:
                                molekularsturktur[positionen[a]].append([2, positionen[b]])
                                molekularsturktur[positionen[b]].append([2, positionen[a]])
                                Doppelbindungsequivalenz -= 2
                                positionen.pop(b)
                                break
                        if self.Anzahloffeneverbindungen(molekularsturktur[positionen[a]]) <= 1:
                            positionen.pop(a)
                        else:
                            a += 1

                # die einfachen STellen werden am schluss auch noch versucht zu schliessen
                positionen = []
                for count, wert in enumerate(restlicheoffenestellen):
                    if wert >= 1:
                        positionen.append(count)

                if positionen != []:

                    a = 0
                    while a < len(positionen):
                        for b in range(a + 1, len(positionen)):

                            if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],
                                                                         positionen[b]) == 0:
                                molekularsturktur[positionen[a]].append([1, positionen[b]])
                                molekularsturktur[positionen[b]].append([1, positionen[a]])
                                Doppelbindungsequivalenz -= 1
                                if self.Anzahloffeneverbindungen(molekularsturktur[positionen[b]]) == 0:
                                    positionen.pop(b)
                                break
                            elif self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],
                                                                           positionen[b]) <= 2:

                                for c in range(1, len(molekularsturktur[positionen[a]])):
                                    if molekularsturktur[positionen[a]][c][1] == positionen[b]:
                                        molekularsturktur[positionen[a]][c][0] += 1

                                for c in range(1, len(molekularsturktur[positionen[b]])):
                                    if molekularsturktur[positionen[b]][c][1] == positionen[a]:
                                        molekularsturktur[positionen[b]][c][0] += 1

                                if self.Anzahloffeneverbindungen(molekularsturktur[positionen[b]]) == 0:
                                    positionen.pop(b)
                                break

                        if self.Anzahloffeneverbindungen(molekularsturktur[positionen[a]]) == 0:
                            positionen.pop(a)
                        else:
                            a += 1

                # es kann leider am Schluss noch sein, dass eine Molekül noch zwei freie Verbindungen hat, diese werden jetzt noch mit einer Art der Elektophilen Additon geschlossen
                for count, atom in enumerate(molekularsturktur):
                    if (self.Anzahloffeneverbindungen(atom)) == 2:
                        # Nun wird versucht eine Doppelbindung oder dreifachbindung um eins zu reduzieren und dann mit den neuen zwei offnene stellen das Atom damit zu verbinden
                        breakbool = False
                        for a in range(len(molekularsturktur)):
                            for b in range(1, len(molekularsturktur[a])):
                                if molekularsturktur[a][b][0] >= 2 and a != count:
                                    breakbool = True
                                    molekularsturktur[a][b][0] -= 1
                                    for c in range(1, len(molekularsturktur[molekularsturktur[a][b][1]])):
                                        if molekularsturktur[molekularsturktur[a][b][1]][c][1] == a:
                                            molekularsturktur[molekularsturktur[a][b][1]][c][0] -= 1
                                            break
                                    Doppelbindungsequivalenz -= 1

                                    if self.AnzahlverbindungenzwischenzweiAtomen(atom, a) == 0:
                                        molekularsturktur[count].append([1, a])
                                        molekularsturktur[a].append([1, count])
                                    else:
                                        for c in range(1, len(molekularsturktur[a])):
                                            if molekularsturktur[a][c][1] == count:
                                                molekularsturktur[a][c][0] += 1
                                                break
                                    if self.AnzahlverbindungenzwischenzweiAtomen(atom, molekularsturktur[a][b][1]) == 0:
                                        molekularsturktur[count].append([1, molekularsturktur[a][b][1]])
                                        molekularsturktur[molekularsturktur[a][b][1]].append([1, count])

                                    else:
                                        for c in range(1, len(molekularsturktur[molekularsturktur[a][b][1]])):
                                            if molekularsturktur[a][c][1] == count:
                                                molekularsturktur[a][c][0] += 1
                                                break

                            if breakbool:
                                break

            #if Doppelbindungsequivalenz // 1 == Doppelbindungsequivalenz and Doppelbindungsequivalenz >= 0:
                #for atom in molekularsturktur:
                    #if self.Anzahloffeneverbindungen(atom) != 0:
                        #print("Fehler bei Molekularstrukturflicken")

            return molekularsturktur



















            return molekularstruktur

        anzahlversuchteKinder = 15 #Es gibt immer die doppelte Anzahl an Kinder, als der Wert dieses Integer
        """
        print("Molekularstrukturen der Eltern")
        print(molekularstrukturVater)
        print(molekularstrukturMutter)
        """

        zerschnittenemolekularstrukturenMutter = []
        zerschnittenemolekularstrukturenVater = []

        count = 0
        for a in range(1000): # Molekülstrukturzerschneiden() liefert nicht immer ein Resultat, muss man es manchmal mehrmals durchlaufen lassen für eine zerschnittene Molekülstruktur. Es funktioniert etwa in 80% der fällen
            manipuliertemolekularstruktur = Molekularstrukturzerschneider(molekularstrukturVater)


            if  manipuliertemolekularstruktur != False:
                zerschnittenemolekularstrukturenVater.append(copy.deepcopy(manipuliertemolekularstruktur))
                count += 1
                if count == anzahlversuchteKinder:
                    break

        count = 0
        for a in range(1000):
            manipuliertemolekularstruktur = Molekularstrukturzerschneider(molekularstrukturMutter)
            if manipuliertemolekularstruktur != False:
                zerschnittenemolekularstrukturenMutter.append(copy.deepcopy(manipuliertemolekularstruktur))
                count += 1
                if count == anzahlversuchteKinder:
                    break


        #die beiden Vater und Mutter Listen werden gleich lang gemacht
        while(len(zerschnittenemolekularstrukturenVater) != len(zerschnittenemolekularstrukturenMutter)):
            if len(zerschnittenemolekularstrukturenVater) > len(zerschnittenemolekularstrukturenMutter):
                zerschnittenemolekularstrukturenVater.pop(0)
            else:
                zerschnittenemolekularstrukturenMutter.pop(0)



        zerschnittenemolekularstrukturenVater.sort(key=lambda x: x[1])
        zerschnittenemolekularstrukturenMutter.sort(key=lambda x: x[1])
        """
        for test in zerschnittenemolekularstrukturenVater:
            print(test[1])

        for test in zerschnittenemolekularstrukturenMutter:
            print(test[1])
        """

        neuemolekularstrukturen = []
        for a in range(len(zerschnittenemolekularstrukturenVater)):
            childs = Molekularstrukturzusammenfügen(zerschnittenemolekularstrukturenMutter[a][0],zerschnittenemolekularstrukturenVater[a][0])
            neuemolekularstrukturen.append(childs[0])
            neuemolekularstrukturen.append(childs[1])






        # Alle childs werden jetzt geflickt und dann anschliessend noch geschaut ob sie wirkliches Molekül sind und falls nicht, wird es nicht gespeichert
        children = []
        for molekularstruktur in neuemolekularstrukturen:
            child = Molekularstrukturflicken((molekularstruktur))
            self.molekularstruktur = child
            if self.isligit():
                children.append(child)


        if children == []:
            self.molekularstruktur = None
            print("es gab kein Child")
            return False

        #print("Es gibt ein Kind")


        # Mithilfe der Heurist wird jetzt das beste Kind ausgewählt und wird zur neuen molekularstruktur
        kleinstepunktzahl = 100000000000000
        positionchild = 0
        for a,child in enumerate(children):
            self.CalcHeuristik()
            if self.heuristikwert < kleinstepunktzahl:
                kleinstepunktzahl = self.heuristikwert
                positionchild = a

        self.molekularstruktur = copy.deepcopy(children[positionchild])
        """
        print("endgültige Molekularstruktur")
        print(self.molekularstruktur)
        """
        return True








    def __init__(self, gruppenkonfiguration = None, anzahlmutationen = 0, molekularstrukturvater = None, molekularstrukturmutter = None, elemente = None,molekularstruktur = None):
        if molekularstrukturvater == None and molekularstrukturvater == None and gruppenkonfiguration == None: # Erschaffung eines Dublikanten und Mutation desen
            self.molekularstruktur = molekularstruktur
            self.elemente = elemente
            self.elementgruppengrenzen = [0, 1, 2, 5]
            self.elementmassen = [12, 13, 14, 28, 16, 15, 17, 29, 45]
            for a in range(anzahlmutationen):
                self.Muation()


        elif molekularstrukturmutter == None and molekularstrukturvater == None: # Startinitialisation
        # 0 = C / 1 = CH / 2 = CH2 / 3 = Keton / 4 = Ether / 5 = CH3 / 6 = OH / 7 = Aldehyd / 8 = Carbonsäure
            self.gruppenkonfiguration = gruppenkonfiguration
            self.elemente = [0] * len(gruppenkonfiguration)
            self.elemente[0] = gruppenkonfiguration[0]
            self.elemente[1] = gruppenkonfiguration[1]
            self.elemente[2] = gruppenkonfiguration[2]
            self.elemente[3] = gruppenkonfiguration[6]
            self.elemente[4] = gruppenkonfiguration[8]
            self.elemente[5] = gruppenkonfiguration[3]
            self.elemente[6] = gruppenkonfiguration[4]
            self.elemente[7] = gruppenkonfiguration[5]
            self.elemente[8] = gruppenkonfiguration[7]
            self.elementgruppengrenzen = [0,1,2,5]
            self.elementmassen = [12,13,14,28,16,15,17,29,45]

            Doppelbindungsequivalenz = (self.elemente[0]*2 + self.elemente[1] -self.elemente[5] - self.elemente[6] -self.elemente[7] -self.elemente[8])/2 + 1
            if not (Doppelbindungsequivalenz // 1 == Doppelbindungsequivalenz and Doppelbindungsequivalenz >= 0):
                print(self.elemente)
                print("FEHLERMELDUNG, dieses Molekül exisistiert aufgrund einer nicht möglichen Doppelbindungswquivalenz nicht.")
                return


            self.Strukturintertialgenerator(anzahlmutationen)

        else: # Erschaffung durch zwei Eltern
            self.elemente = elemente
            self.elementgruppengrenzen = [0, 1, 2, 5]
            self.elementmassen = [12, 13, 14, 28, 16, 15, 17, 29, 45]
            self.ChildauszweiParents(molekularstrukturvater, molekularstrukturmutter)

        self.smistring = None



        #print("überprüfung :" +str(self.isligit()))
        #print(self.molekularstruktur)

    def testfunktion(self, molekularstruktur1, molekularstruktur2):
        for a in range(len(molekularstruktur1)):
            if molekularstruktur1[a] != molekularstruktur2[a]:
                print("nichtgleich")
                print(a)
                print(molekularstruktur1[a])
                print(molekularstruktur2[a])

    def testfunktiondersymetriebeidenBindungen(self,molekularstruktur):

        exist = False

        for position,atom in enumerate(molekularstruktur):
            for a in range(1,len(atom)):
                atom2 = molekularstruktur[atom[a][1]]
                for b in range(1,len(atom2)):
                    if position == atom2[b][1] and atom[a][0] == atom2[b][0]:
                        exist = True
                        break
                if not exist:
                    return False
                exist = False
        return True



















