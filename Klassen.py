import subprocess
import random
import copy

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

    def __init__(self):
        # Ein mutables Objekt wie eine Liste wird in allen Instanzen geteilt, wenn es im Konstruktor nicht explizit initialisiert wird
        # self.CarbonID = [0] * 5
        self.gruppenkonfiguration = None

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



    def isligit(self, anzahlverbindungenüberspringen = False, anazahlAtomeüberspringen = False, topologischzusammenhängendüberspringen = False):

        def AlleAtomehabenalleverbindungenverbraucht():
            for atom in self.molekularstruktur:
                if self.Anzahloffeneverbindungen(atom) != 0:
                    print("Error nicht alle Atome haben die richtige Anzahl an Verbindungen" + str(atom))
                    return False
            return True

        def AlleAtomekommeninderrichtigneAnzahlvor():
            elemente = [0 for a in range(len(self.elemente))]
            for atom in self.molekularstruktur:
                elemente[atom[0]] += 1
            if elemente == self.elemente:
                return True

            print("Error es kommen nicht alle Atome in der richtigen Anzahl vor")
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
                    print("Error das Atom ist nicht zusammenhängend")
                    return False

            return True





        if not AlleAtomehabenalleverbindungenverbraucht() and not anzahlverbindungenüberspringen :
            return False
        if not AlleAtomekommeninderrichtigneAnzahlvor() and not anazahlAtomeüberspringen:
            return False
        if not Topologischzusammenhängend() and not topologischzusammenhängendüberspringen:
            return False
        return True



    def CalcHeuristik(self):
        return False
    def Lokalesuche(self):
        return False
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

            print("liste aller möglichkeiten")
            print(listeallerMehrfachbindungen)
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

            print("Substitutionsmöglichkeiten")
            print(substitutionsmöglichkeiten)

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
            print(substitutionsmöglichkeiten)



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

                    self.Deletverbindung(positionatom, positionatomeinfachbindung)
                    self.Deletverbindung(positionatom, positionatomdoppelbindung)
                    self.Creatverbindung(positionatomeinfachbindung, positionatomdoppelbindung, 1)
                    self.Creatverbindung(positionatomdoppelbindung, positionatom, 1)

                    positionverbleibendesatom = self.molekularstruktur[positionatom][1][
                        1]  # es kann sein, dass mit diesem Atom eine neue Bindung eingeganen werden soll, dann muss es aber anstatt einer neuer einfachbindung eine doppelbindung geben

                    # das herausgelöste Atom wird nun wieder an bestimmten ort wieder eingesetzt
                    positionsubstitiontsatom1 = möglichesubstitutionen[substitutionsverbindung][0]
                    positionsubstitiontsatom2 = self.molekularstruktur[positionsubstitiontsatom1][möglichesubstitutionen[substitutionsverbindung][1]][1]
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
            print("MutationzweiergleichwertigerAtome")
            for a in range(10):
                if MutationzweiergleichwertigerAtome():
                    return True



        elif randommutationszeiger < 0.7:
            print("MutationAtomherausreisenundanderswoneueinsetzen ")
            for a in range(10):
                if MutationAtomherausreisenundanderswoneueinsetzen():
                    return True
        else:
            print("MutationverschiebungMehrfachbindung")
            for a in range(2):
                zwischenspeichermolekularstruktur = copy.deepcopy(self.molekularstruktur)
                if MutationreduzierungMehrfachbindungzuanderemOrt():
                    # Es gibt leider sehr viele Spezialfälle, bei welchem diese Mutation nicht geht. Diese sind aber meistens so speziel, dass es sich nicht lohnt, diese manuel einzelt herauszufiltern. Daher überprüfe ich das Molekül und wenn es nicht mehr zusammenhängend ist und falls dies der Fall ist, muss ich alles wieder rückgängig machen.
                    if self.isligit():
                        return True
                    else:
                        self.molekularstruktur = copy.deepcopy(zwischenspeichermolekularstruktur)



        return False

        
        
        
        

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
                    return "(C(O))"
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
                return 6

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
        ortdergemachtenAtomgruppen[0] = [0,LängederÜbersetzung(0)]
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


    def DarstellungMolekülinSMI(self, mitWasserstoff = False, mitBild = True):
        if self.smistring != None:

            mol = Chem.MolFromSmiles(self.smistring)
            if mitWasserstoff:
                mol = Chem.AddHs(mol)
            rdkit.Chem.AllChem.Compute2DCoords(mol)
            if mitBild:
                img = Draw.MolToImage(mol, size=(2000, 2000), kekulize=True, bgcolor=(255, 255, 255))
                img.show()
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


        return molekularsturktur

        return False

    def __init__(self, gruppenkonfiguration, anzahlmutationen):

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

        Doppelbindungsequivalenz = (self.elemente[0]*2 + self.elemente[1] -self.elemente[5] - self.elemente[6] -self.elemente[7] -self.elemente[8])/2 + 1
        if not (Doppelbindungsequivalenz // 1 == Doppelbindungsequivalenz and Doppelbindungsequivalenz >= 0):
            print(self.elemente)
            print("FEHLERMELDUNG, dieses Molekül exisistiert aufgrund einer nicht möglichen Doppelbindungswquivalenz nicht.")
            return


        self.molekularstruktur = self.Strukturintertialgenerator(anzahlmutationen)
        self.smistring = None



        print("überprüfung :" +str(self.isligit()))
        print(self.molekularstruktur)

    def testfunktion(self):
        for atom in self.molekularstruktur:
            wertigkeitatom = self.Anzahlverbindungen(atom[0])

            if wertigkeitatom == 3:
                print(len(atom))






test = individuum([2,6 , 2, 2, 2, 2, 2, 0, 0],0)


test.SMilestransformator()
test.DarstellungMolekülinSMI(False, True )


for a in range(10000):
    test.Muation()
    if not test.isligit():
        print("FEHLER")
        print(test.molekularstruktur)
        raise ValueError("Ungültige Molekularstruktur")


test.SMilestransformator()
test.DarstellungMolekülinSMI(False, True)
#test.testfunktion()






#for a in range(100000  ):
    #test = individuum([random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20)], 0)
        
print("finish")
# testfehlschlänge

#[5,0,0,0,0,0,0,0,0]












