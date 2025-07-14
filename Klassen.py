import subprocess
import random

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


    def isligit(self, anzahlverbindungenüberspringen = False, anazahlAtomeüberspringen = False, topologischzusammenhängendüberspringen = False):

        def AlleAtomehabenalleverbindungenverbraucht():
            for atom in self.molekularstruktur:
                if self.Anzahloffeneverbindungen(atom) != 0:
                    return False
            return True

        def AlleAtomekommeninderrichtigneAnzahlvor():
            elemente = [0 for a in range(len(self.elemente))]
            for atom in self.molekularstruktur:
                elemente[atom[0]] += 1
            if elemente == self.elemente:
                return True
            return False

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
                    return False

            return True





        if not AlleAtomehabenalleverbindungenverbraucht() or anzahlverbindungenüberspringen :
            return False
        if not AlleAtomekommeninderrichtigneAnzahlvor() or anazahlAtomeüberspringen:
            return False
        if not Topologischzusammenhängend() or topologischzusammenhängendüberspringen:
            return False
        return True



    def CalcHeuristik(self):
        return False
    def Lokalesuche(self):
        return False
    def Muation(self):
        #Alle Mutationsfunktionen returnen True wenn eine Mutation durchgeführt werden konnte und False, wenn keine Mutation durchgeführt wurde. Die Manipulation wird direkt am self.molekularstruktur abgeändert


        def MutationzweiergleichwertigerAtome(spezifischeposition = None):

            for a in range(10):
                if spezifischeposition == None:
                    positionatom1 = random.randint(self.elementgruppengrenzen[1], len(self.molekularstruktur)-1)
                else:
                    positionatom1 = spezifischeposition


                wertigkeitatom1 = self.Anzahlverbindungen(self.molekularstruktur[positionatom1][0])
                gleichwertigeatome = []
                for count, atom in enumerate(self.molekularstruktur):
                    if atom[0] != self.molekularstruktur[positionatom1][0] and self.Anzahlverbindungen(atom[0]) == wertigkeitatom1:
                        gleichwertigeatome.append(count)

                if gleichwertigeatome != []:
                    positionatom2 = gleichwertigeatome[random.randint(0,len(gleichwertigeatome)-1)]
                    atom1 = self.molekularstruktur[positionatom1][0]
                    atom2 = self.molekularstruktur[positionatom2][0]
                    self.molekularstruktur[positionatom1][0] = atom2
                    self.molekularstruktur[positionatom2][0] = atom1
                    return True
            return False

        def MutationreduzierungMehrfachbindungzuanderemOrt():
            return False

        def MutationreduzierungMehrfachbindungzueinemCyclo():
            return False

        def Muationentfernungcyclofüreinemehrfachbindung():
            return False

        def MutationAtomherausreisenundanderswoneueinsetzen():
            for a in range(5):
                positionatom = random.randint(0,len(self.molekularstruktur)-1)
                wertigkeitatom = self.Anzahlverbindungen(self.molekularstruktur[positionatom][0])

                if wertigkeitatom == 1:
                    if MutationzweiergleichwertigerAtome(1):
                        return True
                if wertigkeitatom == 2:
                    #Suche eine einfachbindung und schmuggle dich rein
                    return False
                if wertigkeitatom == 3:
                    return False

                if wertigkeitatom == 4:
                    return False




        MutationzweiergleichwertigerAtome()
        
        
        
        

    def SMilestransformator(self):
        return False

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
                        if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]], molekularsturktur[positionen[b]]) == 0:

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
                        if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],molekularsturktur[positionen[b]]) == 0:
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
                        if self.AnzahlverbindungenzwischenzweiAtomen(molekularsturktur[positionen[a]],molekularsturktur[positionen[b]]) == 0:
                            molekularsturktur[positionen[a]].append([1, positionen[b]])
                            molekularsturktur[positionen[b]].append([1, positionen[a]])
                            Doppelbindungsequivalenz -= 1
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
                    print(self.grppenkonfiguration)

        if False:
            for atom in molekularsturktur:
                print(Anzahloffeneverbindungen(atom))
            print("Doppelbindungsequivalenz:" + str(Doppelbindungsequivalenz))
            print(molekularsturktur)


        return molekularsturktur














        return False

    def __init__(self, gruppenkonfiguration, anzahlmutationen):

        # 0 = C / 1 = CH / 2 = CH2 / 3 = Keton / 4 = Ether / 5 = CH3 / 6 = OH / 7 = Aldehyd / 8 = Carbonsäure
        self.grppenkonfiguration = gruppenkonfiguration
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


        self.molekularstruktur = self.Strukturintertialgenerator(anzahlmutationen)
        self.smistring = None

        print("überprüfung :" +str(self.isligit()))



test = individuum([4, 2, 3, 4, 0, 2, 4, 0, 2],0)


#for a in range(10000):
    #test = individuum([random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20)], 0)
        
print("finish")
# testfehlschlänge

#[5,0,0,0,0,0,0,0,0]




