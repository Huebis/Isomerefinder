import subprocess

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
    def CalcHeuristik(self):
        return False
    def Lokalesuche(self):
        return False
    def Muation(self):
        return False
    def SMilestransformator(self):
        return False
    def Strukturintertialgenerator(self,anzahlmutationen):
        molekularsturktur = [[0] * 1 for _ in range(sum(self.elemente))]

        def Anzahlverbindungen(elementnummer):
            #print(elementnummer)
            if elementnummer >=  self.elementgruppengrenzen[3]:
                return 1
            if elementnummer >= self.elementgruppengrenzen[2]:
                return 2
            if elementnummer >= self.elementgruppengrenzen[1]:
                return 3
            if elementnummer >= self.elementgruppengrenzen[0]:
                return 4

        def Anzahloffeneverbindungen(Atomgruppe):
            anzahlmöglicheverbindungen = Anzahlverbindungen(Atomgruppe[0])


            anzahlverbrauchteverbindungen = 0
            for a in range(1,len(Atomgruppe),1):
                anzahlverbrauchteverbindungen += Atomgruppe[a][0]

            return anzahlmöglicheverbindungen - anzahlverbrauchteverbindungen


        count = 0
        for a in range(len(self.elemente)):
            for b in range(self.elemente[a]):
                molekularsturktur[count][0] = a
                count += 1


        #jetzt kann das wirkliche zusammenbauen eines provisorischen Moleküles beginnen

        offeneverbindungen = Anzahloffeneverbindungen(molekularsturktur[0])
        Doppelbindungsequivalenz = (self.elemente[0]*2 + self.elemente[1] -self.elemente[5] - self.elemente[6] -self.elemente[7] -self.elemente[8])/2

        for a in range(1,len(molekularsturktur),1):
            maxoffeneverbindungen = 0
            ortmax = 0
            for b in range(a):
                if Anzahloffeneverbindungen(molekularsturktur[b]) > maxoffeneverbindungen:
                    maxoffeneverbindungen = Anzahloffeneverbindungen(molekularsturktur[b])
                    ortmax = b


            #print(a)
            #print(Doppelbindungsequivalenz)


            if Anzahlverbindungen(molekularsturktur[a][0]) == 4:
                if maxoffeneverbindungen >= 3 and Doppelbindungsequivalenz >= 2 and offeneverbindungen >= 3:
                    molekularsturktur[a].append([3,ortmax])
                    molekularsturktur[ortmax].append([3,a])
                    offeneverbindungen -= 2
                    Doppelbindungsequivalenz -= 2

                elif maxoffeneverbindungen >= 2 and Doppelbindungsequivalenz >= 1 and offeneverbindungen:
                    molekularsturktur[a].append([2, ortmax])
                    molekularsturktur[ortmax].append([2, a])
                    Doppelbindungsequivalenz -= 1
                    #Anzahl der offenenverbindungen bleibt gleich
                elif maxoffeneverbindungen >= 1:
                    molekularsturktur[a].append([1, ortmax])
                    molekularsturktur[ortmax].append([1, a])
                    offeneverbindungen += 2

            if Anzahlverbindungen(molekularsturktur[a][0]) == 3:
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

            if Anzahlverbindungen(molekularsturktur[a][0]) == 2:

                if maxoffeneverbindungen >= 2 and Doppelbindungsequivalenz >= 1 and offeneverbindungen >= 3:
                    molekularsturktur[a].append([2, ortmax])
                    molekularsturktur[ortmax].append([2, a])
                    Doppelbindungsequivalenz -= 1
                    offeneverbindungen -= 2

                elif maxoffeneverbindungen >= 1:
                    molekularsturktur[a].append([1, ortmax])
                    molekularsturktur[ortmax].append([1, a])

            if Anzahlverbindungen(molekularsturktur[a][0]) == 1:
                molekularsturktur[a].append([1, ortmax])
                molekularsturktur[ortmax].append([1, a])
                offeneverbindungen -= 1


            if a == 5:
                print(offeneverbindungen)
                print(molekularsturktur)
                print(Doppelbindungsequivalenz)
                print(maxoffeneverbindungen)
                print(ortmax)
                print(Anzahloffeneverbindungen(molekularsturktur[0]))
                print(Anzahlverbindungen(molekularsturktur[0][0]))



        print(molekularsturktur)
        return molekularsturktur














        return False

    def __init__(self, gruppenkonfiguration, anzahlmutationen):

        # 0 = C / 1 = CH / 2 = CH2 / 3 = Keton / 4 = Ether / 5 = CH3 / 6 = OH / 7 = Aldehyd / 8 = Carbonsäure
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

test = individuum([1,2,5,0,0,0,2,0,0],0)


        
        






