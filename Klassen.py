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
        molekularsturktur = [0]* sum(self.elemente)
        verfuegbareelemente = self.elemente

        dbe = (self.elemente[0]*2 + self.elemente[1] - self.elemente[5] - self.elemente[6] - self.elemente[7] - self.elemente[8])
        maxdreifachbindungen = dbe // 2



        for element in self.elemente:
            for a in range(element):
                molekularsturktur.append(element)

        #Einsetzen der dreifachbindungen
        anzahlKohlenstoffatome = self.elemente[0]

        if anzahlKohlenstoffatome // 2 > maxdreifachbindungen:
            anzahlKohlenstoffatome = maxdreifachbindungen

        for a in range(anzahlKohlenstoffatome // 2):
            molekularsturktur[2*a].append([3,2*a+1])
            molekularsturktur[2*a +1].append([3,2*a])

        anzahl









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




        
        






