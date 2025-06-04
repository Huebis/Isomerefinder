import subprocess

class Molekuelinfo():

    # Anzahl C / Anzahl 0/ Anzahl H
    isomere = [0,0,0]
    # Anzahl Alkohole / Anzahl Aldehyde / Anzahl Ketone / Anzahl Carbonsäuren / Anzahl Ester + Carbonsäuren / Anzahl Carbonsäuren + Alkohole
    oxygeniumsubstitution = [0,0,0,0,0,0]
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

        result = subprocess.run('java -jar /home/eliahh/Workspace/Matura/test/MAYGEN/target/MAYGEN-1.8.jar -f "'+summenformelstring +'" -setElements -v -t -smi -o /home/eliahh/Workspace/Matura/test/Isomergruppen', shell=True, capture_output=True, text=True)

        print(result)
        
        

        
        


    def test(self):
        print("isomere:", self.isomere)
        print("oxygenium:", self.oxygeniumsubstitution)
        print("Carbonsubstitutionsgrad:", self.Carbonsubstitutionsgrad)
        print("gruppenkonfiguration", self.gruppenkonfiguration)
