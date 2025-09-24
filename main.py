import Cdept
import Formatumwandler
import MS
import Klassen
import NMR
import Math
import Testcase
import GenetischerAlgorythmus
isomerearray = Formatumwandler.Readfile()

#Formatumwandler.AddWasserstoff2("test")

#for isomere in isomerearray:
    #print(isomere)
    #Formatumwandler.AddWasserstoff2(isomere)


Klassen.Molekuelinfo.msmainpeak = 122
Klassen.Molekuelinfo.msdata = [
    [18.0, 2.8],
    [26.0, 1.2],
    [27.0, 3.6],
    [28.0, 1.4],
    [37.0, 3.6],
    [37.5, 1.4],
    [38.0, 6.2],
    [38.5, 2.4],
    [39.0, 7.2],
    [44.0, 1.2],
    [45.0, 3.2],
    [47.0, 1.2],
    [49.0, 2.6],
    [50.0, 19.9],
    [51.0, 36.5],
    [52.0, 6.2],
    [52.5, 2.0],
    [53.0, 1.4],
    [61.0, 1.3],
    [63.0, 1.4],
    [65.0, 2.7],
    [66.0, 2.1],
    [73.0, 1.9],
    [74.0, 7.2],
    [75.0, 4.1],
    [76.0, 5.9],
    [77.0, 75.3],
    [78.0, 8.1],
    [94.0, 3.2],
    [104.0, 1.1],
    [105.0, 100.0],
    [106.0, 8.0],
    [122.0, 80.9],
    [123.0, 6.3],
]
Klassen.Molekuelinfo.cNRMdaten = [172.77, 133.83, 130.28, 129.44, 128.49]
Klassen.Molekuelinfo.cdeptdaten = [[133.83, 1],[130.28, 1],[129.44, 1],[128.49, 1]]


Testcase.Case2()
Klassen.Molekuelinfo.VorbereitungNMRdatenfürCalcheuristik()


# Programm startet

if Klassen.Molekuelinfo.cNRMdaten != None or Klassen.Molekuelinfo.cdeptdaten != None:
    Cdept.Cmindestanzahl()
if Klassen.Molekuelinfo.nmrdaten != None and Klassen.Molekuelinfo.d20nmrdaten != None:
    NMR.Spektralanalyse_verschwundeD2Opeaks()


#MS.Summenformelerkennung()


Klassen.Molekuelinfo.isomere = [10,14,4] # Eingriff da es noch nicht wirklich funktioniert
Klassen.Molekuelinfo.Printinfo()




Cdept.CdeptIDanalyse()

Klassen.Molekuelinfo.Printinfo()

isomergruppen = Math.FunktionelleGruppensubstitution()

print("Isomere werden gebildet")

for gruppe in isomergruppen:
    print(gruppe.gruppenkonfiguration)
    individum = GenetischerAlgorythmus.Evolution(gruppe.gruppenkonfiguration,200,100000,240)
    individum.SMilestransformator()
    individum.DarstellungMolekülinSMI(False,True)
    #gruppe.EntwicklungIsomerelist()

print(len(isomergruppen))
print("HAllo")

#isomergruppen[0].EntwicklungIsomerelist()


#Klassen.Molekuelinfo.Printinfo()



#print(len(Formatumwandler.countWasserstoff(isomerearray,[2,0,5],1)))

#print(isomerearray)







