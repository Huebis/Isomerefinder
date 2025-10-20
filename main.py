import Cdept
import Formatumwandler
import MS
import Klassen
import NMR
import Math
import Testcase
import GenetischerAlgorythmus


#Hier kann ein Testcase angegeben werden, dieser wird aus dem File Testcase.py genommen.
Testcase.Case6()



# Programm startet

if Klassen.Molekuelinfo.cNRMdaten != None or Klassen.Molekuelinfo.cdeptdaten != None:
    Cdept.Cmindestanzahl()
if Klassen.Molekuelinfo.nmrdaten != None and Klassen.Molekuelinfo.d20nmrdaten != None:
    NMR.Spektralanalyse_verschwundeD2Opeaks()


#MS.Summenformelerkennung()


Klassen.Molekuelinfo.Printinfo()




Cdept.CdeptundSymetrieIDanalyse()

Klassen.Molekuelinfo.Printinfo()

Gruppenkonfigurationen = Math.FunktionelleGruppensubstitution()

print("Isomere werden gebildet")

print("Anzahl Gruppenkonfigurationen: " + str(len(Gruppenkonfigurationen)))
print("Alle Gruppenkonfigurationen")
for gruppe in Gruppenkonfigurationen:
    print(gruppe.gruppenkonfiguration)


gewinner = []
for gruppe in Gruppenkonfigurationen:
    print(gruppe.gruppenkonfiguration)
    gewinner.append(GenetischerAlgorythmus.Evolution(gruppe.gruppenkonfiguration,60,1000,60))
    print(gewinner[-1].molekularstruktur)
    gewinner[-1].SMilestransformator()
    #gewinner[-1].DarstellungMolekülinSMI(True,True) #Hier kann ausgewählt werden ob jeder Gewinner aus dem Genetischen Algorithmus ausgegeben werden soll (als Bild)
    print("Gewinnerheurist: "+ str(gewinner[-1].CalcHeuristik()))
    #gruppe.EntwicklungIsomerelist()


besteheuristik = gewinner[0].heuristikwert
position = 0
for pos,gewinner_Individiuum in enumerate(gewinner):
    print("Gewinnerheurist: " + str(gewinner_Individiuum.heuristikwert))
    if gewinner_Individiuum.heuristikwert < besteheuristik:
        besteheuristik = gewinner_Individiuum.heuristikwert
        position = pos


gewinner[position].SMilestransformator()
gewinner[position].DarstellungMolekülinSMI(True,True)




for gruppe in Gruppenkonfigurationen:
    print(gruppe.gruppenkonfiguration)
print(len(Gruppenkonfigurationen))
print("HAllo")

#isomergruppen[0].EntwicklungIsomerelist()


#Klassen.Molekuelinfo.Printinfo()



#print(len(Formatumwandler.countWasserstoff(isomerearray,[2,0,5],1)))

#print(isomerearray)







