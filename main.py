import Cdept
import Formatumwandler
import MS
import Klassen
import NMR
import Math
import Testcase
import GenetischerAlgorythmus


Testcase.Case2()



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

print(Gruppenkonfigurationen)
gewinner = []
for gruppe in Gruppenkonfigurationen:
    print(gruppe.gruppenkonfiguration)
    gewinner_individum = GenetischerAlgorythmus.Evolution(gruppe.gruppenkonfiguration,60,10000,60)
    print(gewinner_individum.molekularstruktur)
    gewinner_individum.SMilestransformator()
    gewinner_individum.DarstellungMolekülinSMI(False,True)
    print(gewinner_individum.CalcHeuristik())
    #gruppe.EntwicklungIsomerelist()


for gruppe in Gruppenkonfigurationen:
    print(gruppe.gruppenkonfiguration)
print(len(Gruppenkonfigurationen))
print("HAllo")

#isomergruppen[0].EntwicklungIsomerelist()


#Klassen.Molekuelinfo.Printinfo()



#print(len(Formatumwandler.countWasserstoff(isomerearray,[2,0,5],1)))

#print(isomerearray)







