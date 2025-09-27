import Cdept
import Formatumwandler
import MS
import Klassen
import NMR
import Math
import Testcase
import GenetischerAlgorythmus


Testcase.Case3()



# Programm startet

if Klassen.Molekuelinfo.cNRMdaten != None or Klassen.Molekuelinfo.cdeptdaten != None:
    Cdept.Cmindestanzahl()
if Klassen.Molekuelinfo.nmrdaten != None and Klassen.Molekuelinfo.d20nmrdaten != None:
    NMR.Spektralanalyse_verschwundeD2Opeaks()


#MS.Summenformelerkennung()


Klassen.Molekuelinfo.Printinfo()




Cdept.CdeptIDanalyse()

Klassen.Molekuelinfo.Printinfo()

isomergruppen = Math.FunktionelleGruppensubstitution()

print("Isomere werden gebildet")

print(isomergruppen)

for gruppe in isomergruppen:
    print(gruppe.gruppenkonfiguration)
    individum = GenetischerAlgorythmus.Evolution(gruppe.gruppenkonfiguration,200,100,240)
    individum.SMilestransformator()
    individum.DarstellungMolekülinSMI(False,True)
    print(individum.CalcHeuristik())
    #gruppe.EntwicklungIsomerelist()

print(len(isomergruppen))
print("HAllo")

#isomergruppen[0].EntwicklungIsomerelist()


#Klassen.Molekuelinfo.Printinfo()



#print(len(Formatumwandler.countWasserstoff(isomerearray,[2,0,5],1)))

#print(isomerearray)







