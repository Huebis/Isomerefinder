import Cdept
import Formatumwandler
import MS
import Klassen
import NMR
import Math
import Testcase
import GenetischerAlgorythmus


#Hier kann ein Testcase angegeben werden, dieser wird aus dem File Testcase.py genommen.
#Testcase.Case1() weist noch manchmal einen Fehler auf, bis jetzt habe ich noch nicht herausgefunden warum. Bei allen anderen Cases sind mir keine Bugs bekannt
#Falls ein Bug bzw. eine Fehlermeldung bei einem Testcase auftaucht, kann man mir gerne Schreiben (oder Issue aufmachen).
Testcase.Case4()



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
    #Erste Parameter gibt an, wie gross die Startpopulation ist
    #Zweiter Parameter gibt an, wie viele Generationen durchgemacht werden.
    #Dritter Parameter gibt an, um wie viele Individuen die Population wachsen kann, bis sie wieder zur grösse der Startpopulation dezimiert wird.
    gewinner.append(GenetischerAlgorythmus.Evolution(gruppe.gruppenkonfiguration,60,10000,60))
    gewinner[-1].SMilestransformator()
    #gewinner[-1].DarstellungMolekülinSMI(True,True) #Hier kann ausgewählt werden ob jeder Gewinner aus dem Genetischen Algorithmus ausgegeben werden soll (als Bild)
    print("Gewinnerheurist: "+ str(gewinner[-1].CalcHeuristik()))






besteheuristik = gewinner[0].heuristikwert
position = 0
for pos,gewinner_Individiuum in enumerate(gewinner):
    print("Gewinnerheurist: " + str(gewinner_Individiuum.heuristikwert))
    if gewinner_Individiuum.heuristikwert < besteheuristik:
        besteheuristik = gewinner_Individiuum.heuristikwert
        position = pos

print("BESTES INDIVDIUM:")
print("Heuristik: " + str(gewinner[position].heuristikwert))
print("Molekülstruktur: " + str(gewinner[position].molekularstruktur))
gewinner[position].SMilestransformator()
#Erster Parameter gibt an, ob die Wasserstoffatome auch gezeichnet werden sollen.
#Zweiter Parameter gibt an, ob das Bild ausgegeben werden soll (sollte auf dem Bildschirm erscheinen).
#Dritter Parameter gibt an, ob der SMILES-String in der Console ausgegeben werden soll.
gewinner[position].DarstellungMolekülinSMI(True,True,True)



print("Nochmals alle Durchgerechneten Gruppenkonfigurationen")
for gruppe in Gruppenkonfigurationen:
    print(gruppe.gruppenkonfiguration)
print("-FINISH-")

#isomergruppen[0].EntwicklungIsomerelist()


#Klassen.Molekuelinfo.Printinfo()



#print(len(Formatumwandler.countWasserstoff(isomerearray,[2,0,5],1)))

#print(isomerearray)







