from itertools import repeat

import Testcase
import random
import Klassen
import Math
import itertools

zusatz = list(itertools.product([0,1], repeat=2))
print(list())
Math.allemöglichenKombinationenundjedermussmindestenseinalvorkommen(8,5)
"""
hallo = [[2.3, 1, [1, 1]], [2.3, 1, [0, 1], [4, 1]], [2.3, 1, [3, 1]], [2.3, 1, [2, 1], [7, 1]], [2.3, 1, [1, 1], [5, 1]], [2.3, 1, [4, 1], [6, 1]], [2, 2, [5, 1]], [2, 2, [3, 1]]]

test = Klassen.individuum([2, 6, 2, 2, 2, 2, 2, 0, 0], 0)
test.molekularstruktur = [[8, [1, 2]], [1, [2, 2], [1, 3]], [1, [2, 1], [1, 0]], [1, [1, 1], [1, 4], [1, 5]], [1, [2, 6], [1, 3]], [1, [2, 7], [1, 3]], [1, [2, 4], [1, 8]], [2, [2, 5]], [2, [1, 6], [1, 9]], [8, [1, 8]]]
test.SMilestransformator()
test.DarstellungMolekülinSMI(False, True)
"""
hallo = "C((C(C(=C(C(=C))))(C(=C(C(C(=O))))))=O)"

Testcase.Case2()
print(Klassen.Molekuelinfo.gruppierted20nmrdaten)



Testcase.Case2()
test = Klassen.individuum([2, 6, 2, 0, 4, 0, 0, 0, 0], 0)

test.SMilestransformator()
test.DarstellungMolekülinSMI(False, False)
test.CalcHeuristik()


print("Heuristkwert: " + str(test.heuristikwert))



test.SMilestransformator()
#test.DarstellungMolekülinSMI(False, True)


#child.SMilestransformator()
#child.DarstellungMolekülinSMI(False,True)

print("finish")



"""               #nun werden die NMR daten Analysiert
                nmrwerte = copy.deepcopy(Molekuelinfo.gruppierted20nmrdaten)

                #jeder approximationswert wird mit jedem nmrwert ausprobiert und verglichen
                def nmrwertemitapproximationvergleich(nmrwert,aproximierterwert):
                    unterschied = (nmrwert[0]-aproximierterwert[0])**2
                    unterschied += abs(len(nmrwert) - len(aproximierterwert))*2
                    #Idee ob es ohne besser funktionert
                    if len(nmrwert) > len(aproximierterwert):
                        unterschied -= len(aproximierterwert)*4
                    else:
                        unterschied -= len(nmrwert) * 4

                    return unterschied



                vergleiche = []

                for a in range(len(approximation)):
                    vergleiche.append([])
                    for b in range(len(nmrwerte)):
                        vergleiche[a].append(nmrwertemitapproximationvergleich(nmrwerte[b],approximation[a]))



                transformation = [None for a in range(len(approximation))] # jeder approximationswert bekommt ein NMRwert zugewisen

                #Es darf auch mehrfachbelegungen geben (aber in einer genau bestimmen Anzahl)
                anzahlmehrfachbelegungen =  len(approximation) - len(nmrwerte)

                for temp in range(len(approximation)):
                    # da der Wert möglichst klein sein soll , wird jetzt das Atom ausgewält, bei welchem das beste, das schlechteste ist
                    besterwertjederspalte = [None for a in range(len(approximation))]
                    for a in range(len(approximation)):
                        spalenwerte = []
                        for b in range(len(nmrwerte)):
                            if vergleiche[a][b] != None:
                                spalenwerte.append(vergleiche[a][b])
                        if spalenwerte != []:
                            besterwertjederspalte[a] = min(spalenwerte)

                    print("besterwertjederspalte")
                    print(besterwertjederspalte)
                    print("vergleiche")
                    print(vergleiche)
                    print("Transformation")
                    print(transformation)
                    print("Übrige mehrfachbelegungen")
                    print(anzahlmehrfachbelegungen)

                    schlechtesterWert = max([wert for wert in besterwertjederspalte if wert is not None])

                    breakbool = False
                    for a,wert in enumerate(besterwertjederspalte):

                        if wert == schlechtesterWert:
                            for b in range(len(nmrwerte)):
                                if vergleiche[a][b] == wert:
                                    if b in transformation:
                                        anzahlmehrfachbelegungen -= 1

                                    transformation[a] = b
                                    for c in range(len(nmrwerte)):
                                        vergleiche[a][c] = None

                                    if anzahlmehrfachbelegungen <= 0:
                                        for spalte in vergleiche:
                                            spalte[b] = None

                                    breakbool = True
                                    break
                        if breakbool:
                            break

                print("nmrwerte")
                print(nmrwerte)
                print("approximationen")
                print(approximation)

                # Noch nicht final aber jetzt kann mal alles getestet werden
                summe = 0
                for a in range(len(approximation)):
                    summe += nmrwertemitapproximationvergleich(nmrwerte[transformation[a]], approximation[a])"""






"""


#Versuch um zu erkennen ob Cyclogrösse Aromate wirklich erkennt (er erkennts)

import Klassen




test = Klassen.individuum([0, 6, 0, 0, 0, 0, 0, 0, 0], 0)


test.SMilestransformator()
test.DarstellungMolekülinSMI(True, True)

print(test.Cyclogrösse(0))

"""


""" Abgespeichert, falls nochmals gebraucht 

 def BewertungvonCH2undCH1gruppenn():
                approximation = [None for a in range(len(self.molekularstruktur))]

                for position,atom in enumerate(self.molekularstruktur):
                    if atom[0] == self.elementgruppengrenzen[3]:  # ist ein CH3
                        approximation[position] = [1]
                        #spezifikationen
                        if atom[1][1] == self.elementgruppengrenzen[2] + 2:
                            approximation[position][0] += 2

                    elif atom[0] == self.elementgruppengrenzen[2]:  # ist ein CH2
                        if len(atom) == 2: # es ist ein Olephin
                            approximation[position] = [5.5]
                        else:
                            approximation[position] = [2]


                    elif atom[0] == self.elementgruppengrenzen[1]:  # ist ein CH
                        if self.Iszyklisch(position):
                            if self.Cyclogrösse(position) == 0: # Spezialfall wenn es ein Aromat ist, ansonsten hat es ca den wert 2.3
                                approximation[position] = [7.3]

                            else:
                                approximation[position] = [2.3]

                        else:
                            approximation[position] = [2.3]


                    elif atom[0] == self.elementgruppengrenzen[3] + 2:  # ist ein C=O(H)
                        approximation[position] = [9.5]


                for a in range(len(approximation)):
                    if approximation[a] != None:
                        kopplungen = self.PositionendernachbarnwelcheeineWasserstoffkopplungeingehen(a)
                        if kopplungen != None:
                            for kopplung in kopplungen:
                                approximation[a].append(kopplung)

                #Nun hat es in der Liste noch gewisse Stellen mit None. Diese Müssen jetzt gelöscht werden und alle Positionen müssen angepasst werden

                #zuerst müssen noch lehre [] gelöscht werden, welche wege self.PositionendernachbarnwelcheeineWasserstoffkopplungeingehen




                for a in range(len(approximation)-1,-1,-1):
                    if approximation[a] == None:
                        for b in range(len(approximation)):
                            if approximation[b] != None:
                                for c in range(1,len(approximation[b])):
                                    if approximation[b][c][0] > a:
                                        approximation[b][c][0] -= 1
                        approximation.pop(a)
             
                print("approximation")
                print(approximation)
                print(self.molekularstruktur)
                

                #nun werden die NMR daten Analysiert
                nmrwerte = copy.deepcopy(Molekuelinfo.gruppierted20nmrdaten)

                #jeder approximationswert wird mit jedem nmrwert ausprobiert und verglichen
                def nmrwertemitapproximationvergleich(nmrwert,aproximierterwert):
                    unterschied = (nmrwert[0]-aproximierterwert[0])**2
                    unterschied += abs(len(nmrwert) - len(aproximierterwert))*2
                     #Idee ob es ohne besser funktionert
                    if len(nmrwert) > len(aproximierterwert):
                        unterschied -= len(aproximierterwert)*4
                    else:
                        unterschied -= len(nmrwert) * 4

                    return unterschied



                vergleiche = []

                for a in range(len(approximation)):
                    vergleiche.append([])
                    for b in range(len(nmrwerte)):
                        vergleiche[a].append(nmrwertemitapproximationvergleich(nmrwerte[b],approximation[a]))



                transformation = [None for a in range(len(approximation))] # jeder approximationswert bekommt ein NMRwert zugewisen

                #Es darf auch mehrfachbelegungen geben (aber in einer genau bestimmen Anzahl)
                anzahlmehrfachbelegungen =  len(approximation) - len(nmrwerte)

                for temp in range(len(approximation)):
                    # da der Wert möglichst klein sein soll , wird jetzt das Atom ausgewält, bei welchem das beste, das schlechteste ist
                    besterwertjederspalte = [None for a in range(len(approximation))]
                    for a in range(len(approximation)):
                        spalenwerte = []
                        for b in range(len(nmrwerte)):
                            if vergleiche[a][b] != None:
                                spalenwerte.append(vergleiche[a][b])
                        if spalenwerte != []:
                            besterwertjederspalte[a] = min(spalenwerte)
                          
                    print("besterwertjederspalte")
                    print(besterwertjederspalte)
                    print("vergleiche")
                    print(vergleiche)
                    print("Transformation")
                    print(transformation)
                    print("Übrige mehrfachbelegungen")
                    print(anzahlmehrfachbelegungen)
                  
                    schlechtesterWert = max([wert for wert in besterwertjederspalte if wert is not None])

                    breakbool = False
                    for a,wert in enumerate(besterwertjederspalte):

                        if wert == schlechtesterWert:
                            for b in range(len(nmrwerte)):
                                if vergleiche[a][b] == wert:
                                    if b in transformation:
                                        anzahlmehrfachbelegungen -= 1

                                    transformation[a] = b
                                    for c in range(len(nmrwerte)):
                                        vergleiche[a][c] = None

                                    if anzahlmehrfachbelegungen <= 0:
                                        for spalte in vergleiche:
                                            spalte[b] = None

                                    breakbool = True
                                    break
                        if breakbool:
                            break

                print("nmrwerte")
                print(nmrwerte)
                print("approximationen")
                print(approximation)
                # Noch nicht final aber jetzt kann mal alles getestet werden
                summe = 0
                for a in range(len(approximation)):
                    summe += nmrwertemitapproximationvergleich(nmrwerte[transformation[a]], approximation[a])

                return summe
            print("Methylgruppen")
"""