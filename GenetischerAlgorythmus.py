import copy
import time
import Klassen
import Testcase
import random
"""
Testcase.Case1()
test = Klassen.individuum([2, 6, 2, 2, 2, 2, 2, 0, 0], 0)

test.SMilestransformator()
test.DarstellungMolekülinSMI(False, False)
test.CalcHeuristik()


for a in range(10):
    test.Muation()
    if not test.isligit():
        print("FEHLER")
        print(test.molekularstruktur)
        raise ValueError("Ungültige Molekularstruktur")


test2 = Klassen.individuum([2, 6, 2, 2, 2, 2, 2, 0, 0], 0)


test.SMilestransformator()
#test.DarstellungMolekülinSMI(False, True)

test2.SMilestransformator()
#test2.DarstellungMolekülinSMI(False, True)

child = Klassen.individuum(None,0,test.molekularstruktur,test2.molekularstruktur,test.elemente)

#child.SMilestransformator()
#child.DarstellungMolekülinSMI(False,True)

print("finish")
"""

def Evolution(gruppenkonfiguration,grössePopulation, anzahlgenerationen, anzahlpopulationenzuwachchs_bis_schlechte_sterben):
    # Die grösse der Population muss mindestens 4 sein, aber mindestens 50 sind empfohlen
    individuen = []
    populationszuwachs = 0

    print("Gruppenkonfiguration")
    print(gruppenkonfiguration)
    individuen.append(Klassen.individuum(gruppenkonfiguration,0))
    individuen[0].CalcHeuristik()

    for a in range(grössePopulation-1):
        individuen.append(Klassen.individuum(gruppenkonfiguration,100))
        individuen[a+1].CalcHeuristik()



    while anzahlgenerationen != 0:
        anzahlgenerationen -= 1
        if anzahlgenerationen % 16 == 0:
            print("Noch " + str(anzahlgenerationen) + " Generationen")

        #vier Eltern werden aus gewählt und je zwei "kämpfen" gegeneinander
        # Es ist beabsichtigt, dass es auch überschneidungen geben kann
        positionVater1 = random.randint(0,grössePopulation + populationszuwachs-1)
        positionVater2 = random.randint(0, grössePopulation + populationszuwachs - 1)
        positionMutter1 = random.randint(0, grössePopulation + populationszuwachs- 1)
        positionMutter2 = random.randint(0, grössePopulation + populationszuwachs- 1)

        #vier Eltern werden aus gewählt und je zwei "kämpfen" gegeneinander

        #start = time.time()

        if individuen[positionVater1].heuristikwert > individuen[positionVater2].heuristikwert:
            if individuen[positionMutter1].heuristikwert > individuen[positionMutter2].heuristikwert:
                individuen.append(Klassen.individuum(None,None,copy.deepcopy(individuen[positionVater2].molekularstruktur),copy.deepcopy(individuen[positionMutter2].molekularstruktur),copy.deepcopy(individuen[positionMutter1].elemente)))
            else:
                individuen.append(Klassen.individuum(None, None, copy.deepcopy(individuen[positionVater2].molekularstruktur),copy.deepcopy(individuen[positionMutter1].molekularstruktur),copy.deepcopy(individuen[positionMutter1].elemente)))
        else:
            if individuen[positionMutter1].heuristikwert > individuen[positionMutter2].heuristikwert:
                individuen.append(Klassen.individuum(None, None, copy.deepcopy(individuen[positionVater1].molekularstruktur),copy.deepcopy(individuen[positionMutter2].molekularstruktur),copy.deepcopy(individuen[positionMutter1].elemente)))
            else:
                individuen.append(Klassen.individuum(None, None, copy.deepcopy(individuen[positionVater1].molekularstruktur),copy.deepcopy(individuen[positionMutter1].molekularstruktur),copy.deepcopy(individuen[positionMutter1].elemente)))

        #ende = time.time()

        #print("Zeit für das genrieren eines Kindes: " + str(ende - start))
        #überprüfen ob es wirklich ein Kind gegeben hat


        if individuen[-1].molekularstruktur == None:
            individuen.pop(-1)
        else:
            #start = time.time()
            individuen[-1].CalcHeuristik()
            populationszuwachs += 1
            #print(individuen[-1].heuristikwert)

            #ende = time.time()
            #print("Zeit zum berechnen der Heuristik " + str(ende - start))

        #eine gewisse Anzahl an Individuen wird dubliziert und durch Mutationen verändert und der Population hinzugefügt.
        #start = time.time()


        sum = 0
        for a in range(2):
            pos_ausgewähltes_Individuum = random.randint(0,len(individuen)-1)
            #print(individuen[pos_ausgewähltes_Individuum].molekularstruktur)
            individuen.append(Klassen.individuum(None,2,None,None,copy.deepcopy(individuen[pos_ausgewähltes_Individuum].elemente),copy.deepcopy(individuen[pos_ausgewähltes_Individuum].molekularstruktur)))
            individuen[-1].CalcHeuristik()
            populationszuwachs += 1
        #ende = time.time()

        #print("Zeit zur Mutationen" + str(ende - start))

        #start = time.time()
        if anzahlpopulationenzuwachchs_bis_schlechte_sterben == populationszuwachs:
            individuen.sort(key=lambda a: a.heuristikwert)
            while populationszuwachs != 0:
                individuen.pop(-1)
                populationszuwachs -= 1

        #individuen.sort(key=lambda a: a.heuristikwert)
        # ende = time.time()
        # print("Zeit zum dezimieren der Population" + str(ende - start))
        #print("Anzahl Moleküle :" +str(len(individuen)))
    return individuen[0]









