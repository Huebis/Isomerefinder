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

    individuen.append(Klassen.individuum(gruppenkonfiguration,0))
    individuen[0].CalcHeuristik()

    for a in range(grössePopulation-1):
        individuen.append(Klassen.individuum(gruppenkonfiguration,100))
        individuen[a+1].CalcHeuristik()



    while anzahlgenerationen != 0:
        anzahlgenerationen -= 1

        #vier Eltern werden aus gewählt und je zwei "kämpfen" gegeneinander
        # Es ist beabsichtigt, dass es auch überschneidungen geben kann
        positionVater1 = random.randint(0,grössePopulation + populationszuwachs-1)
        positionVater2 = random.randint(0, grössePopulation + populationszuwachs - 1)
        positionMutter1 = random.randint(0, grössePopulation + populationszuwachs- 1)
        positionMutter2 = random.randint(0, grössePopulation + populationszuwachs- 1)

        #vier Eltern werden aus gewählt und je zwei "kämpfen" gegeneinander

        if individuen[positionVater1].heuristikwert > individuen[positionVater2].heuristikwert:
            if individuen[positionMutter1].heuristikwert > individuen[positionMutter2].heuristikwert:
                individuen.append(Klassen.individuum(None,None,individuen[positionVater2].molekularstruktur,individuen[positionMutter2].molekularstruktur,individuen[positionMutter1].elemente))
            else:
                individuen.append(Klassen.individuum(None, None, individuen[positionVater2].molekularstruktur,individuen[positionMutter1].molekularstruktur,individuen[positionMutter1].elemente))
        else:
            if individuen[positionMutter1].heuristikwert > individuen[positionMutter2].heuristikwert:
                individuen.append(Klassen.individuum(None, None, individuen[positionVater1].molekularstruktur,individuen[positionMutter2].molekularstruktur,individuen[positionMutter1].elemente))
            else:
                individuen.append(Klassen.individuum(None, None, individuen[positionVater1].molekularstruktur,individuen[positionMutter1].molekularstruktur,individuen[positionMutter1].elemente))


        #überprüfen ob es wirklich ein Kind gegeben hat

        if individuen[-1].molekularstruktur == None:
            individuen.pop(-1)
        else:
            individuen[-1].CalcHeuristik()
            populationszuwachs += 1


        if anzahlpopulationenzuwachchs_bis_schlechte_sterben == populationszuwachs:
            individuen.sort(key=lambda a: a.heuristikwert)
            while populationszuwachs != 0:
                individuen.pop(-1)
                populationszuwachs -= 1

    individuen.sort(key=lambda a: a.heuristikwert)
    return individuen[0]









