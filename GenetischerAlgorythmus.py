import Klassen
import Testcase

Testcase.Case1()
test = Klassen.individuum([2, 6, 2, 2, 2, 2, 2, 0, 0], 0)

test.SMilestransformator()
#test.DarstellungMolekülinSMI(False, True)

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
