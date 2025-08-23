import Klassen
import Testcase

Testcase.Case1()
test = Klassen.individuum([2, 6, 2, 2, 2, 2, 2, 0, 0], 0)

test.SMilestransformator()
test.DarstellungMolekülinSMI(False, True)

for a in range(10000):
    test.Muation()
    test.CalcHeuristik()
    if not test.isligit():
        print("FEHLER")
        print(test.molekularstruktur)
        raise ValueError("Ungültige Molekularstruktur")

test.SMilestransformator()
test.DarstellungMolekülinSMI(False, True)
# test.testfunktion()


# for a in range(100000  ):
# test = individuum([random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20), random.randint(1,20)], 0)

print("finish")
# testfehlschlänge

# [5,0,0,0,0,0,0,0,0]