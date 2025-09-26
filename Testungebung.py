import Testcase
import random
import Klassen


Testcase.Case1()
test = Klassen.individuum([2, 6, 2, 2, 2, 2, 2, 0, 0], 0)

test.SMilestransformator()
test.DarstellungMolekülinSMI(False, True)
test.CalcHeuristik()


print("Heuristkwert" + str(test.heuristikwert))



test.SMilestransformator()
#test.DarstellungMolekülinSMI(False, True)


#child.SMilestransformator()
#child.DarstellungMolekülinSMI(False,True)

print("finish")













#Versuch um zu erkennen ob Cyclogrösse Aromate wirklich erkennt (er erkennts)
"""
import Klassen




test = Klassen.individuum([0, 6, 0, 0, 0, 0, 0, 0, 0], 0)


test.SMilestransformator()
test.DarstellungMolekülinSMI(True, True)

print(test.Cyclogrösse(0))

"""