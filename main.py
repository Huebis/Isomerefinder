import Formatumwandler

isomerearray = Formatumwandler.Readfile()
isomerearray = Formatumwandler.AddWasserstoff(isomerearray)
print(isomerearray)