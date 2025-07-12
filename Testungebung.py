
def Anzahloffeneverbindungen(Atomgruppe):



    anzahlverbrauchteverbindungen = 0
    for a in range(1 ,len(Atomgruppe) ,1):
        anzahlverbrauchteverbindungen += At[a][0]
        print([a][0])

    return anzahlverbrauchteverbindungen


test = [[0, [3, 1], [1, 2]], [1, [3, 0], [1, 3]], [1, [1, 0], [1, 4]], [2, [1, 1]], [2, [1, 2]], [2, [1, 0]], [2], [2], [3], [3]]

print(Anzahloffeneverbindungen(test[0]))