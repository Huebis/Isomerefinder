import random


hallo = []
for a in range(10):
    test = [[random.randint(0,10)],[random.randint(0,10)]]
    hallo.append([test,1])

print(hallo)

daten = [
    ["Apfel", 3],
    ["Banane", 1],
    ["Kirsche", 5]
]


daten.sort(key=lambda x: x[1])


print(daten)



test = [[1,2],[3,4]]
for a in test:
    a[0] += 1

print(test)

hallo = [[0]]

print(hallo[0][0])
hallo[0].append([3,2])
print(hallo[0])


for a in range(2,10,1):
    print(a)

