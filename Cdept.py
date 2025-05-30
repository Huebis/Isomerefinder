
def cdeptIDanalyse(cdeptdaten, lenofcNMR):
    cdeptID = [0,0,0]
    ceptID[0] = lenofcNMR - len(cdeptdaten)
    for a in cdeptdaten:
        if a[1] < 0:
            cdeptID[1] += 1
        else:
            cdeptID[2] += 1
    return cdeptID
