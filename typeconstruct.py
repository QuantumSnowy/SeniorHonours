from gettext import find
import numpy as np
import csv
import scipy as sy
from operator import itemgetter

def readfile(filename,nbeads):
    data = []
    with open(filename) as f:
        row = f.readlines()
        for line in row:
            data.append(line.split(' '))
    L = float(data[5][1].replace('\n',''))*2
    data = data[9::]
    ends = []
    
    shape = np.array(data).shape
    for i in range(0,shape[0]):
        data[i] = data[i][0:9]
        for j in range(0,len(data[i])):
            data[i][j] = float(data[i][j])
        rem = data[i][0]%nbeads
        if rem==1 or rem==0:
            ends.append(data[i])
    return data, ends, L
        
def read_bonds(filename):
    data = []
    with open(filename) as f:
        row = f.readlines()
        for line in row:
            data.append(line.split(' '))
    data = data[9::]
    bonds = []
    shape = np.array(data).shape
    for i in range(0,shape[0]):
        data[i] = data[i][1:3]
        for j in range(0,len(data[i])):
            data[i][j] = float(data[i][j])
        bonds.append(data[i])
    return bonds

def find_molecule(data,index):
    for i in data:
        if i[0] == index:
            return i[1]
        #else:
        #    print(i[0], index)
        #    print('failed to find molecule')

def find_type(data,index):
    for i in data:
        if i[0] == index:
            return i[2]


def find_molecule1(data,index):
     
     return data[int(index-1)][1]

def find_end(endid,nbeads):
    #for i in list(range(len(ends))):
    #    if ends[i][2] == molnum and ends[i][0] != endind:
    #        otherend = ends[i][0]
    if endid%nbeads == 0:
        otherend = endid + 1 - nbeads
    else:
        otherend = endid + nbeads - 1
    return otherend

#def chainfinder(data):


def sortstructures(filename,nbeads,bonds, n, type):
    #print(nbeads)
    stringy = str(type)
    bigdata2 = readfile("../dump.l30n100m50t" + stringy + "/melt.linear.0",nbeads)[0]
    #bigdata = []
    molly = []
    data = []
    bigdata = sorted(bigdata2, key=itemgetter(0))
    #print("yo" + str(bigdata[0::51]))
    
            
    for i in bigdata:
        if i[2] == n and i[1] not in molly:
            molly.append(i[1])
    #newid = 0
    otherdata = readfile(filename, nbeads)[0]
    for i in otherdata:
        if i[1] in molly:
    #        #newid += 1
    #        #[0] = newid
            data.append(i)
   # count = 1
   # for i in data:
   #     i[0] = count
   #     count += 1
   # print(data[0:10])
    #data = bigdata
   # print(molly)
    #print(data[0])
    bondcopy = list(bonds)    
    #print(len(bondcopy))
    #endbonds = []
    #for i in data:
    #    if i[0]%50 == 0 or i[0]%50 ==1:
    #        endbonds.append(i[0])      
    #print(endbonds)
    #bondcheck = []
    #for i in bondcopy:
    #    for j in endbonds:
    #        if i[0] == j or i[1] == j:
    #            bondcheck.append(i)
    #bondcopy = bondcheck
    #print(data[0:10])
    N = range(len(data))
    rings = []
    chains = []
    taken = []
    p = 0
    
    
    for i in data:
        
        if i[2] == n and find_type(data, i[0]) not in taken:
            #print(taken)
            #print(i, data[i])
            line = []
            endid = i[0]
            #print(endid)
            mol_c = find_molecule(data,endid)
            #print(check)
           
            
            otherend = find_end(endid,nbeads)
            #print(endid, otherend)
            #print(chains)
            if find_type(data, int(otherend)) == n:
                #print("Single Chain")
                #print(data[int(otherend)-1])
                line.append(mol_c)
                taken.append(otherend)
            
            else:
                line.append(mol_c)
                #print("Longer chain")
                while find_type(data, int(otherend)) != n:
                    #print('still going')
                    #print(data[int(otherend)-1])
                    #print(otherend)
                    otherend = find_end(endid,nbeads)
                    #print(find_end(otherend, nbeads))
                    #print(otherend, endid)
                #print("Other end:" +str(otherend))                    
                    #print(otherend)
                    if any(otherend in x for x in bondcopy) == True:
                        search = 0
                        #print(search)
                        while otherend not in bondcopy[search]:
                            search += 1
                            #print(search)   
                        #print(bondcopy[search])
                        #print(otherend)
                        bondcopy[search].remove(otherend)
                        #print(bondcopy[search])
                        #print(bondcopy[search-1:search+2])
                        endid = bondcopy[search][0]
                        #print(search)
                        #print(endid[0])


                        mol_c = find_molecule(data,endid)
                    #print(mol_c)
                        
                        bondcopy.remove(bondcopy[search])
                        #print(bondcopy[search-1:search+2])

                        line.append(mol_c)
                    #print(data[int(endid)][1])
                    else:
                        #print(otherend)
                        endid = otherend
                        #print(endid)
                        #print("End!")
                        break
                #print(line)
            #print(line)
            taken.append(endid)
            chains.append(line)
        else:
            #print(i)
            continue

            
        
        
            
    #print("################################# SEARCHING FOR RINGS #################################")
    #print(taken)
    #print(chains)
    #print(len(bondcopy))
   # print(taken)
    while bondcopy:
        #print(i)
        #print(bondcopy[0][1])
        
        mol_i = find_molecule(data,bondcopy[0][0])
        mol_f = find_molecule(data,bondcopy[0][1])
        if mol_i == mol_f:
            rings.append([mol_i])
            bondcopy.remove(bondcopy[0])
            #print("Single ring")
            #print(rings)
        else:
            cycle = []

            endid = bondcopy[0][1]
            cycle.append(mol_i)
            #print(bondcopy[0])
            #print(mol_i)
            #print(mol_f)
            
            while mol_i != mol_f:
                cycle.append(mol_f)
                #print(endid)
                otherend = find_end(endid,nbeads)
                find = 0
                while otherend not in bondcopy[find]:
                    #print(find)
                    find += 1  
                    #print(len(bondcopy))
                    #print(otherend, bondcopy[find]) 
                    #print(bondcopy[search])
                bondcopy[find].remove(otherend)    
                endid = bondcopy[find][0]
                
                #print(endid)
                bondcopy.remove(bondcopy[find])
                mol_f = find_molecule(data,endid)
            bondcopy.remove(bondcopy[0])
            rings.append(cycle)
    #print(rings)


    return chains, rings



def lengths(types, timestep):
    
    filename = timestep
    bigmean = []
    appendage = str(types)
    bigchain = []
    bigring = []
    for i in range(2,types+1):

        ends = readfile("../dump.l30n100m50t" + appendage + "/melt.linear." + filename,50)[0]            ####For my own usage to read files
        #print(ends)
        num = str(i)
        bonds = read_bonds("../dump.l30n100m50t" + appendage + "/bonds" + num + "." + filename)
        #print(bonds)
        chains, rings = sortstructures("../dump.l30n100m50t" + appendage +  "/melt.linear." + filename,50,bonds, i, types)
        bigchain += chains
        bigring += rings
        #print(rings)
        #print(chains)
#Kymo_file_writer(filename, "../dump.l30n100m50.type2/melt.linear.", "../dump.l30n100m50.type2/bonds2.", 50, 30, 2)
        tick = []    
        for i in rings:
           tick.append(50*(len(i)))
        for i in chains:
           tick.append(50*(len(i)))
        bigmean.append((np.mean(np.array(tick))))
    final = np.mean(np.array(bigmean))
    var = np.std(np.array(bigmean))/np.sqrt(np.size(np.array(bigmean)))
    #print(final, var)
    return bigchain, bigring, final, var


woop = np.array([2,3,5,6,11,21,26])
timmy = np.arange(1000000, 6000000, step=1000000, dtype=int)
listy = []
for i in timmy:
    lengthdata = open(str(i) + "data.dat",  'w')
    lengthdata.write("Families" + " " + "Length"+ " " + "Error" + "\n") 
    chaindata = open(str(i) + "chainringdata.dat",  'w')
    chaindata.write("Families" + " " + "Chain"+  " " + "Ring" + "\n")
    step = str(i)
#    
    for j in woop:
        chains, rings, length, var = lengths(j, step)
        print(len(chains), len(rings))
        #print(chains)
        lengthdata.write(str(j-1) + " " + str(length)+ " " + str(var) + "\n")
        chaindata.write(str(j-1) + " " + str(len(chains))+ " " +  str(len(rings)) + "\n")
#    print(j)
#lengths(26, "4000000")
