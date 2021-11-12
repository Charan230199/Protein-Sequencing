"""
Protein Sequencing Project
Name: G.V.S Sai Charan
Roll Number: 2021501008
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    read = open(filename)
    text = read.read().splitlines()
    str ="".join(text)
    return str


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    d=dna[startIndex:]
    list=[]
    string=""
    replaced=["UAG", "UAA","UGA"]
    for i in range(len(d)):
        if len(string)!=3:
            string+=d[i]
        if len(string)==3:
            x= string.replace("T","U")
            if x in replaced:
                list.append(x)
                return list
            else:
                list.append(x)
                string=""
    return list


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename)
    read = json.load(f)
    Dict={}
    for x,y in read.items():
        for i in y:
            Dict[i.replace('T','U')]=x
    return Dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    lst=[]
    if codons[0] =="AUG":
        lst.append("Start")
    for i in range(1,len(codons)):
        if codons[i] in codonD.keys():
            lst.append(codonD[codons[i]])
    return lst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    file = readFile(dnaFilename) 
    Dict = makeCodonDictionary(codonFilename) 
    i=0
    count=0
    temp=[]
    while i < len(file):
        if file[i:i+3] == "ATG":
            dna= dnaToRna(file,i) 
            protien = generateProtein(dna, Dict) 
            temp.append(protien)
            i = i+3*len(dna)
        else:
            i+=1
            count+=1
    return temp


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    temp=[]
    for a in proteinList1:
        for b in proteinList2:
            if a==b and a not in temp:
                temp.append(a)

    return temp


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    lst=[]
    for i in proteinList:
        for word in i:
            lst.append(word)
    return lst


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    Dict={}
    for x in aaList:
        if x not in Dict:
            Dict[x]=1
        else:
            Dict[x]+=1
    return Dict



'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    c1,c2 = combineProteins(proteinList1), combineProteins(proteinList2) 
    dict1,dict2 = aminoAcidDictionary(c1),aminoAcidDictionary(c2) 
    temp,result=[],[]            
    fd1,fd2={},{}  
    for i in dict1:
        fd1[i] = dict1[i]/len(c1)
        if i not in temp and i !="Start" and i!="Stop":
            temp.append(i)
    for j in dict2:
        fd2[j] = dict2[j]/len(c2)
        if j not in temp and j !="Start" and j!="Stop":
            temp.append(j)
    for a in temp:
        freq1,freq2=0,0
        if a in fd1:
            freq1= fd1[a]
        if a in fd2:
            freq2= fd2[a]
        difference = freq2-freq1
        if difference < -cutoff or difference > cutoff  :
            result.append([a , freq1, freq2])
    return result


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("Print the Commnalities!")
    for i in sorted(commonalities):
        commonProteins = []
        let = i[1:len(i)-1] 
        count=0
        for j in let:
            commonProteins+=j   
            count+=1            
            if count !=len(let): 
                commonProteins+="-" 
        print(commonProteins)
    print("Print DNA sequences!")
    for item in differences:
        print(item[0],":",round(item[1]*100,2),"% in seq1,",round(item[2]*100,2),"% in seq2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testSynthesizeProteins()

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
