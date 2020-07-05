
#Pattern count
def patternCount(Text, Pattern):
    count=0
    for index in range(0,len(Text)-len(Pattern)+1):
        if Text[index:index+len(Pattern)]==Pattern:
            count=count+1
    return count

def SymbolArray(Symbol,Genome):
    n=len(Genome)
    ExtendedGenome=Genome+Genome[0:n//2]
    count_dict=dict()
    for index in range(n):
        count_dict[index]=patternCount(ExtendedGenome[index:index+(n//2)],Symbol)
    return count_dict

def FasterSymbolArray(Genome, symbol):
    n=len(Genome)
    ExtendedGenome=Genome+Genome[0:n//2]
    count_dict=dict()
    count_dict[0]= patternCount(ExtendedGenome[0:n // 2], symbol)
    for index in range(1,n):
        count_dict[index]=count_dict[index-1]
        if ExtendedGenome[index-1]==symbol:
            count_dict[index]-=1
        if ExtendedGenome[index+(n//2)-1]==symbol:
            count_dict[index]+=1
    return count_dict

def SkewArray(Genome):
    skew_array=[None]*(len(Genome)+1)
    skew_array[0]=0
    for index in range(0,len(Genome)):
        skew_array[index+1]=skew_array[index]
        if Genome[index]=="G":
            skew_array[index+1]+=1
        elif Genome[index]=="C":
            skew_array[index+1]-=1
    return skew_array

def MinimumSkew(Genome):
    skew_array=SkewArray(Genome)
    optimal_min=min(skew_array)
    min_skew=[]
    for index in range(0,len(skew_array)):
        if optimal_min==skew_array[index]:
            min_skew.append(index)
    return min_skew

def HammingDistance(p, q):
    hammingDistance=0
    for index in range(0,len(p)):
        if p[index]!=q[index]:
            hammingDistance+=1
    return hammingDistance


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for index in range(0,len(Text)-len(Pattern)+1):
        subtext=Text[index:index+len(Pattern)]
        hammingDistance=HammingDistance(Pattern,subtext)
        if hammingDistance<=d:
            positions.append(index)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for index in range(0, len(Text) - len(Pattern) + 1):
        subtext = Text[index:index + len(Pattern)]
        hammingDistance = HammingDistance(Pattern, subtext)
        if hammingDistance <= d:
            count+=1
    return count