import random


# counts the occurances of each "ACTG" in the motfis array
# returns a dictionary with counts for each base in the j-th column of the motifs array/matrix
def count(Motifs):
    count = initializeCountArray(Motifs)
    t = len(Motifs)
    for i in range(t):
        k = len(Motifs[i])
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


# helper method
def initializeCountArray(Motifs):
    count = {}
    length = len(Motifs)
    for index in range(0, length):
        k = len(Motifs[index])
        for symbol in "ACGT":
            count[symbol] = []
            for j in range(0, k):
                count[symbol].append(1)
    return count


# retunrs a profile with the frequency of occurance for each base in the j-th column of the motifs matrix
def profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = count(Motifs)
    for symbol, countArray in count.items():
        profile[symbol] = [x / t for x in countArray]
    return profile


# finds the consensus string from the count matrix for set of motifs
def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    count = CountWithPseudocounts(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    k = len(Motifs[0])
    for j in range(0, k):
        for symbol in "ACGT":
            if symbol != consensus[j]:
                score += count[symbol][j]
    return score


# calculates the probability of Text using the Profile matrix of probabilities
def Pr(Text, Profile):
    pr = 1
    k = len(Text)
    for index in range(0, k):
        for symbol in "ACGT":
            if Text[index] == symbol:
                pr = pr * Profile[symbol][index]
    return pr


def ProfileMostProbableKmer(text, k, profile):
    highestProbability = -1
    mostProbablePattern = ""
    for index in range(0, len(text) - k + 1):
        pattern = text[index:index + k]
        prob = Pr(pattern, profile)
        if prob > highestProbability:
            highestProbability = prob
            mostProbablePattern = pattern
    return mostProbablePattern


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    # get the first k-mers from each DNA to form an intitial array of best motifs
    for index in range(0, t):
        BestMotifs.append(Dna[index][0:k])
    # going trough all the possible k-mers in DNA[O] and forming motifs matrix from the rest dna strings
    # finding the best motifs by caluclating the best score (which is the lowest score computed for the motifs)
    n = len(Dna[0])
    for j in range(0, n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][j:j + k])
        for l in range(1, t):
            profile = profile(Motifs[0:l])
            Motifs.append(ProfileMostProbableKmer(Dna[l], k, profile))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def CountWithPseudocounts(Motifs):
    count = initializeCountArray(Motifs)
    t = len(Motifs)
    for i in range(t):
        k = len(Motifs[i])
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ProfileWithPseudocounts(Motifs, l):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)
    for symbol, countArray in count.items():
        profile[symbol] = [x / (t + l) for x in countArray]
    return profile


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    # get the first k-mers from each DNA to form an intitial array of best motifs
    for index in range(0, t):
        BestMotifs.append(Dna[index][0:k])
    # going trough all the possible k-mers in DNA[O] and forming motifs matrix from the rest dna strings
    # finding the best motifs by caluclating the best score (which is the lowest score computed for the motifs)
    n = len(Dna[0])
    for j in range(0, n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][j:j + k])
        for l in range(1, t):
            profile = ProfileWithPseudocounts(Motifs[0:l], l)
            Motifs.append(ProfileMostProbableKmer(Dna[l], k, profile))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def Motifs(Profile, Dna):
    k = len(Profile["A"])
    n = len(Dna)
    motifs = []
    for index in range(0, n):
        dnaString = Dna[index]
        motifs.append(ProfileMostProbableKmer(dnaString, k, Profile))
    return motifs


def RandomMotifs(Dna, k, t):
    motifs = []
    for index in range(0, t):
        dnaString = Dna[index]
        startMotifIndex = random.randint(1, k)
        endMotifIndex = startMotifIndex + k
        randomMotif = dnaString[startMotifIndex:endMotifIndex]
        motifs.append(randomMotif)
    return motifs


def RandomizedMotifSearch(Dna, k, t):
    motifs = RandomMotifs(Dna, k, t)
    BestMotifs = motifs
    count = 1
    while True:
        profile = ProfileWithPseudocounts(motifs, count)
        motifs = Motifs(profile, Dna)
        if Score(motifs) < Score(BestMotifs):
            BestMotifs = motifs
        else:
            return BestMotifs
        count = +1


def RepeatedRandomizedMotifSearch(Dna, k, t, N):
    BestMotifs = RandomMotifs(Dna, k, t)
    for index in range(0, N):
        motifs = RandomizedMotifSearch(Dna, k, t)
        if Score(motifs) < Score(BestMotifs):
            BestMotifs = motifs
    return BestMotifs


def Normalize(Probabilities):
    total_count = sum(Probabilities.values())
    for key in Probabilities:
        Probabilities[key] = Probabilities[key] / total_count
    return Probabilities


def WeightedDie(Probabilities):
    weight = random.uniform(0, 1)
    sorted_probabilities = dict(sorted(Probabilities.items(), key=lambda x: x[1]))
    for key in sorted_probabilities:
        weight = weight - sorted_probabilities[key]
        if weight < 0:
            return key


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


def GibbsSampler(Dna, k, t, N):
    BestMotifs = []
    random_motifs = RandomMotifs(Dna, k, t)
    BestMotifs = random_motifs
    for j in range(N):
        index = random.randint(0, t - 1)
        deleted_motif = random_motifs.pop(index)
        profile_matrix = ProfileWithPseudocounts(random_motifs, 4)
        kmer = ProfileGeneratedString(deleted_motif, profile_matrix, k)
        random_motifs.insert(index, kmer)
        if Score(random_motifs) < Score(BestMotifs):
            BestMotifs = random_motifs
    return BestMotifs
