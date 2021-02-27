# Reverse DNA strand
def reverse(Pattern):
    str = ""
    for index in range(len(Pattern), 0, -1):
        str += Pattern[index - 1]
    return str


# Complement DNA strand
def complement(Pattern):
    str = ""
    for char in Pattern:
        if char == "A":
            str += "T"
        elif char == "T":
            str += "A"
        elif char == "C":
            str += "G"
        else:
            str += "C"
    return str


def reverseComplement(Pattern):
    complementarty_strand = reverse(Pattern)
    return complement(complementarty_strand)


# Returns the location of Pattern in Genome
def patternMatching(Pattern, Genome):
    location_of_pattern = []
    for index in range(0, len(Genome) - len(Pattern) + 1):
        if Genome[index:index + len(Pattern)] == Pattern:
            location_of_pattern.append(index)
    return location_of_pattern
