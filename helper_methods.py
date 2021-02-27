
# Counts how many time Pattern is contained in Text
def patternCount(Text, Pattern):
    count = 0
    for index in range(0, len(Text) - len(Pattern) + 1):
        if Text[index:index + len(Pattern)] == Pattern:
            count = count + 1
    return count


# Most frequent k-mer in Text
def frequentWord(k, Text):
    dicionary_of_frequent_words = dict()
    for index in range(0, len(Text) - k):
        pattern = Text[index:index + k]
        dicionary_of_frequent_words[pattern] = patternCount(Text, pattern)
    max_frequency = max(dicionary_of_frequent_words.values())
    most_frequent_words = []
    for pattern, count in dicionary_of_frequent_words.items():
        if count == max_frequency:
            most_frequent_words.append(pattern)
    return most_frequent_words


# Frequency map of k-mers in Text
def frequencyMap(Text, k):
    dicionary_of_frequent_words = dict()
    for index in range(0, len(Text) - k + 1):
        pattern = Text[index:index + k]
        dicionary_of_frequent_words[pattern] = patternCount(Text, pattern)
    return dicionary_of_frequent_words


# Reverse DNA strand
def Reverse(Pattern):
    str = ""
    for index in range(len(Pattern), 0, -1):
        str += Pattern[index - 1]
    return str


# Complement DNA strand
def Complement(Pattern):
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
