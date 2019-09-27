#!/usr/bin/env python3
import re,sys

#wrap python regex for easy use just passing in argv PATTERN, STRING(also splitted by shell, re joined with spaces)
if len(sys.argv)<2:
    print("usage PATTERN, STRING")

#join bash tokenization, newline that has separated tokens now are converted inspaces
def regex(pattern,string):
    return re.findall(pattern,string)
def joinRegex(pattern,stringTokenized):
    #join list of string (e.g. tokinzation of bash of redirected input" and apply pattern regex
    stringJoined=" ".join(stringTokenized)
    return regex(pattern,stringJoined)

if __name__=="__main__":
    PATTERN=sys.argv[1]
    STRING_TOKENIZED=sys.argv[2:]
    STRING=" ".join(sys.argv[2:]) 
    res=joinRegex(PATTERN,STRING_TOKENIZED)
    if len(res) ==0:
        print(" searching regex ",PATTERN)
        print(" inside string of len ",len(STRING))
        print("no result founded")
        exit()
    
    for r in res:
        print(r)
