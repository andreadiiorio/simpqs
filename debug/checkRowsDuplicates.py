#!/usr/bin/env python3

#check duplicates in rows using hashmap with K row
import sys
from pythonRegex import *
def check_duplicate_row(rows):
    #duplicates will be printed and number of occurrences returned
    dup=0
    rowDict=dict()
    for r in rows:
        if rowDict.get(r)!=None:
            print("DUPLICATION",r,file=sys.stderr)
            dup+=1
        #allert for null row-> inference row key kind by type of list passed element
        if type(r)==str:
            if "1" in r == False:
                print("NULL ROW",r)
        else:
            if "1" in r[1] == False:
                print("NULL ROW",r[1])
        rowDict[r]=True
    return dup

def tokenizeRows(rowsRaw,onlyRow):
    #build row tokens from rowsRaw list of lines as K [ROW 0/1]
    #if onlyRow is true a list of [ROW 0/1] will be returned
    #otherwise list of [(POL_VALUE,ROW 0/1)] returned
    rows=list()
    for rowRaw in rowsRaw:
        keyRaw=rowRaw.split(" ")
        #print("keyRaw",keyRaw)
        #print("rowRaw",rowRaw)
        if "" in keyRaw:
            keyRaw.remove("")
        if onlyRow:
            key="".join(keyRaw[1:])
        else:
            key=tuple(keyRaw)
        rows.append(key)
    return rows
REGEX_ROW_MATCH="> *([0-9]+ +[0-1]+) "
if __name__=="__main__":
    CHECK_DUPLICAT_ONLY_ON_ROW=False
    rowsRegexParsed=joinRegex(REGEX_ROW_MATCH,sys.argv[1:])
    rowsParsed=tokenizeRows(rowsRegexParsed,CHECK_DUPLICAT_ONLY_ON_ROW)
    #print(rowsParsed)
    duplicates=check_duplicate_row(rowsParsed)
    print(duplicates)

