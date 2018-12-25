
# coding: utf-8

# Find Rice Gene Name V2

# namecheck : Check your input is MSU-ID or RAP-DB number  
def namecheck(name):
      if 'LOC_Os' in name :
            return 1
      elif 'g' in name :
            return 0
      else :
            return 2

# Copy to the chipBoard
import os
def addToClipBoard(text):
    command = 'echo ' + text.strip() + '| clip'
    os.system(command)

# Open file and return a list
def OpenFile(filename):
    with open(filename,'r') as file:
        Lst = file.readlines()
    for i in range(len(Lst)):
        Lst[i] = Lst[i].strip()
        Lst[i] = Lst[i].split('\t')

    return Lst

# Open data and store as list 
# RAP_MSU list
RM_lst = OpenFile("RAP_MSU.txt")
# RAP_Name list
RN_lst = OpenFile("RAP_name.txt")
# RAP_Transcrip evidence 
RT_lst = OpenFile("RAP_Transcript.txt")


# Search MSU, Transcript evidence and Gene name


def Findinlst(RAPID,Lst):
    Lstnew =[]
    for i in range(len(Lst)):
        if Lst[i][0] == RAPID:
            Lstnew.append(Lst[i][1])
    return ",".join(Lstnew)
def quickFindinlst(RAPID):
    Lst = ['RAP','MSU','Transcript','Name']
    Lst[0] = RAPID    
    Lst[1] = Findinlst(RAPID,RM_lst)
    Lst[2] = Findinlst(RAPID,RT_lst)
    Lst[3] = Findinlst(RAPID,RN_lst)
    return Lst


# You can unput MSU-ID、RAPDB number、transcript evidence


def Result(ID):
    ResultLst = ['RAP','MSU','Transcript','Name']
    check = namecheck(ID)
    if check == 0:
        # Your input is RAPDB-number
         ResultLst = quickFindinlst(ID)
         #addToClipBoard(ResultLst[1])
    elif check == 1:
        # Your input is MSU-ID, Exchange to RAP number
        # Exchange it to RAP-DB number
        for i in range(len(RM_lst)):
            if RM_lst[i][1] == ID:
                ID = RM_lst[i][0]
        ResultLst = quickFindinlst(ID)
        # addToClipBoard(ResultLst[0])
    else:
        # Your input is Transctipt Evidence
        # Exchange it to RAP-DB number
        for i in range(len(RT_lst)):
            if RT_lst[i][1] == ID:
                ID = RT_lst[i][0]
        ResultLst = quickFindinlst(ID)
        # addToClipBoard(ResultLst[1])       
    return ResultLst






