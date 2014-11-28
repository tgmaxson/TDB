from __future__ import division
from fractions import gcd
import re
import os

def DRange(start, stop, step):
    '''
    Helper generator that can handle float ranges
    '''
    r = start
    while r < stop:
        yield r
        r += step
        
def CountList(lst, a):
    '''
    Helper function that can count the occurrences of an object in a list
    '''
    x = 0
    for i in lst:
        if i == a:
            x+=1
    return x

def BreakdownFormula(x, red=True):
    '''
    Used by FormatForm to break down and reduce the formula
    '''
    
    form = {}
    for i in re.findall('([A-Z][a-z]{0,2})(\d*)', x):
        if len(i)>0:
            if not i[0] in form:
                form[i[0]]=0
            if not i[1]=="":
                form[i[0]]+=int(i[1])
            else:
                form[i[0]]+=1
    if not red:
        return form
    if len(form)==1:
        for i in form:
            form[i]=1
            return form
    n = []
    for x in form:
        n.append(form[x])
    GCD = gcd(n[0],n[1])
    if len(form) > 2:
        for i in range(0,len(n)-1):
            GCD = gcd(GCD,n[i])
    for i in form:
        form[i]/=GCD    
    return form

def FormatFormula(x, formula_prefix, formula_suffix, red=True):
    '''
    Formats a formula to be ordered properly and reduced by default
    '''
    
    if x == "":
        return ""
    form = BreakdownFormula(x, red)
    ret = ""
    for prefix in formula_prefix:
        try:
            if form[prefix]==1:
                ret+=prefix
                form[prefix]=0
            elif form[prefix]>1:
                ret+=prefix+str(int(form[prefix]))
                form[prefix]=0
        except:
            pass
    for current in form:
        if current in formula_suffix:
            continue
        if form[current]==1:
            ret+=current
        elif form[current]>0:
            ret+=current+str(int(form[current]))
    for suffix in formula_suffix:
        try:
            if form[suffix]==1:
                ret+=suffix
                form[suffix]=0
            elif form[suffix]>1:
                ret+=suffix+str(int(form[suffix]))
                form[suffix]=0
        except:
            pass
    return ret
    
def ReportMessage(mess, ex=True):
    '''
    Used to exit the script in case of an error of some sort.  
    Prints a message,  exits and also gives the directory being worked on.
    '''
    
    print mess
    os.system("echo '"+mess+"' >> "+os.path.expanduser("~/db.log"))
    if ex:
        exit()