# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 16:17:35 2017

@author: HumayunMunazza
"""
import numpy as np

#Finite Differences

def pointGen(n):
    if n > 0:
        temp2 = "+ " + str(int(n))
    else:
        temp2 = "- " + str(np.abs(int(n)))        
    
    temp1 = "f(x "
    temp3 = "h)"
    return str(temp1 + temp2 + temp3)

def taylorSeries(maxorder = 3):
    lhs = "f(x)"
    equals = " = "
    rhs = ""
    for i in range(0, maxorder + 1):
        term = " + a" + str(i)
        termPart2 = "(x - a)^" + str(i)
        rhs = rhs + term + termPart2
    final = lhs + equals + rhs
    return final

def taylorSeriesExpanded(maxorder = 3):
    lhs = pointGen(maxorder)
    order = list(lhs)
    if order[4] == "-":
        coef = order[4] + order[6]
    else:
        coef = order[6]    
    equals = " = "
    rhs = ""
    firstterm = "f(x)"
    for i in range(1, np.abs(maxorder) + 1):
        term = " + {f" + "'" * i
        termPart2 = "}(" + str(coef) + "h)^" + str(i)
        termPart3 = "[1/" + str(i) + "!]"
        rhs = rhs + term + termPart3 + termPart2
    final = lhs + equals + firstterm + rhs
    return final

# This expands the point of interest. For example 'x + 2h' or 'x + 1h'. pointNum
# controls the coefficient of 'h' and accuracy expands the taylor series to the 
# desired residual term. So if pointNum is 1 and accuracy is 2, then f(x + h) 
# should be expanded to the 2nd derivative. note that the residual term is the
# 3rd term.

def pointExpand(pointNum = 1, accuracy = 2):
    maxorder = accuracy
    term = np.zeros(maxorder)
    lhs = pointGen(pointNum)
    order = list(lhs)
    if order[4] == "-":
        coef = order[4] + order[6]
    else:
        coef = order[6]
    coef = float(coef)
    for i in range(0, maxorder):
        term[i] = pow(coef, i) * (1.0 / np.math.factorial(float(i)))
    return term

def finiteDiff(D = 1, StencilLen = 3, mode = "central"):
    if StencilLen % 2 == 0:
        StencilLen = StencilLen + 1
        
    order = StencilLen
    if D < StencilLen:
        pointNum = np.floor((StencilLen - 1)  / 2)
        CoefMat = np.zeros((StencilLen, StencilLen))
        pointNumInd = np.linspace(-1 * pointNum, pointNum, StencilLen, dtype = int)
        num = -1
        for i in pointNumInd:
            num = num + 1
            taylorExpand = pointExpand(i, StencilLen)
            np.copyto(CoefMat[num], taylorExpand)
            
        invCoefMat = np.linalg.inv(CoefMat.T)
        CoefVec = np.zeros(order)
        CoefVec[D] = 1
        result = np.dot(invCoefMat, CoefVec.T)
        return result
    else:
        StencilInc = D - StencilLen
        StencilLen = StencilLen + StencilInc
        Err = "Stencil Length must be greater than order of Derivative" 
        return Err + "use" + StencilLen

print finiteDiff(4, 5)



    













