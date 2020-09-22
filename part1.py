import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def H(x):
    if x <= 0:
        return 0
    else:
        return 1


def coef4(a, x=0):
    return H(x-a)*((x-a)**4)/(24*EI)


def coef3(a, x=0):
    return H(x-a)*((x-a)**3)/(6*EI)


def coef2(a, x=0):
    return H(x-a)*((x-a)**2)/(2*EI)


def coef1(a, x=0):
    return H(x-a)*(x-a)/EI


def coef(a, x=0):
    return H(x-a)/EI


class Pivot:
    def __init__(self, crdnt=0):
        self.crdnt = crdnt


class Torque:
    def __init__(self, crdnt=0, value=0):
        self.crdnt = crdnt
        self.value = value


class Force:
    def __init__(self, crdnt=0, value=0):
        self.crdnt = crdnt
        self.value = value


class DistributedPower:
    def __init__(self, lborder=0, rborder=0, value=0):
        self.lborder = lborder
        self.rborder = rborder
        self.value = value


class ElasticForce:
    def __init__(self, crdnt=0, k_coef=0):
        self.crdnt = crdnt
        self.k_coef = k_coef


#def BorderDeterminator(type=0):
#    if type == 1:
#        W0 = 0
#        Theta0 = 1
#        M0 = 0
#        R0 = 1
#    elif type == 2:



E = 2*10**11
I = 0.2 * (0.1 ** 3) / 12
EI = E * I

inFile = open('input.txt', 'r', encoding='utf8')
inFileLines = inFile.readlines()
lbordertype, rbordertype = map(int, inFileLines[0].split())
#nmbrOfPivots = int(inFileLines[0])
# R = [Pivot()]*nmbrOfPivots
R = []
pivotsCrdnt = list(map(float, inFileLines[1].split()))
i = 0
#if lbordertype == 3:
#    R.append(Pivot(0))
while i < len(pivotsCrdnt):
    R.append(Pivot(pivotsCrdnt[i]))
    i += 1
M = list()
i = 0
torqueCrdnt = list(map(float, inFileLines[2].split()))
while i < len(torqueCrdnt):
    M.append(Torque(torqueCrdnt[i], torqueCrdnt[i+1]))
    i += 2
P = list()
i = 0
forceCrdnt = list(map(float, inFileLines[3].split()))
while i < len(forceCrdnt):
    P.append(Force(forceCrdnt[i], forceCrdnt[i+1]))
    i += 2
q = list()
i = 0
dispowerCrdnt = list(map(float, inFileLines[4].split()))
while i < len(dispowerCrdnt):
    q.append(DistributedPower(dispowerCrdnt[i], dispowerCrdnt[i+1],dispowerCrdnt[i+2]))
    i += 3
Ru = list()
i = 0
elasticforceCrdnt = list(map(float, inFileLines[5].split()))
while i < len(elasticforceCrdnt):
    Ru.append(ElasticForce(elasticforceCrdnt[i], elasticforceCrdnt[i+1]))
    i += 2
if rbordertype == 4:
    end4 = int(inFileLines[6])
# можно поставить счетчик недостающих уравнений и вычитать при добавлении строк в А
A = []
C = []
for i in range(len(R)):  # W = 0 in pivots
    if R[i].crdnt != 0:
        Temp = []
        tempC = 0
        if (rbordertype == 1) or (rbordertype == 2):
            for j in range(len(R) - 1):
                Temp.append(coef3(R[j].crdnt, R[i].crdnt))
        else:
            for j in range(len(R)):
                Temp.append(coef3(R[j].crdnt, R[i].crdnt))
        for j in range(len(Ru)):
            Temp.append(coef3(Ru[j].crdnt, R[i].crdnt))
        Temp.append(1)
        Temp.append(coef1(0, R[i].crdnt))
        Temp.append(coef2(0, R[i].crdnt))
        A.append(Temp)
        for j in range(len(M)):
            tempC -= M[j].value * coef2(M[j].crdnt, R[i].crdnt)
        for j in range(len(P)):
            tempC -= P[j].value * coef3(P[j].crdnt, R[i].crdnt)
        for j in range(len(q)):
            tempC += q[j].value * (coef4(q[j].rborder, R[i].crdnt) - coef4(q[j].lborder, R[i].crdnt))
        C.append(tempC)

for i in range(len(Ru)):  # W = 0 in pivots
    Temp = []
    tempC = 0
    if (rbordertype == 1) or (rbordertype == 2):
        for j in range(len(R) - 1):
            Temp.append(Ru[i].k_coef * coef3(R[j].crdnt, Ru[i].crdnt))
    else:
        for j in range(len(R)):
            Temp.append(Ru[i].k_coef * coef3(R[j].crdnt, Ru[i].crdnt))
    for j in range(len(Ru)):
        if j == i:
            Temp.append(Ru[i].k_coef * coef3(Ru[j].crdnt, Ru[i].crdnt) - 1)  # =0
        else:
            Temp.append(Ru[i].k_coef * coef3(Ru[j].crdnt, Ru[i].crdnt))
    Temp.append(Ru[i].k_coef)
    Temp.append(Ru[i].k_coef * coef1(0, Ru[i].crdnt))
    Temp.append(Ru[i].k_coef * coef2(0, Ru[i].crdnt))
    A.append(Temp)
    for j in range(len(M)):
        tempC -= Ru[i].k_coef * M[j].value * coef2(M[j].crdnt, R[i].crdnt)
    for j in range(len(P)):
        tempC -= Ru[i].k_coef * P[j].value * coef3(P[j].crdnt, R[i].crdnt)
    for j in range(len(q)):
        tempC += Ru[i].k_coef * q[j].value * (coef4(q[j].rborder, R[i].crdnt) - coef4(q[j].lborder, R[i].crdnt))
    C.append(tempC)

#  остановился на добавлении столбцов для тета0 омега0 и м0 в предыдущие два блока формирования матрицы А

#Temp = [] # добавление уравнение R=kW
#tempC = 0
#for j in range(len(R)-1):
#    Temp.append(Ru[0].k_coef * coef3(R[j].crdnt, Ru[0].crdnt))
##for j in range(len(Ru)):
##    Temp.append(coef3(Ru[j].crdnt, R[i].crdnt))
#Temp.append(Ru[0].k_coef * coef3(Ru[0].crdnt, Ru[0].crdnt) - 1)  # = 0
#Temp.append(Ru[0].k_coef * coef1(0, Ru[0].crdnt))
#A.append(Temp)
#tempC = Ru[0].k_coef * (- M[0].value * coef2(M[0].crdnt, Ru[0].crdnt) -
#                        P[0].value * coef3(P[0].crdnt, Ru[0].crdnt) +
#                        q[0].value * (coef4(q[0].rborder, Ru[0].crdnt) - coef4(q[0].lborder, Ru[0].crdnt)))
#C.append(tempC)

# добавление строк начальных условий

if lbordertype == 1:  # W0=0; M0=0
    Temp = []
    tempC = 0
    for i in range(len(R) + len(Ru) - 1 + H(rbordertype - 2)):
        Temp.append(0)
    Temp.append(1)
    Temp.append(0)
    Temp.append(0)
    A.append(Temp)
    C.append(0)
    Temp = []
    tempC = 0
    for i in range(len(R) + len(Ru) - 1 + H(rbordertype - 2)):
        Temp.append(0)
    Temp.append(0)
    Temp.append(0)
    Temp.append(1)
    A.append(Temp)
    C.append(0)
elif lbordertype == 2:  # M0=0
    Temp = []
    tempC = 0
    for i in range(len(R) + len(Ru) - 1 + H(rbordertype - 2)):
        Temp.append(0)
    Temp.append(1)
    Temp.append(0)
    Temp.append(0)
    A.append(Temp)
    C.append(0)
    Temp = []
    tempC = 0
    for i in range(len(R) + len(Ru) - 1 + H(rbordertype - 2)):
        Temp.append(0)
    Temp.append(0)
    Temp.append(1)
    Temp.append(0)
    A.append(Temp)
    C.append(0)
elif lbordertype == 3:  # W0=0; Theta0=0
    Temp = []
    tempC = 0
    for i in range(len(R) + len(Ru) - 1 + H(rbordertype - 2)):
        Temp.append(0)
    Temp.append(0)
    Temp.append(0)
    Temp.append(1)
    A.append(Temp)
    C.append(0)
elif lbordertype == 4:  # M0=0
    Temp = []
    tempC = 0
    for i in range(len(R) + len(Ru) - 1):
        Temp.append(0)
    Temp.append(0)
    Temp.append(0)
    Temp.append(1)
    A.append(Temp)
    C.append(0)


if rbordertype == 1:  # M(20)=0; 20 - sharnir
    Temp = []
    tempC = 0
    if rbordertype == 1:
        for j in range(len(R) - 1):
            Temp.append(coef1(R[j].crdnt, R[-1].crdnt))
    else:
        for j in range(len(R)):
            Temp.append(coef1(R[j].crdnt, R[-1].crdnt))
    for j in range(len(Ru)):
        Temp.append(coef1(Ru[j].crdnt, R[-1].crdnt))
    Temp.append(0)
    Temp.append(0)
    Temp.append(coef(0, R[-1].crdnt))
    A.append(Temp)
    for j in range(len(M)):
        tempC -= M[j].value * coef(M[j].crdnt, R[-1].crdnt)
    for j in range(len(P)):
        tempC -= P[j].value * coef1(P[j].crdnt, R[-1].crdnt)
    for j in range(len(q)):
        tempC += q[j].value * (coef2(q[j].rborder, R[-1].crdnt) - coef2(q[j].lborder, R[-1].crdnt))
    C.append(tempC)
elif rbordertype == 2:  # Theta(20)=0; 20- консоль
    Temp = []
    tempC = 0
    if (rbordertype == 1) or (rbordertype == 2):
        for j in range(len(R) - 1):
            Temp.append(coef2(R[j].crdnt, R[-1].crdnt))
    else:
        for j in range(len(R)):
            Temp.append(coef2(R[j].crdnt, R[-1].crdnt))
    for j in range(len(Ru)):
        Temp.append(coef2(Ru[j].crdnt, R[-1].crdnt))
    Temp.append(0)
    Temp.append(coef(0, R[-1].crdnt))
    Temp.append(coef1(0, R[-1].crdnt))
    A.append(Temp)
    for j in range(len(M)):
        tempC -= M[j].value * coef1(M[j].crdnt, R[-1].crdnt)
    for j in range(len(P)):
        tempC -= P[j].value * coef2(P[j].crdnt, R[-1].crdnt)
    for j in range(len(q)):
        tempC += q[j].value * (coef3(q[j].rborder, R[-1].crdnt) - coef3(q[j].lborder, R[-1].crdnt))
    C.append(tempC)
elif rbordertype == 3:  # M(20)=0, R=kW
    Temp = []
    tempC = 0
    for j in range(len(R)):
        Temp.append(coef1(R[j].crdnt, Ru[-1].crdnt))
    for j in range(len(Ru)):
        Temp.append(coef1(Ru[j].crdnt, Ru[-1].crdnt))
    Temp.append(0)
    Temp.append(0)
    Temp.append(coef(0, Ru[-1].crdnt))
    A.append(Temp)
    for j in range(len(M)):
        tempC -= M[j].value * coef(M[j].crdnt, Ru[-1].crdnt)
    for j in range(len(P)):
        tempC -= P[j].value * coef1(P[j].crdnt, Ru[-1].crdnt)
    for j in range(len(q)):
        tempC += q[j].value * (coef2(q[j].rborder, Ru[-1].crdnt) - coef2(q[j].lborder, Ru[-1].crdnt))
    C.append(tempC)
    Temp = []
    tempC = 0
    for j in range(len(R)):
        Temp.append(coef(R[j].crdnt, Ru[-1].crdnt) - Ru[-1].k_coef * coef2(R[j].crdnt, Ru[-1].crdnt))
    for j in range(len(Ru)):
        Temp.append(coef(Ru[j].crdnt, Ru[-1].crdnt) - Ru[-1].k_coef * coef2(Ru[j].crdnt, Ru[-1].crdnt))
    Temp.append(0)
    Temp.append(0)
    Temp.append(0)
    A.append(Temp)
    for j in range(len(M)):
        tempC -= M[j].value * (- Ru[-1].k_coef * coef1(M[j].crdnt, Ru[-1].crdnt))
    for j in range(len(P)):
        tempC -= P[j].value * (coef(P[j].crdnt, Ru[-1].crdnt) - Ru[-1].k_coef * coef2(P[j].crdnt, Ru[-1].crdnt))
    for j in range(len(q)):
        tempC += q[j].value * (coef1(q[j].rborder, Ru[-1].crdnt) - coef1(q[j].lborder, Ru[-1].crdnt) -\
                               Ru[-1].k_coef * coef3(q[j].rborder, Ru[-1].crdnt) +\
                               Ru[-1].k_coef * coef3(q[j].lborder, Ru[-1].crdnt))
    C.append(tempC)
elif rbordertype == 4:  # M()=0, Q()=0
    Temp = []
    tempC = 0
    for j in range(len(R)):
        Temp.append(coef1(R[j].crdnt, end4))
    for j in range(len(Ru)):
        Temp.append(coef1(Ru[j].crdnt, end4))
    Temp.append(0)
    Temp.append(0)
    Temp.append(coef(0, end4))
    A.append(Temp)
    for j in range(len(M)):
        tempC -= M[j].value * coef(M[j].crdnt, end4)
    for j in range(len(P)):
        tempC -= P[j].value * coef1(P[j].crdnt, end4)
    for j in range(len(q)):
        tempC += q[j].value * (coef2(q[j].rborder, end4) - coef2(q[j].lborder, end4))
    C.append(tempC)
    Temp = []
    tempC = 0
    for j in range(len(R)):
        Temp.append(coef(R[j].crdnt, end4))
    for j in range(len(Ru)):
        Temp.append(coef(Ru[j].crdnt, end4))
    Temp.append(0)
    Temp.append(0)
    Temp.append(0)
    A.append(Temp)
    for j in range(len(M)):
        tempC -= M[j].value * 0  #coef(M[j].crdnt, end4)
    for j in range(len(P)):
        tempC -= P[j].value * coef(P[j].crdnt, end4)
    for j in range(len(q)):
        tempC += q[j].value * (coef1(q[j].rborder, end4) - coef1(q[j].lborder, end4))
    C.append(tempC)


#Temp = []
#tempC = 0
#for j in range(len(R) - 1):
#    Temp.append(coef2(R[j].crdnt, R[3].crdnt))
#for j in range(len(Ru)):
#    Temp.append(coef3(Ru[j].crdnt, R[i].crdnt))
#Temp.append(coef(0, R[3].crdnt))
#A.append(Temp)
#tempC = - M[0].value * coef1(M[0].crdnt, R[3].crdnt) - \
#        P[0].value * coef2(P[0].crdnt, R[3].crdnt) + \
#        q[0].value * (coef3(q[0].rborder, R[3].crdnt) - coef3(q[0].lborder, R[3].crdnt))
#C.append(tempC)

A = np.array(A)
C = np.array(C)

S = list(range(len(R) + len(Ru) + 2))


S = np.linalg.solve(A, C)


def W1(x):
    temp = 0
    j = 0
    if rbordertype == 1:
        for i in range(len(R) - 1):
            temp += S[j] * coef3(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(R)):
            temp += S[j] * coef3(R[i].crdnt, x)
            j += 1
    if rbordertype == 2:
        for i in range(len(Ru) - 1):
            temp += S[j] * coef3(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(Ru)):
            temp += S[j] * coef3(Ru[i].crdnt, x)
            j += 1
    temp += S[j]
    temp += S[j+1] * coef1(0, x)
    temp += S[j+2] * coef2(0, x)
    for i in range(len(q)):
        temp += q[i].value * (coef4(q[i].lborder, x) - coef4(q[i].rborder, x))
    for i in range(len(P)):
        temp += P[i].value * coef3(P[i].crdnt, x)
    for i in range(len(M)):
        temp += M[i].value * coef2(M[i].crdnt, x)
    return temp


def dW1(x):
    temp = 0
    j = 0
    if rbordertype == 1:
        for i in range(len(R) - 1):
            temp += S[j] * coef2(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(R)):
            temp += S[j] * coef2(R[i].crdnt, x)
            j += 1
    if rbordertype == 2:
        for i in range(len(Ru) - 1):
            temp += S[j] * coef2(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(Ru)):
            temp += S[j] * coef2(Ru[i].crdnt, x)
            j += 1
    temp += 0
    temp += S[j+1] * coef(0, x)
    temp += S[j+2] * coef1(0, x)
    for i in range(len(q)):
        temp += q[i].value * (coef3(q[i].lborder, x) - coef3(q[i].rborder, x))
    for i in range(len(P)):
        temp += P[i].value * coef2(P[i].crdnt, x)
    for i in range(len(M)):
        temp += M[i].value * coef1(M[i].crdnt, x)
    return temp


def d2W1(x):
    temp = 0
    j = 0
    if rbordertype == 1:
        for i in range(len(R) - 1):
            temp += S[j] * coef1(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(R)):
            temp += S[j] * coef1(R[i].crdnt, x)
            j += 1
    if rbordertype == 2:
        for i in range(len(Ru) - 1):
            temp += S[j] * coef1(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(Ru)):
            temp += S[j] * coef1(Ru[i].crdnt, x)
            j += 1
    temp += 0
    temp += 0
    temp += S[j+2] * coef(0, x)
    for i in range(len(q)):
        temp += q[i].value * (coef2(q[i].lborder, x) - coef2(q[i].rborder, x))
    for i in range(len(P)):
        temp += P[i].value * coef1(P[i].crdnt, x)
    for i in range(len(M)):
        temp += M[i].value * coef(M[i].crdnt, x)
    return temp


def d3W1(x):
    temp = 0
    j = 0
    if rbordertype == 1:
        for i in range(len(R) - 1):
            temp += S[j] * coef(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(R)):
            temp += S[j] * coef(R[i].crdnt, x)
            j += 1
    if rbordertype == 2:
        for i in range(len(Ru) - 1):
            temp += S[j] * coef(R[i].crdnt, x)
            j += 1
    else:
        for i in range(len(Ru)):
            temp += S[j] * coef(Ru[i].crdnt, x)
            j += 1
    temp += 0
    temp += 0
    temp += 0
    for i in range(len(q)):
        temp += q[i].value * (coef1(q[i].lborder, x) - coef1(q[i].rborder, x))
    for i in range(len(P)):
        temp += P[i].value * coef(P[i].crdnt, x)
#    for i in range(len(M)):
#        temp += 0
    return temp

#def dW1(x):
#    return S[0] * coef2(R[0].crdnt, x) + S[1] * coef2(R[1].crdnt, x) + S[2] * coef2(R[2].crdnt, x) +\
#           S[3] * coef2(R[3].crdnt, x) + S[4] * coef(0, x) +\
#           q[0].value * (coef3(q[0].lborder, x) - coef3(q[0].rborder, x)) +\
#           P[0].value * coef2(P[0].crdnt, x) + M[0].value * coef1(M[0].crdnt, x)


#def d2W1(x):
#    return S[0] * coef1(R[0].crdnt, x) + S[1] * coef1(R[1].crdnt, x) + S[2] * coef1(R[2].crdnt, x) +\
#            S[3] * coef1(R[3].crdnt, x) +\
#           q[0].value * (coef2(q[0].lborder, x) - coef2(q[0].rborder, x)) +\
#           P[0].value * coef1(P[0].crdnt, x) + M[0].value * coef(M[0].crdnt, x)


#def d3W1(x):
#    return S[0] * coef(R[0].crdnt, x) + S[1] * coef(R[1].crdnt, x) + S[2] * coef(R[2].crdnt, x) +\
#           S[3] * coef(R[3].crdnt, x) +\
#           q[0].value * (coef1(q[0].lborder, x) - coef1(q[0].rborder, x)) +\
#           P[0].value * coef(P[0].crdnt, x)


fig1, fig1_axes = plt.subplots(ncols=1, nrows=4)
GRID_STEP = 0.1


def create_graph1(func1, N, name):
    y = []
    GRID_STEP2 = GRID_STEP/10
    x = np.arange(0, 20, GRID_STEP2)
    for i in x:
        y.append(func1(i))
    fig1_axes[N].set_title(name)
    fig1_axes[N].plot(x, y, c='green')
    fig1_axes[N].xaxis.set_major_locator(MultipleLocator(base=1))
    fig1_axes[N].grid()


create_graph1(W1, 0, 'W(x) - прогиб')
create_graph1(dW1, 1, 'Theta(x) - угол поворота')
create_graph1(d2W1, 2, 'M(x) - изгибающий момент')
create_graph1(d3W1, 3, 'Q(x) - перерезывающая сила')


plt.savefig('figure1.png')
pass
