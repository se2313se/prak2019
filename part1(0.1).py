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


E = 2*10**11
I = 0.2 * (0.1 ** 3) / 12
EI = E * I

inFile = open('input.txt', 'r', encoding='utf8')
inFileLines = inFile.readlines()
nmbrOfPivots = int(inFileLines[0])
# R = [Pivot()]*nmbrOfPivots
R = []
pivotsCrdnt = list(map(float, inFileLines[1].split()))
for i in range(nmbrOfPivots):
    R.append(Pivot(pivotsCrdnt[i]))
M = list()
M.append(Torque(list(map(float, inFileLines[2].split()))[0], list(map(float, inFileLines[2].split()))[1]))
P = list()
P.append(Force(list(map(float, inFileLines[3].split()))[0], list(map(float, inFileLines[3].split()))[1]))
q = list()
q.append(DistributedPower(list(map(float, inFileLines[4].split()))[0],
                          list(map(float, inFileLines[4].split()))[1],
                          list(map(float, inFileLines[4].split()))[2]))
Ru = list()
Ru.append(ElasticForce(list(map(float, inFileLines[5].split()))[0], list(map(float, inFileLines[5].split()))[1]))
# можно поставить счетчик недостающих уравнений и вычитать при добавлении строк в А
A = []
C = []
for i in range(1, len(R)):  # W = 0 in pivots
    Temp = []
    tempC = 0
    for j in range(len(R) - 1):
        Temp.append(coef3(R[j].crdnt, R[i].crdnt))
    Temp.append(coef3(Ru[0].crdnt, R[i].crdnt))
    Temp.append(coef1(0, R[i].crdnt))
    A.append(Temp)
    tempC = - M[0].value * coef2(M[0].crdnt, R[i].crdnt) - \
            P[0].value * coef3(P[0].crdnt, R[i].crdnt) + \
            q[0].value * (coef4(q[0].rborder, R[i].crdnt) - coef4(q[0].lborder, R[i].crdnt))
    C.append(tempC)

Temp = []
tempC = 0
for j in range(len(R) - 1):
    Temp.append(coef2(R[j].crdnt, R[3].crdnt))
Temp.append(coef2(Ru[0].crdnt, R[3].crdnt))
Temp.append(coef(0, R[3].crdnt))
A.append(Temp)
tempC = - M[0].value * coef1(M[0].crdnt, R[3].crdnt) - \
        P[0].value * coef2(P[0].crdnt, R[3].crdnt) + \
        q[0].value * (coef3(q[0].rborder, R[3].crdnt) - coef3(q[0].lborder, R[3].crdnt))
C.append(tempC)

Temp = []
tempC = 0
for j in range(len(R)-1):
    Temp.append(Ru[0].k_coef * coef3(R[j].crdnt, Ru[0].crdnt))
Temp.append(Ru[0].k_coef * coef3(Ru[0].crdnt, Ru[0].crdnt) - 1)  # = 0
Temp.append(Ru[0].k_coef * coef1(0, Ru[0].crdnt))
A.append(Temp)
tempC = Ru[0].k_coef * (- M[0].value * coef2(M[0].crdnt, Ru[0].crdnt) -
                        P[0].value * coef3(P[0].crdnt, Ru[0].crdnt) +
                        q[0].value * (coef4(q[0].rborder, Ru[0].crdnt) - coef4(q[0].lborder, Ru[0].crdnt)))
C.append(tempC)

A = np.array(A)
C = np.array(C)

S = list(range(5))


S = np.linalg.solve(A, C)


def W1(x):
    return S[0] * coef3(R[0].crdnt, x) + S[1] * coef3(R[1].crdnt, x) + S[2] * coef3(R[2].crdnt, x) +\
           S[3] * coef3(R[3].crdnt, x) + S[4] * coef1(0, x) +\
           q[0].value * (coef4(q[0].lborder, x) - coef4(q[0].rborder, x)) +\
           P[0].value * coef3(P[0].crdnt, x) + M[0].value * coef2(M[0].crdnt, x)


def dW1(x):
    return S[0] * coef2(R[0].crdnt, x) + S[1] * coef2(R[1].crdnt, x) + S[2] * coef2(R[2].crdnt, x) +\
           S[3] * coef2(R[3].crdnt, x) + S[4] * coef(0, x) +\
           q[0].value * (coef3(q[0].lborder, x) - coef3(q[0].rborder, x)) +\
           P[0].value * coef2(P[0].crdnt, x) + M[0].value * coef1(M[0].crdnt, x)


def d2W1(x):
    return S[0] * coef1(R[0].crdnt, x) + S[1] * coef1(R[1].crdnt, x) + S[2] * coef1(R[2].crdnt, x) +\
           S[3] * coef1(R[3].crdnt, x) +\
           q[0].value * (coef2(q[0].lborder, x) - coef2(q[0].rborder, x)) +\
           P[0].value * coef1(P[0].crdnt, x) + M[0].value * coef(M[0].crdnt, x)


def d3W1(x):
    return S[0] * coef(R[0].crdnt, x) + S[1] * coef(R[1].crdnt, x) + S[2] * coef(R[2].crdnt, x) +\
           S[3] * coef(R[3].crdnt, x) +\
           q[0].value * (coef1(q[0].lborder, x) - coef1(q[0].rborder, x)) +\
           P[0].value * coef(P[0].crdnt, x)


fig1, fig1_axes = plt.subplots(ncols=1, nrows=4)
GRID_STEP = 0.1


def create_graph(func1, func2, N, name):
    y = []
    y1 = []
    GRID_STEP2 = GRID_STEP/10
    x = np.arange(0, 20, GRID_STEP2)
    for i in x:
        y.append(func1(i))
        y1.append(func2(i))
    fig1_axes[N].set_title(name)
    fig1_axes[N].plot(x, y1, c='green')
    fig1_axes[N].plot(x, y, c='purple')
    fig1_axes[N].xaxis.set_major_locator(MultipleLocator(base=1))
    fig1_axes[N].grid()


create_graph(W1, W1, 0, 'W(x) - прогиб')
create_graph(dW1, dW1, 1, 'Theta(x) - угол поворота')
create_graph(d2W1, d2W1, 2, 'M(x) - изгибающий момент')
create_graph(d3W1, d3W1, 3, 'Q(x) - перерезывающая сила')


plt.savefig('figure1.png')
pass
