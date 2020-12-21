# Nama : Theodorus Riyanto
# NIM : 12217017

import matplotlib.pyplot as plt
import csv, math, copy
import numpy as np
from scipy.linalg import lu_solve, lu_factor
from scipy.linalg import solve

data = []
pvtwater = []
densitystc = []
fluidrock = []
grid = []
pvtgas = []
reservoirrock = []
time_rate = []
well = []

def readdata():
    global wellint, Nw, wlx, wly, wlz
    with open("data.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            data.append(row)

    with open("pvtwater.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            pvtwater.append(row)

    with open("densitystc.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            densitystc.append(row)

    with open("fluidrock.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            fluidrock.append(row)

    with open("grid.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            grid.append(row)

    with open("pvtgas.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            pvtgas.append(row)

    with open("reservoirrock.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            reservoirrock.append(row)

    with open("Time_Rate.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            time_rate.append(row)

    with open("Well.csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            well.append(row)

    wellint = [int(i) for i in well[0]]
    Nw = wellint[0]
    wlx = [wellint[1]-1, wellint[4]-1]
    wly = [wellint[2]-1, wellint[5]-1]
    wlz = [wellint[3]-1, wellint[6]-1]

    return

def lookupoil(lookout,lookin,n):
    for i in range (len(data)):
        if (data[i][lookin-1]==n):
            return (data[i][lookout-1])
        elif (i == len(data)-1):
            if(data[0][lookin-1]<data[i][lookin-1] and n>data[i][lookin-1]):
                return ((n-data[i-1][lookin-1])/(data[i][lookin-1]-data[i-1][lookin-1])*(data[i][lookout-1]-data[i-1][lookout-1])+data[i-1][lookout-1])
            elif(data[0][lookin-1]>data[i][lookin-1] and n<data[i][lookin-1]):
                return ((n-data[i-1][lookin-1])/(data[i][lookin-1]-data[i-1][lookin-1])*(data[i][lookout-1]-data[i-1][lookout-1])+data[i-1][lookout-1])
            else:
                return ((n-data[0][lookin-1])/(data[1][lookin-1]-data[0][lookin-1])*(data[1][lookout-1]-data[0][lookout-1])+data[0][lookout-1])
        elif ((data[i][lookin-1]<n<data[i+1][lookin-1]) or (data[i+1][lookin-1]<n<data[i][lookin-1])):
            return ((n-data[i][lookin-1])/(data[i+1][lookin-1]-data[i][lookin-1])*(data[i+1][lookout-1]-data[i][lookout-1])+data[i][lookout-1])
    return()

def lookupgas(lookout,lookin,n):
    for i in range (len(pvtgas)):
        if (pvtgas[i][lookin-1]==n):
            return (pvtgas[i][lookout-1])
        elif (i == len(pvtgas)-1):
            if(pvtgas[0][lookin-1]<pvtgas[i][lookin-1] and n>pvtgas[i][lookin-1]):
                return ((n-pvtgas[i-1][lookin-1])/(pvtgas[i][lookin-1]-pvtgas[i-1][lookin-1])*(pvtgas[i][lookout-1]-pvtgas[i-1][lookout-1])+pvtgas[i-1][lookout-1])
            elif(pvtgas[0][lookin-1]>pvtgas[i][lookin-1] and n<pvtgas[i][lookin-1]):
                return ((n-pvtgas[i-1][lookin-1])/(pvtgas[i][lookin-1]-pvtgas[i-1][lookin-1])*(pvtgas[i][lookout-1]-pvtgas[i-1][lookout-1])+pvtgas[i-1][lookout-1])
            else:
                return ((n-pvtgas[0][lookin-1])/(pvtgas[1][lookin-1]-pvtgas[0][lookin-1])*(pvtgas[1][lookout-1]-pvtgas[0][lookout-1])+pvtgas[0][lookout-1])
        elif ((pvtgas[i][lookin-1]<n<pvtgas[i+1][lookin-1]) or (pvtgas[i+1][lookin-1]<n<pvtgas[i][lookin-1])):
            return ((n-pvtgas[i][lookin-1])/(pvtgas[i+1][lookin-1]-pvtgas[i][lookin-1])*(pvtgas[i+1][lookout-1]-pvtgas[i][lookout-1])+pvtgas[i][lookout-1])
    return()

def lookupfr(lookout,lookin,n):
    #Sw = 1
    #Krw = 2
    #Kro = 3
    #Pcow = 4
    for i in range(len(fluidrock)):
        if (fluidrock[i][lookin-1]==n):
            return (fluidrock[i][lookout-1])
        elif (i == len(fluidrock)-1):
            if(fluidrock[0][lookin-1]<fluidrock[i][lookin-1] and n>fluidrock[i][lookin-1]):
                return ((n-fluidrock[i-1][lookin-1])/(fluidrock[i][lookin-1]-fluidrock[i-1][lookin-1])*(fluidrock[i][lookout-1]-fluidrock[i-1][lookout-1])+fluidrock[i-1][lookout-1])
            elif(fluidrock[0][lookin-1]>fluidrock[i][lookin-1] and n<fluidrock[i][lookin-1]):
                return ((n-fluidrock[i-1][lookin-1])/(fluidrock[i][lookin-1]-fluidrock[i-1][lookin-1])*(fluidrock[i][lookout-1]-fluidrock[i-1][lookout-1])+fluidrock[i-1][lookout-1])
            else:
                return ((n-fluidrock[0][lookin-1])/(fluidrock[1][lookin-1]-fluidrock[0][lookin-1])*(fluidrock[1][lookout-1]-fluidrock[0][lookout-1])+fluidrock[0][lookout-1])
        elif ((fluidrock[i][lookin-1]<n<fluidrock[i+1][lookin-1]) or (fluidrock[i+1][lookin-1]<n<fluidrock[i][lookin-1])):
            return ((n-fluidrock[i][lookin-1])/(fluidrock[i+1][lookin-1]-fluidrock[i][lookin-1])*(fluidrock[i+1][lookout-1]-fluidrock[i][lookout-1])+fluidrock[i][lookout-1])
    return()

def fluidprop(lookout,p):
    #Bo
    if (lookout==1):
        return(lookupoil(3,2,p))
    #Visc,o
    elif (lookout==2):
        return(lookupoil(4,2,p))
    #Rs
    elif (lookout==3):
        return(lookupoil(1,2,p)/5.6146*1000)
    #Bg
    elif (lookout==4):
        return(lookupgas(2,1,p)*5.6146/1000)
    #Vg
    elif (lookout==5):
        return(lookupgas(3,1,p))
    #Bw
    elif (lookout==7):
        x=pvtwater[0][2]*(p-pvtwater[0][0])
        bw=pvtwater[0][1]/(1+x+x**2/2)
        return(bw)
    #Visc,w
    elif (lookout==8):
        y=-pvtwater[0][4]*(p-pvtwater[0][0])
        mw=pvtwater[0][3]/(1+y+y**2/2)
        return(mw)
    return()

def fluidprop2(lookout2,p):
    #Rho O
    if (lookout2==1):
        ro=(densitystc[0][0]+fluidprop(3,p)*densitystc[0][2])/fluidprop(1,p)
        return(ro/144)
    #Rho G
    elif (lookout2==2):
        rg=densitystc[0][2]/fluidprop(4,p)
        return(rg/144)
    #Rho W
    elif (lookout2==3):
        rw=densitystc[0][1]/fluidprop(7,p)
        return(rw/144)
    #Phi
    elif (lookout2==4):
        por=reservoirrock[0][3]*math.exp(reservoirrock[0][4]*(p-reservoirrock[0][2]))
        return(por)
    return()

def deriv(lookout,p):
    h=0.001
    y1=fluidprop(lookout,p)
    y2=fluidprop(lookout,p-h)
    dydp=(y1-y2)/h
    return(dydp)

def deriv2(lookout2,p):
    h=0.001
    y1=fluidprop2(lookout2,p)
    y2=fluidprop2(lookout2,p-h)
    dydp=(y1-y2)/h
    return(dydp)

def deriv3(lookout3,sw):
    h=0.001
    y1=lookupfr(lookout3,1,sw)
    y2=lookupfr(lookout3,1,sw-h)
    dydp=(y1-y2)/h
    return(dydp)

def initpres():
    global pg3d, sw3d, ngx, ngy, ngz, dltz, tpx, tpy, tpz, vb, owip, ogip, ooip
    dltx = grid[0][3]/grid[0][0]
    dlty = grid[0][4]/grid[0][1]
    dltz = grid[0][5]/grid[0][2]
    vb = dltx*dlty*dltz

    tpx = 6.3283*(10**(-3))*grid[0][6]*dlty*dltz/dltx
    tpy = 6.3283*(10**(-3))*grid[0][7]*dltx*dltz/dlty
    tpz = 6.3283*(10**(-3))*grid[0][8]*dlty*dltx/dltz

    Pz=[]
    Pz.append(reservoirrock[0][0])
    # print("Pz[0] = "+str(round(Pz[0],3)))
    for i in range (1,5):
        e=1
        Pa=Pz[i-1]
        Pb=Pa
        while abs(e)>0.000001:
            pmid=0.5*(Pa+Pb)
            fpb=Pa+fluidprop2(1,pmid)*dltz-Pb
            fdpb=deriv2(1,pmid)*dltz/2-1
            e=fpb/fdpb
            Pb=Pb-e
        Pz.append(Pb)
        # print("Pz["+str(i)+"] = "+str(round(Pz[i],3)))

    ngx=int(grid[0][0])
    ngy=int(grid[0][1])
    ngz=int(grid[0][2])
    pg3d = [[[0 for k in range(ngx)] for j in range(ngy)] for i in range(ngz)]
    sw3d = [[[0 for k in range(ngx)] for j in range(ngy)] for i in range(ngz)]

    sumoil = 0
    sumgas = 0
    sumwat = 0

    for i in range (ngx):
        for j in range (ngy):
            for k in range (ngz):
                pg3d[i][j][k]=Pz[k]
                sw3d[i][j][k]=reservoirrock[0][1]
                sumoil+=vb*fluidprop2(4,pg3d[i][j][k])*(1-sw3d[i][j][k])/fluidprop(1,pg3d[i][j][k])
                sumgas+=vb*fluidprop(3,pg3d[i][j][k])*fluidprop2(4,pg3d[i][j][k])*(1-sw3d[i][j][k])/fluidprop(1,pg3d[i][j][k])
                sumwat+=vb*fluidprop2(4,pg3d[i][j][k])*sw3d[i][j][k]/fluidprop(7,pg3d[i][j][k])
    ooip = round(sumoil/5.6146, 3)
    ogip = round(sumgas, 3)
    owip = round(sumwat/5.6146, 3)
    # print("OOIP (MMSTB) : " + str(round(sumoil/(5.6146*1000000),3)))
    # print("OGIP (BSCF)  : " + str(round(sumgas/1000000000,3)))
    # print("OWIP (MMSTB) : " + str(round(sumwat/(5.6146*1000000),3)))
    

def calc_rem():
    global roip, rwip
    sumoil = 0
    sumwat = 0

    for i in range(0, ngx):
        for j in range(0, ngy):
            for k in range(0, ngz):
                sumoil += vb*fluidprop2(4,pg3d[i][j][k])*(1-sw3d[i][j][k])/fluidprop(1,pg3d[i][j][k])
                sumwat += vb*fluidprop2(4,pg3d[i][j][k])*sw3d[i][j][k]/fluidprop(7,pg3d[i][j][k])
    roip = round(sumoil/5.6146, 3)
    rwip = round(sumwat/5.6146, 3)


def calc_rate(time):
    global qw, qo, dqwds, dqods, dqwdp, dqodp
    rate_var = int(wellint[7])

    qw = [0 for i in range(Nw)]
    qo = [0 for i in range(Nw)]
    dqwds = [0 for i in range(Nw)]
    dqods = [0 for i in range(Nw)]
    dqwdp = [0 for i in range(Nw)]
    dqodp = [0 for i in range(Nw)]

    for i in range (Nw):
        pwell = pg3d[wlx[i]][wly[i]][wlz[i]]
        swwell = sw3d[wlx[i]][wly[i]][wlz[i]]

        x = 0
        while time > time_rate[x][0]:
            x += 1
        if i == 0:
            # print("Untuk injector :")
            qtot = time_rate[x][1] * 5.6146

            fw = 1
        elif i == 1:
            # print("Untuk producer :")
            qtot = time_rate[x][2] * 5.6146
            mobility = (fluidprop(2,pwell) * fluidprop(1,pwell) / lookupfr(3, 1, swwell))/(fluidprop(8,pwell) * fluidprop(7,pwell) / lookupfr(2, 1, swwell))
            fw = mobility/(1+mobility)

        qw[i] = qtot * (fw)
        qo[i] = qtot * (1-fw)
        
        # print(f"Qw adalah {qw[i]} ft3/day")
        # print(f"Qo adalah {qo[i]} ft3/day")

        dqwds[i] = qw[i]*(1-fw)*((deriv3(2, swwell)/lookupfr(2,1,swwell))+(deriv3(3, swwell)/lookupfr(3,1,swwell)))
        dqods[i] = -dqwds[i]
        dqwdp[i] = qw[i]*(1-fw)*((deriv(2, pwell)/fluidprop(2, pwell))+(deriv(1, pwell)/fluidprop(1, pwell))-(deriv(8, pwell)/fluidprop(8, pwell))-(deriv(7, pwell)/fluidprop(7, pwell)))
        dqodp[i] = -dqwdp[i]
        # print(f"dQwdSw adalah {dqwds[i]}")
        # print(f"dQodSw adalah {dqods[i]}")
        # print(f"dQwdP adalah {dqwdp[i]}")
        # print(f"dQodP adalah {dqodp[i]}")
        # print()


def potent():
    global ibw, icw, idw, iew, ifw, igw, ifo, igo
    ibw = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    icw = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    idw = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    iew = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    ifw = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    igw = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    ifo = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]
    igo = [[[0 for i in range(ngx)] for j in range(ngy)] for k in range(ngz)]

    for i in range(ngx):
        for j in range(ngy):
            for k in range(ngz):
                ps = pg3d[i][j][k]
                #ibw
                if i > 0:
                    pn = pg3d[i-1][j][k]
                    potw = pn - ps
                    if potw > 0:
                        ibw[i][j][k] = 1
                #icw
                if i < ngx-1:
                    pn = pg3d[i+1][j][k]
                    potw = pn - ps
                    if potw > 0:
                        icw[i][j][k] = 1
                #idw
                if j > 0:
                    pn = pg3d[i][j-1][k]
                    potw = pn - ps
                    if potw > 0:
                        idw[i][j][k] = 1
                #iew
                if j < ngy-1:
                    pn = pg3d[i][j+1][k]
                    potw = pn - ps
                    if potw > 0:
                        iew[i][j][k] = 1
                #ifo & ifw (top)
                if k > 0:
                    pn = pg3d[i][j][k-1]
                    pm = 0.5 * (pn + ps)
                    potw = pn - ps + fluidprop2(3, pm) * dltz
                    poto = pn - ps + fluidprop2(1, pm) * dltz
                    if potw > 0:
                        ifw[i][j][k] = 1
                    if poto > 0:
                        ifo[i][j][k] = 1
                #igo & igw (bottom)
                if k < ngz-1:
                    pn = pg3d[i][j][k+1]
                    pm = 0.5 * (pn + ps)
                    potw = pn - ps - fluidprop2(3, pm) * dltz
                    poto = pn - ps - fluidprop2(1, pm) * dltz
                    if potw > 0:
                        igw[i][j][k] = 1
                    if poto > 0:
                        igo[i][j][k] = 1


def derive(sws, swn, ps, pn, dz, ijkw, ijko, pgeo):
    # print(sws, swn, ps, pn, dz, ijkw, ijko, pgeo)
    global dT, fTw, fTo
    pmid = 0.5 *(pn + ps)
    if ijkw == 1:
        krwn = lookupfr(2, 1, swn)
        dkrwn = deriv3(2, swn)
    else:
        krwn = lookupfr(2, 1, sws)
        dkrwn = 0

    if ijko == 1:
        kron = lookupfr(3, 1, swn)
        dkron = deriv3(3, swn)
    else:
        kron = lookupfr(3, 1, sws)
        dkron = 0

    if dz == 0:
        dno = 0
        dnw = 0
        ddno = 0
        ddnw = 0
    else:
        dno = fluidprop2(1, pmid) * dz
        dnw = fluidprop2(3, pmid) * dz
        ddno = deriv2(1, pmid) / 2 * dz
        ddnw = deriv2(3, pmid) / 2 * dz

    Tw = krwn / (fluidprop(8, pmid) * fluidprop(7, pmid)) * pgeo
    To = kron / (fluidprop(2, pmid) * fluidprop(1, pmid)) * pgeo

    dT = [0] * 4
    dT[0] = dkron * pgeo / (fluidprop(2, pmid) * fluidprop(1, pmid)) * (pn-ps-dno)
    dT[1] = -To / 2 * (deriv(2, pmid) / fluidprop(2, pmid) + deriv(1, pmid) / fluidprop(1, pmid)) * (pn-ps-dno) + To * (1-ddno)
    dT[2] = dkrwn * pgeo / (fluidprop(8, pmid) * fluidprop(7, pmid)) * (pn-ps-dnw)
    dT[3] = -Tw / 2 * (deriv(8, pmid) / fluidprop(8, pmid) + deriv(7, pmid) / fluidprop(7, pmid)) * (pn-ps-dnw) + Tw * (1-ddnw)
    fTo = To * (pn-ps-dno)
    fTw = Tw * (pn-ps-dnw)


def jacob():
    global ja, jb, jc, jd, je, jf, jg, Fw, Fo
    ja = [[[[0 for i in range(4)] for j in range(ngz)] for k in range(ngy)] for l in range(ngx)]
    jb = [[[[0 for i in range(4)] for j in range(ngz)] for k in range(ngy)] for l in range(ngx+1)]
    jc = [[[[0 for i in range(4)] for j in range(ngz)] for k in range(ngy)] for l in range(ngx+1)]
    jd = [[[[0 for i in range(4)] for j in range(ngz)] for k in range(ngy+1)] for l in range(ngx)]
    je = [[[[0 for i in range(4)] for j in range(ngz)] for k in range(ngy+1)] for l in range(ngx)]
    jf = [[[[0 for i in range(4)] for j in range(ngz+1)] for k in range(ngy)] for l in range(ngx)]
    jg = [[[[0 for i in range(4)] for j in range(ngz+1)] for k in range(ngy)] for l in range(ngx)]

    flux_o = [[[0 for i in range(ngz)] for j in range(ngy)] for k in range(ngx)]
    flux_w = [[[0 for i in range(ngz)] for j in range(ngy)] for k in range(ngx)]

    Fo = np.zeros((ngx, ngy, ngz), dtype=float)
    Fw = np.zeros((ngx, ngy, ngz), dtype=float)

    for i in range(ngx):
        for j in range(ngy):
            for k in range(ngz):
                ps = pg3d[i][j][k]
                sws = sw3d[i][j][k]

                #nb
                if i==0:
                    bfo = 0
                    bfw = 0
                else:
                    pn = pg3d[i-1][j][k]
                    swn = sw3d[i-1][j][k]
                    d = 0
                    derive(sws, swn, ps, pn, d, ibw[i][j][k], ibw[i][j][k], tpx)
                    bfo = fTo
                    bfw = fTw
                    for l in range(4):
                        jb[i][j][k][l] = dT[l]

                #nc
                if i==ngx-1:
                    cfo = 0
                    cfw = 0
                else:
                    pn = pg3d[i+1][j][k]
                    swn = sw3d[i+1][j][k]
                    d = 0
                    derive(sws, swn, ps, pn, d, icw[i][j][k], icw[i][j][k], tpx)
                    cfo = fTo
                    cfw = fTw
                    for l in range(4):
                        jc[i][j][k][l] = dT[l]

                #nd
                if j==0:
                    dfo = 0
                    dfw = 0
                else:
                    pn = pg3d[i][j-1][k]
                    swn = sw3d[i][j-1][k]
                    d = 0
                    derive(sws, swn, ps, pn, d, idw[i][j][k], idw[i][j][k], tpy)
                    dfo = fTo
                    dfw = fTw
                    for l in range(4):
                        jd[i][j][k][l] = dT[l]

                #ne
                if j==ngy-1:
                    efo = 0
                    efw = 0
                else:
                    pn = pg3d[i][j+1][k]
                    swn = sw3d[i][j+1][k]
                    d = 0
                    derive(sws, swn, ps, pn, d, iew[i][j][k], iew[i][j][k], tpy)
                    efo = fTo
                    efw = fTw
                    for l in range(4):
                        je[i][j][k][l] = dT[l]

                #nf
                if k==0:
                    ffo = 0
                    ffw = 0
                else:
                    pn = pg3d[i][j][k-1]
                    swn = sw3d[i][j][k-1]
                    d = -dltz
                    derive(sws, swn, ps, pn, d, ifw[i][j][k], ifo[i][j][k], tpz)
                    ffo = fTo
                    ffw = fTw
                    for l in range(4):
                        jf[i][j][k][l] = dT[l]

                #ng
                if k==ngz-1:
                    gfo = 0
                    gfw = 0
                else:
                    pn = pg3d[i][j][k+1]
                    swn = sw3d[i][j][k+1]
                    d = dltz
                    derive(sws, swn, ps, pn, d, igw[i][j][k], igo[i][j][k], tpz)
                    gfo = fTo
                    gfw = fTw
                    for l in range(4):
                        jg[i][j][k][l] = dT[l]

                flux_w[i][j][k] = bfw+cfw+dfw+efw+ffw+gfw
                flux_o[i][j][k] = bfo+cfo+dfo+efo+ffo+gfo

    acc = [0]*4

    for i in range(ngx):
        for j in range(ngy):
            for k in range(ngz):
                ps = pg3d[i][j][k]
                pn = P[i][j][k]
                sws = sw3d[i][j][k]
                swn = S[i][j][k]
                accw = -vb/dt*(fluidprop2(4,ps)*sws/fluidprop(7,ps)-fluidprop2(4,pn)*swn/fluidprop(7,pn))
                acco = -vb/dt*(fluidprop2(4,ps)*(1-sws)/fluidprop(1,ps)-fluidprop2(4,pn)*(1-swn)/fluidprop(1,pn))

                acc[0] = vb/dt*(fluidprop2(4,ps)/fluidprop(1,ps))
                acc[1] = -vb/dt*(1-sws)*((deriv2(4,ps)*fluidprop(1,ps)-fluidprop2(4,ps)*deriv(1,ps))/(fluidprop(1,ps)**2))
                acc[2] = -vb/dt*(fluidprop2(4,ps)/fluidprop(7,ps))
                acc[3] = -vb/dt*sws*((deriv2(4,ps)*fluidprop(7,ps)-fluidprop2(4,ps)*deriv(7,ps))/(fluidprop(7,ps)**2))

                sso = 0
                ssw = 0
                ss = [0]*4
                for l in range(Nw):
                    if (i==wlx[l] and j==wly[l] and k==wlz[l]):
                        sso = qo[l]
                        ssw = qw[l]
                        ss[0] = dqods[l]
                        ss[1] = dqodp[l]
                        ss[2] = dqwds[l]
                        ss[3] = dqwdp[l]

                for l in range(4):
                    ja[i][j][k][l] = -jb[i+1][j][k][l]-jc[i-1][j][k][l]-jd[i][j+1][k][l]
                    ja[i][j][k][l] = ja[i][j][k][l]-je[i][j-1][k][l]-jf[i][j][k+1][l]-jg[i][j][k-1][l]
                    ja[i][j][k][l] = ja[i][j][k][l] + acc[l]-ss[l]

                Fo[i][j][k] = flux_o[i][j][k] + acco - sso
                Fw[i][j][k] = flux_w[i][j][k] + accw - ssw


def jm_positioner():
    global jp

    # Absensi Jacobian
    jp = np.zeros((ngx, ngy, ngz, 7), dtype=int)
    count = 0
    for k in range(0, ngz):
        for j in range(0, ngy):
            for i in range(0, ngx):
                # Location A relative to i,j,k
                jp[i][j][k][0] = count
                count +=1
    # Neighbour Coordinator
    for k in range(0, ngz):
        for j in range(0, ngy):
            for i in range(0, ngx):
                # Location F relative of i,j,k
                # Up
                if(k!=0):
                    jp[i][j][k][1] = jp[i][j][k-1][0]
                else:
                    jp[i][j][k][1] = -1
                # Location D relative of i,j,k
                # Back
                if(j!=0):
                    jp[i][j][k][2] = jp[i][j-1][k][0]
                else:
                    jp[i][j][k][2] = -1
                # Location B relative of i,j,k
                # Left
                if(i!=0):
                    jp[i][j][k][3] = jp[i-1][j][k][0]
                else:
                    jp[i][j][k][3] = -1
                # Location C relative of i,j,k
                # Right
                if(i!=ngx-1):
                    jp[i][j][k][4] = jp[i+1][j][k][0]
                else:
                    jp[i][j][k][4] = -1
                # Location E relative of i,j,k
                # Front
                if(j!=ngy-1):
                    jp[i][j][k][5] = jp[i][j+1][k][0]
                else:
                    jp[i][j][k][5] = -1
                # Location G relative of i,j,k
                # Down
                if(k!=ngz-1):
                    jp[i][j][k][6] = jp[i][j][k+1][0]
                else:
                    jp[i][j][k][6] = -1
    return


def jm_constructor():
    global jm, jmm
    jm = np.zeros((ngx*ngy*ngz*2, 2*ngx*ngy*ngz), dtype=float)
    n = 0
    for k in range(0, ngz):
        for j in range(0, ngy):
            for i in range(0, ngx):
                # 2-Rows per grid
                for h in range(0, 2):   # h={0,1}
                    # 7 Derivate Members
                    for m in range(0, 7):
                        # if(jp[i][j][k][m]!=-1):
                        for mm in range(0, 2):
                            if(m==0 and jp[i][j][k][m]!=-1):   # A
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 1111
                                jm[n][jp[i][j][k][m] * 2 + mm] = ja[i][j][k][h * 2+mm]
                            elif(m==1 and jp[i][j][k][m]!=-1): # F
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 6666
                                jm[n][jp[i][j][k][m] * 2 + mm] = jf[i][j][k][h * 2+mm]
                            elif(m==2 and jp[i][j][k][m]!=-1): # D
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 4444
                                jm[n][jp[i][j][k][m] * 2 + mm] = jd[i][j][k][h * 2+mm]
                            elif(m==3 and jp[i][j][k][m]!=-1): # B
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 2222
                                jm[n][jp[i][j][k][m] * 2 + mm] = jb[i][j][k][h * 2+mm]
                            elif(m==4 and jp[i][j][k][m]!=-1): # C
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 3333
                                jm[n][jp[i][j][k][m] * 2 + mm] = jc[i][j][k][h * 2+mm]
                            elif(m==5 and jp[i][j][k][m]!=-1): # E
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 5555
                                jm[n][jp[i][j][k][m] * 2 + mm] = je[i][j][k][h * 2+mm]
                            elif(m==6 and jp[i][j][k][m]!=-1):  # G
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 7777
                                jm[n][jp[i][j][k][m] * 2 + mm] = jg[i][j][k][h * 2+mm]
                    n+=1
    nrow = 0
    jmm = np.zeros(2*ngx*ngy*ngz, dtype=float)
    for k in range(0, ngz):
        for j in range(0, ngy):
            for i in range(0, ngx):
                jmm[nrow] = -Fo[i][j][k]
                nrow+=1
                jmm[nrow] = -Fw[i][j][k]
                nrow+=1


def jm_constructor2(Ngx, Ngy, Ngz, Ja, Jb, Jc, Jd, Je, Jf, Jg, Fo, Fw):
    global jm, jmm
    jm = np.zeros((Ngx*Ngy*Ngz*2, 2*Ngx*Ngy*Ngz), dtype=float)
    n = 0
    for k in range(0, Ngz):
        for j in range(0, Ngy):
            for i in range(0, Ngx):
                # 2-Rows per grid
                for h in range(0, 2):   # h={0,1}
                    # 7 Derivate Members
                    for m in range(0, 7):
                        # if(jp[i][j][k][m]!=-1):
                        for mm in range(0, 2):
                            if(m==0 and jp[i][j][k][m]!=-1):   # A
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 1111
                                jm[n][jp[i][j][k][m] * 2 + mm] = Ja[i][j][k][h * 2+mm]
                            elif(m==1 and jp[i][j][k][m]!=-1): # F
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 6666
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jf[i][j][k][h * 2+mm]
                            elif(m==2 and jp[i][j][k][m]!=-1): # D
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 4444
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jd[i][j][k][h * 2+mm]
                            elif(m==3 and jp[i][j][k][m]!=-1): # B
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 2222
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jb[i][j][k][h * 2+mm]
                            elif(m==4 and jp[i][j][k][m]!=-1): # C
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 3333
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jc[i][j][k][h * 2+mm]
                            elif(m==5 and jp[i][j][k][m]!=-1): # E
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 5555
                                jm[n][jp[i][j][k][m] * 2 + mm] = Je[i][j][k][h * 2+mm]
                            elif(m==6 and jp[i][j][k][m]!=-1):  # G
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 7777
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jg[i][j][k][h * 2+mm]
                    n+=1
    nrow = 0
    jmm = np.zeros(2*Ngx*Ngy*Ngz, dtype=float)
    for k in range(0, Ngz):
        for j in range(0, Ngy):
            for i in range(0, Ngx):
                jmm[nrow] = -Fo[i][j][k]
                nrow+=1
                jmm[nrow] = -Fw[i][j][k]
                nrow+=1


# Main Program
aTIME = []
aDT = []
aWATINJ = []
aOILPROD = []
aWATPROD = []
aWC = []
aWOR = []
aCUMINJ = []
aCUMOPROD = []
aCUMWPROD = []
aPWBINJ = []
aPWBPROD = []
aMB_ERR_OIL = []
aMB_ERR_WAT = []

print("Subprogram:Readdata/running")
readdata()
print("Subprogram:Readdata/success")
print("")

print("Subprogram:Initpres/running")
initpres()
print("Subprogram:Initpres/success")
print("")

print("Subprogram:jm_positioner/running")
jm_positioner()
print("Subprogram:jm_positioner/success")
print("")

E_s = float(0.001)  # batasan dSw untuk dianggap konvergen
E_p = float(0.1)    # batasan dP utk dianggap konvergen
E_fo = float(1)
E_fw = float(5)

dSLIM = 0.02
dPLIM = 50

t = 0
dt = 3
tmax = 7500
# tmax = 2000
cum_oilprod = 0
cum_watprod = 0
cum_watinj = 0
cum_oilinj = 0

while t<tmax:
    P = np.zeros((ngx, ngy, ngz), dtype=float)
    S = np.zeros((ngx, ngy, ngz), dtype=float)

    t = t + dt
    for k in range(0, ngz):
        for j in range(0, ngy):
            for i in range(0, ngx):
                P[i][j][k]=pg3d[i][j][k]
                S[i][j][k]=sw3d[i][j][k]

    print("Subprogram:Potent/running")
    potent()
    print("Subprogram:Potent/success")
    print("")

    c = 0
    niter = 0
    itermax = 100
    while c==0:
        niter += 1
        print("Subprogram:Calc_rate/running")
        calc_rate(t)
        print("Subprogram:Calc_rate/success")
        print("")

        print("Subprogram:Jacob/running")
        jacob()
        print("Subprogram:Jacob/success")
        print("")

        print("Subprogram:jm_creator/running")
        # jm_constructor()
        jm_constructor2(ngx, ngy, ngz, ja, jb, jc, jd, je, jf, jg, Fo, Fw)
        print("Subprogram:jm_creator/success")
        print("")

        print("Subprogram:gauss/running")
        # sol = solve(jm, jmm)  # gauss
        lu, piv = lu_factor(jm) # LU Decomposition
        sol = lu_solve((lu, piv), jmm)
        print("time: ", t)
        print("iter: ", niter)
        print("Subprogram:gauss/success")
        print("")

        # Update Values
        # Separate Solution to Sw & P
        x_dsw = np.zeros((ngx, ngy, ngz), dtype=float)
        x_dp = np.zeros((ngx, ngy, ngz), dtype=float)
        dr = 0
        for k in range(0, ngz):
            for j in range(0, ngy):
                for i in range(0, ngx):
                    x_dsw[i][j][k] = sol[dr]
                    dr+=1
                    x_dp[i][j][k] = sol[dr]
                    dr+=1
        for k in range(0, ngz):
            for j in range(0, ngy):
                for i in range(0, ngx):
                    sw3d[i][j][k] = sw3d[i][j][k]+x_dsw[i][j][k]
                    pg3d[i][j][k] = pg3d[i][j][k]+x_dp[i][j][k]
        x_dsw_max = np.amax(abs(x_dsw))
        x_dp_max = np.amax(abs(x_dp))
        fo_max = np.amax(abs(Fo))
        fw_max = np.amax(abs(Fw))

        # if(x_dp_max<E_p and x_dsw_max<E_s):
        #     c = 1
        if(fo_max<E_fo and fw_max<E_fw and x_dp_max<E_p and x_dsw_max<E_s):
            c = 1
        else:
            if(niter>itermax):
                # t = tmax
                dt = dt*0.5
                t=t-dt
            # if dt<10**-6:
            #     t = tmax

    print("Subprogram:Calc_Rem/running")
    calc_rem()
    print("Subprogram:Calc_Rem/success")
    print("")

    for i in range(0, Nw):
        if qw[i]>0:
            Qw = qw[i]
            Qo = qo[i]
            cum_watprod += Qw*dt
            # cum_watprod += Qw*dt/5.6146
            cum_oilprod += Qo*dt
            # cum_oilprod += Qo*dt/5.6146
        if qw[i]<0:
            Qi = abs(qw[i])
            cum_watinj += abs(qw[i])*dt
            # cum_watinj += abs(qw[i])*dt/5.6146 (bbl)
            cum_oilinj += abs(qo[i])*dt
            # cum_oilinj += abs(qo[i])*dt/5.6146

    mbew = (owip-rwip-cum_watprod+cum_watinj)/owip
    mbeo = (ooip-roip-cum_oilprod+cum_oilinj)/ooip

    watcut = Qw/(Qo+Qw)
    wor = Qw/Qo

    aTIME.append(t)
    aDT.append(dt)
    aWATINJ.append(Qi)
    aOILPROD.append(Qo)
    aWATPROD.append(Qw)
    aWC.append(watcut)
    aWOR.append(wor)
    aCUMINJ.append(cum_watinj/1000)
    aCUMOPROD.append(cum_oilprod/1000)
    aCUMWPROD.append(cum_watprod/1000)
    aPWBINJ.append(pg3d[0][0][4])
    aPWBPROD.append(pg3d[4][4][4])
    aMB_ERR_OIL.append(mbeo)
    aMB_ERR_WAT.append(mbew)

    dPMAX = np.amax(abs(pg3d-P))
    dSMAX = np.amax(abs(sw3d-S))

    dtold = dt
    dT_new_p = dPLIM/dPMAX
    dT_new_s = dSLIM/dSMAX
    dt = dt*min([dT_new_s, dT_new_p])
    if(dt/dtold>2):
        dt = dtold*2
    if(dt>30):
        dt = 30
    if(t<tmax and t+dt>tmax):
        dt = tmax - t

listy = [[] for i in range(len(aTIME)+2)]
listy[0] = ["TIME","DT","WATINJ","OILPROD","WATPROD","WC","WOR","CUMINJ","CUMOPROD","CUMWPROD","PWBINJ","PWBPROD","MBEO","MBEW"]
listy[1] = ["Days","Days","SCF/D","SCF/D","SCF/D","dec.","SCF/SCF","MSCF","MSCF","MSCF","psia","psia","dec.","dec."]

for x in range(0, len(aTIME)):
    listy[x+2].append(aTIME[x])
    listy[x+2].append(aDT[x])
    listy[x+2].append(aWATINJ[x])
    listy[x+2].append(aOILPROD[x])
    listy[x+2].append(aWATPROD[x])
    listy[x+2].append(aWC[x])
    listy[x+2].append(aWOR[x])
    listy[x+2].append(aCUMINJ[x])
    listy[x+2].append(aCUMOPROD[x])
    listy[x+2].append(aCUMWPROD[x])
    listy[x+2].append(aPWBINJ[x])
    listy[x+2].append(aPWBPROD[x])
    listy[x+2].append(aMB_ERR_OIL[x])
    listy[x+2].append(aMB_ERR_WAT[x])

s = [[str(e) for e in row] for row in listy]
lens = [max(map(len, col)) for col in zip(*s)]
fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
table = [fmt.format(*row) for row in s]
with open("Output Case 5.txt", "w+") as file:
    file.write('\n'.join(table))

print("Simulation Run Completed")
