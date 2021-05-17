import math
import numpy as np

def generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin):
    opts = varargin

    w0 = w0 * math.pi/180
    w1 = w1 * math.pi/180
    a = a * math.pi/180
    ph1 = ph1*math.pi/180
    dph0 = dph0*math.pi/180

    if len(cr) == chains - 1:
        #cr = cr.insert(0, 1)   #insert 1 at front
        cr = np.insert(cr, 0, 1, axis = 0)
    if len(dph0) == chains - 1:
        dph0 = np.insert(dph0, 0, 0, axis = 0)
    if len(zoff) == chains - 1:
        #zoff.insert(0, 0)
        zoff = np.insert(zoff, 0, 0, axis = 0)
    if cr[0] == 0:
        print("The first entry of the orientation vector can not be 0.") 
    if len(ph1) == 1:
        ph1 = [ph1]*chains
    elif len(ph1) != chains:
        print("%d helical phases specified for %d chains" % (len(ph1), chains))
    
    XYZ = np.zeros((chL * chains, 3))
    t =np.arange(0, chL)
    for i in range(0, chains):
        if cr[i] == 0:
            x, y, z = crickEQ(r0, r1, -w0, -w1, a, 0, -ph1[i], t)
            if opts == "apNNzoff":
                zoff[i] = zoff[i] + XYZ[i*chL, 3] - z[-1]
            elif opts == "registerzoff":
                zo = XYZ[0, 2] - z[-1]
                dz = absoluteToRegisterZoff(zo, r0, w0, a, w1, ph1[0], ph1[i], cr[i])
                zo = zo - dz
                zoff[i] = zoff[i] + zo
            elif opts == "zoffaa":
                zo = XYZ[0, 3] - z[-1]
                dz = absoluteToZoff_aa(zo, r0, w0, a, r1, w1, ph1[0], ph1[i], cr[i])
                zo = zo - dz
                zoff[i] = zoff[i] + zo
            T = np.array([ 
                 [math.cos(dph0[i] - zoff[i]*math.tan(a)/r0) , math.sin(dph0[i] - zoff[i] * math.tan(a)/r0) , 0],
                 [-math.sin(dph0[i] - zoff[i]*math.tan(a)/r0), math.cos(dph0[i] - zoff[i] * math.tan(a)/r0) , 0],
                 [0, 0, 1]
            ])
        else:
            x, y, z = crickEQ(r0, r1, w0, w1, a, 0, ph1[i], t)
            if opts == "registerzoff":
                zo = 0
                dz = absoluteToRegisterZoff(zo, r0, w0, a, w1, ph1[0], ph1[i], cr[i])
                zo = zo - dz
                zoff[i] = zoff[i] + zo
            elif opts == "zoffaa":
                zo = 0
                dz = absoluteToZoff_aa(zo, r0, w0, a, r1, w1, ph1[0], ph1[i], cr[i])
                zo = zo - dz
                zoff[i] = zoff[i] + zo
            T = np.array([
                 [math.cos(dph0[i] - zoff[i]*math.tan(a)/r0), math.sin(dph0[i] - zoff[i]*math.tan(a)/r0), 0],
                 [-math.sin(dph0[i] - zoff[i]*math.tan(a)/r0), math.cos(dph0[i] - zoff[i] * math.tan(a)/r0), 0],
                 [0, 0, 1 ] 
            ])  #please check *np.ones(1)

        # if len(r0) == chL:
        #     xyz = [x, y, z]
        #     for k in range(0, len(r0)):
        #         xyz[:, k] = T[:, k::chL]* xyz[:, k] 
        # else:
        #    xyz = T*[x, y, z]
        xyz = np.matmul(T, np.array([x, y, z])) 

        XYZ[i*chL : (i+1)*chL, 0] = (xyz[0, :]).conj().T 
        XYZ[i*chL : (i+1)*chL, 1] = (xyz[1, :]).conj().T 
        XYZ[i*chL : (i+1)*chL, 2] = (xyz[2, :]).conj().T + zoff[i]

    return XYZ



def absoluteToRegisterZoff(zoff, R0, w0, a, w1, ph1_1, ph1_2, p_ap):
    if p_ap != 1 and p_ap != 0:
        print("p_ap flag expected to be either 1 or -1")
    
    aa1 = 2*math.pi/w1
    b1 = (math.pi - ph1_1)/w1
    z1 = (R0*w0/math.tan(a))*(aa1*1 + b1)

    if p_ap == 0:
        aa2 = 2*math.pi/(-w1)
        b2 = (math.pi + ph1_2)/(-w1)
        n = ((z1 - zoff)/(-R0*w0/math.tan(a)) - b2)/aa2
        dz = -(R0*w0/math.tan(a))*(aa2*math.floor(n) + b2) + zoff - z1
        dz1 = -(R0*w0/math.tan(a))*(aa2*math.ceil(n) + b2) + zoff - z1
    else:
        b2 = (math.pi - ph1_2)/w1
        n = ((z1-zoff)/(R0*w0/math.tan(a)) - b2)/aa1
        dz = -(R0*w0/math.tan(a))*(aa1*math.floor(n) + b2) + zoff -z1
        dz1 = (R0*w0/math.tan(a))*(aa1*math.ceil(n) + b2) + zoff - z1
    
    if abs(dz1) < abs(dz):
        rzoff = dz1
    else:
        rzoff = dz
    
    return rzoff

def absoluteToZoff_aa(zoff, R0, w0, a, R1, w1, ph1_1, ph1_2, p_ap):
    if p_ap != 1 and p_ap != 0:
        print("p_ap flag expected to be either 1 or -1")

    rng = range(0, 7)
    

def angleDiff(a, b):
    d = (a % (2*math.pi) - b % (2*math.pi)) % 2*math.pi
    _in = [i for i,v in enumerate(d) if v > math.pi]
    for i in _in:
        d[i] = d[i] - 2*math.pi
    return d

def canonicalPhases(ind):
    # canonical phases are:
    # [41.0, 95.0, 146.0, 197.0, 249.0, 300.0, 351.0]
    # corresponding to {'c', 'g', 'd', 'a', 'e', 'b', 'f'}, respectively
    median_phases = [197.0, 300.0, 41.0, 146.0, 249.0, 351.0, 95.0]
    return median_phases[ind]

def crickEQ(r0, r1, w0, w1, a, ph0, ph1, t):
    x = r0*np.cos(w0*t + ph0) + r1*np.cos(w0*t + ph0)*np.cos(w1*t + ph1) - r1*np.cos(a)*np.sin(w0*t + ph0)*np.sin(w1*t + ph1)
    y = r0*np.sin(w0*t + ph0) + r1*np.sin(w0*t + ph0)*np.cos(w1*t + ph1) + r1*np.cos(a)*np.cos(w0*t + ph0)*np.sin(w1*t + ph1)
    z = w0*t*r0/np.tan(a) - r1*np.sin(a)*np.sin(w1*t + ph1)
    return [x, y, z]

def getHeptadPos(ph1, varargin):
    meds = canonicalPhases(slice(0, 7))*pi/180
    hps = ['a', 'b', 'c', 'd', 'e', 'f', 'g']

    si = np.argsort(meds).tolist()
    meds = [meds[i] for i in si]
    hps = [hps[i] for i in si]

    ph1 = ph1%(2*math.pi)
    for i in range(1, len(meds)+1):
        pin = i - 1
        if i == 1:
            pin = len(meds)
        nin = i + 1
        if i == len(meds):
            nin = 1
        lb = (angleDiff(meds[pin - 1], meds[i-1])/2 + meds[i-1] ) % (2*math.pi)
        ub = (angleDiff(meds[nin - 1], meds[i-1])/2 + meds[i-1] ) % (2*math.pi)
        
        if angleDiff(ph1, lb) > 0 and angleDiff(ub, ph1) > 0:
            if varargin == 1:
                hp = i - 1    
            else:
                hp = hps[i-1]
    return hp

def default_pdb_ca_coords():
    chains = 4
    chL = 28
    r0 = 7.36
    r1 = 2.26
    w0 = -2.45
    w1 = 102.857
    a = -12.01
    ph1 = np.array([-9, -9, -9, -9])
    cr = np.array([0, 1, 0])
    zoff = np.array([0.0, 0.0, 0.0])
    dph0 = np.array([90, 180, 270])
    varargin = 'registerzoff'

    XYZ = generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)

    np.savetxt('XYZ_py.txt', XYZ, delimiter='\t')