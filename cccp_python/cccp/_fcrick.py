import math
import numpy as np
from _readCA import readCA as rc

def sign(x):
    if x >0:
        hand = 1
    elif x <0:
        hand = -1
    else:
        hand =0
    return hand

#def frcick(file, chains, parType, coorType, outType, iIP, iLB, iUB, imask):

file = '/mnt/e/GitHub_Design/cccp/XYZ.txt'
chains = 4
parType = 'GENERAL'
coorType = 0
outType = 'out'
iIP = []
iLB = []
iUB = []
imask = []

#check input params
paralist = ['GENERAL', 'SYMMETRIC', 'ZOFF-SYMM', 'DPH0-SYMM', 'GENERAL-HLXPH', 'SYMMETRIC-HLXPH', 'ZOFF-SYMM-HLXPH', 'DPH0-SYMM-HLXPH']
if not any(parType == s for s in paralist):
    print('Unknown parameterization type ', parType)

if not iIP:
    IP = [5.0, 2.26, -2*math.pi/100, 4*math.pi/7, -12*math.pi/180, 0]
elif len(iIP[0])*len(iIP[1]) != 6:
    print('Unexpected number of initial parameters!')
else:
    IP = iIP

M = rc(file, coorType)

n = len(M[:, 0])
if n%chains != 0:
    print('number of coordinates ', n, ' is not divisible by the number of chains ', chains)
olig = chains
sym_type = 'GENERAL'
shw = 0
mask = imask


### determine chain order (clock-wise)
#if olig > 2:

nc = int(n/olig)
mask_zero = [i for i, e in enumerate(mask) if e == 0]
ind_init = np.setdiff1d(list(range(1, nc+1)), mask_zero)
c = [M[round(len(ind_init)/2)-1, :]]
pz = M[ind_init[-1]-1, :] - M[ind_init[0]-1, :]
for i in range(1, olig):
    ind = np.setdiff1d(list(range(i*nc+1, (i+1)*nc+1)) , mask_zero)

    pow_arr = np.power((M[ind-1, :] - c[0]), 2)
    sum_arr = np.sum((pow_arr), axis=1)
    mi = np.argmin(sum_arr)
    #m = sum_arr[mi]
    c = np.vstack([c, M[ind[0] + mi -1, :]])
cc = np.mean(c, axis=0)
an = [0]
for i in range(1, olig):
    rot_an= np.dot((cc - c[0, :]), (cc - c[i, :]))/np.linalg.norm(c[0, :]- cc)/np.linalg.norm(cc - c[i, :])          
    an.append(math.acos(rot_an))
    if np.dot(np.cross((cc - c[0, :]), (cc - c[i, :])) ,pz) < 0:
        an[i] = 2*math.pi - an[i]
co = np.argsort(an).tolist()
an = [an[i] for i in co]

#else:
#    co = list(range(1:olig+1))

### determine chain orientation and handedness (sign of crossing angle)
#if  olig >= 2:
nc = int(n/olig)
mask_zero = [i for i, e in enumerate(mask) if e == 0]
ind_init = np.setdiff1d(list(range(1, nc+1)), mask_zero)
pz = M[ind_init[-1]-1, :] - M[ind_init[0]-1, :]
crsg = [0]*(olig-1)

cr = [1]
for i in range(1, olig):
    ind = np.setdiff1d(list(range(i*nc+1, (i+1)*nc+1)) , mask_zero)
    pz_curr = M[ind[-1]-1, :] - M[ind[0]-1, :]
    cr.append(1 if np.dot(pz_curr, pz) > 0 else 0)
    crsg[i-1] = crossingAngle(M[ind_init-1, :], M[ind-1, :], cr[i])

if len(iIP) == 0:
    hand = sign(np.mean(crsg))
    if hand==0:
        hand = 1
    IP[2] = abs(IP[2])*hand 
    IP[4] = abs(IP[4])*hand
# else:
#     cr = [1]

ap = [i for i in range(len(cr)) if cr[i] == 0]
nap = len(ap)



### Minimization
p0 =[]
#[val, name, LB, UB, pri]
p0.append([IP[0], 'R0 (A)', 0, 30, 1])
p0.append([IP[1], 'R1 (A)', 2, 3, 1])
p0.append([IP[2], 'w0 (rad/res)', -10*np.pi/100, 10*np.pi/100, 1])
p0.append([IP[3], 'w1 (rad/res)', 1*np.pi/100, 8*np.pi/100, 1])
p0.append([IP[4], 'alpha (rad)', -np.pi/6, np.pi/6, 1])

if not 'HLXPH' in parType:
    p0.append([IP[5], 'ph1 (rad)', -np.pi, np.pi, 1])
else:
    for i in range(olig):
        p0.append([IP[5], ('ph1 for chain %d (rad)'%(i+1)), -np.pi, np.pi, 1])

if 'SYMMETRIC' in parType or 'ZOFF-SYMM' in parType:
    if len(ap)==0:
        p0.append([max(M[:, 2])-min(M[:, 2]), 'absolute ap zoff (A)', -max(M[:, 2])+min(M[:, 2]), max(M[:, 2])-min(M[:, 2]), 4])
else:
    for i in range(1, olig):
        if cr[i]:
            p0.append([0, ('absolute zoff_%d (A)' % i), -max(M[:, 2])+min(M[:, 2]), max(M[:, 2])-min(M[:, 2]), 4])
        else:
            p0.append([max(M[:, 2])-min(M[:, 2]), ('absolute zoff_%d (A)'%(i+1)), -max(M[:, 2])+min(M[:, 2]), max(M[:, 2])-min(M[:, 2]), 4])

if not 'SYMMETRIC' in parType and not 'DPHO-SYMM' in parType:
    for i in range(1, olig):
        p0.append([-(co[i])*2*np.pi/olig, ('dph0_%d (rad)'%(i+1)), -2*np.pi, 0, 1])
elif 'SYMMETRIC' in parType and all([1 for c in co[::2] if cr[c]== 1]) and all([1 for c in co[1::2] if cr[c]== 0]):
    p0.append([-2*np.pi/olig, 'dph0_p_ap (rad)', -2*np.pi, 0, 1])

err = np.Infinity
x0 = [p[0] for p in p0]
p0_name = [p[1] for p in p0]

#for it in range(100):
it = 0
cparType = 'single%s'%parType
ordi = [ip for ip in range(len(p0_name)) if 'zoff' in p0_name[ip]]
ordi.extend([ip for ip in range(len(p0_name)) if 'dph0' in p0_name[ip]])
ordi.extend([ip for ip in range(len(p0_name)) if 'ph1' in p0_name[ip]])
_r1_w1 = [ip for ip in range(len(p0_name)) if 'R1' in p0_name[ip]]
_r1_w1.extend([ip for ip in range(len(p0_name)) if 'w1' in p0_name[ip]])
ordi_c = ordi.copy()
ordi_c.extend(_r1_w1)
ordi.extend( list(set(range(len(x0))) - set(ordi_c)))



for i in range(len(ordi)):
    extr_vary = ordi[i]



#def crickssd(p,  M, x0, olig, co, cr, cparType, shw, extr_vary, mask):
    '''
    % Returns the error and its gradient between the coordinate set and the
    % ideal backbone given Crick parameters
    '''
i = 0
extr_vary = ordi[i]
p = x0[ordi[i]]
if 'SINGLE' in str.upper(cparType):
    ptemp = p 
    p = x0
    p[extr_vary] = ptemp

if 'GENERAL' in str.upper(cparType):
    r0, p = shift(p)
    r1, p = shift(p)
    w0, p = shift(p)
    w1, p = shift(p)
    a, p = shift(p)
    if 'HLXPH' in str.upper(cparType):
        ph1, p = shift(p, [olig])
    else:
        ph1, p = shift(p)
    zoff, p = shift(p, [olig -1])
    dph0, p = shift(p, [olig -1])
    
if 'SYMMETRIC' in str.upper(cparType):
    r0, p = shift(p)
    r1, p = shift(p)
    w0, p = shift(p)
    w1, p = shift(p)
    a, p = shift(p)
    if 'HLXPH' in str.upper(cparType):
        ph1, p = shift(p, [olig])
    else:
        ph1, p = shift(p)
    #TO DO: not completed

if 'ZOFF-SYMM' in str.upper(cparType):
    r0, p = shift(p)
    r1, p = shift(p)
    w0, p = shift(p)
    w1, p = shift(p)
    a, p = shift(p)
    if 'HLXPH' in str.upper(cparType):
        ph1, p = shift(p, [olig])
    else:
        ph1, p = shift(p)
    #TO DO: not completed

if 'DPH0-SYMM' in str.upper(cparType):
    r0, p = shift(p)
    r1, p = shift(p)
    w0, p = shift(p)
    w1, p = shift(p)
    a, p = shift(p)
    if 'HLXPH' in str.upper(cparType):
        ph1, p = shift(p, [olig])
    else:
        ph1, p = shift(p)
    #TO DO: not completed

a = abs(a)*sign(w0)

if 'GENERAL' in str.upper(sym_type):
    n = M.shape[0]
    nc = int(n/olig)
    XYZ = generateCrickBBrad()
else if 'NONE' in str.upper(sym_type):
    print('TO DO.')
else if 'C2' in str.upper(sym_type):
    print('TO DO.')
else if 'C3' in str.upper(sym_type):
    print('TO DO.')
else if 'C' in str.upper(sym_type):
    print('TO DO.')
else if 'GENERAL.old' in str.upper(sym_type):
    print('TO DO.')

def generateCrickBBrad(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin = []):
    if len(varargin):
        opts = varargin[0]
    else
        opts = []

    if len(cr) == chains-1:
        cr = [1] +  cr 
    if len(dph0) == chains-1:
        dph0 = [0] + dph0
    if len(zoff) == chains -1:
        zoff = [0] + zoff
    if cr[0] ==0:
        print('The first entry of the orientation vector can not be 0.')
    if len(ph1) ==1:
        ph1 = ph1*np.ones(chains)
    elif len(ph1)!=chains:
        print('helical phases specified')
    
    XYZ = np.zeros((chL*chains, 3))
    t = list(range(chL))
    for i in range(chains):
        if cr[0] == 0:
            x, y, z = crickEQ(r0, r1, -w0, -w1, a, 0, -ph1[i], t)
        if 'apNNzoff' in opts:
            zoff[i] = zoff[i] + XYZ[i*chL, 2] - z[-1]
        elif 'registerzoff' in opts:
            zo = XYZ[0, 2] - z[-1]
            dz = absoluteToRegisterZoff(zo, r0, w0, a, w1, ph1[0], ph1[i], cr[i])
            zo = zo-dz
            zoff[i] = zoff[i] + zo
        elif 'zoffaa' in opts:
            zo = XYZ[0, 2] - z[-1]
            dz = absoluteToZoff_aa(zo, r0, w0, a, w1, ph1[0], ph1[i], cr[i])
            zo = zo-dz
            zoff[i] = zoff[i] + zo
        




def shift(x, varargin = []):
    n = 1
    if len(varargin) > 0:
        n = varargin[0]
    if n > len(x):
        print('vector does not have %d values'%n)
    v = x[:n] if n>1 else x[0]
    x = x[n:]
    return v, x

def crossingAngle(A, B, pap):
    '''
    % Computes the crossing angle between two chains of any size. First looks
    % for a pair of segment of 8 residues, one on each chain, such that they
    % are closest to each other. Then computes the crossing angle for these two
    % segments
    '''
    if  len(A[0, :]) != 3 or len(B[0, :]) != 3:
        print('Unexpected matrix size.')
    if pap == 0:
        B = np.flip(B, 0)
    if min(A.shape[0], B.shape[0]) < 3:
        a = 0
        return
    axsA = helicalAxisPoints(A)
    axsB = helicalAxisPoints(B)

    a = -dihe(axsA[0], axsA[-1], axsB[-1], axsB[0])
    # p1 = axsA[0]
    # p2 = axsA[-1]
    # p3 = axsB[-1]
    # p4 = axsB[0]

    return a

def helicalAxisPoints(H):
    axs = []
    if H.shape[1] < 3:
        return
    for i in range(1, H.shape[0] - 1):
        r = (H[i-1, :] - H[i, :]) + (H[i+1, :] - H[i, :])
        r = 2.26*r/np.linalg.norm(r)    
        axs.append(2.26*r/np.linalg.norm(r) + H[i, :])    
    return axs
    

def dihe(p1, p2, p3, p4):
    '''
    % DIHE calculates the dihedral angle given four points
    % DIHE(p1, p2, p3, p4) returns the dihedral angle p1p2p3p4.
    '''
    dim = 1
    v12 = p1 - p2
    v23 = p2 - p3
    v43 = p4 - p3

    px1 = np.cross(v12, v23)
    px1 = px1/(np.ones(dim)*math.sqrt(sum(px1*px1)))

    px2 = np.cross(v43, v23)
    px2 = px2/(np.ones(dim)*math.sqrt(sum(px2*px2)))

    dp12 = np.dot(px1, px2)
    sin2 = 1 - dp12*dp12

    d = np.pi/2.0 - math.atan(dp12/math.sqrt(sin2))

    px3 = np.cross(px1, px2)
    ind = np.dot(px3, v23) > 0
    if ind:
        #d[ind] = -d[ind]
        d = -d
    return d