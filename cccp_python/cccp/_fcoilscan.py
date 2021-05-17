from os
import readCA as rc
import frick as fk

#This is for accommodation index calculation. Impletement it in the future.
#The function is not complete.
def fcoilscan(file, coorType, Nchains, win, Ph_ideal, outDir, xlOut = 0):

    M = rc.readCA(file)

    if not os.path.exists(outDir):
        os.mkdir(outDir)

    if len(M) % Nchains != 0:
        sys.exit('Total number of bundle coords is not divisible by number of chains!')

    N = len(M) / Nchains
    output = []

    for a in range(1, N + 2 - win):
        print('\nfitting window %d/%d :REPORT', a, (N - (win-1))))
        Calphas = []
        for m in range(1, Nchains + 1):
            Calphas.append(M[N*(m - 1) + a:N*(m - 1) + a + win - 1][:])
        
        err, XYZ, pret = fcrick(Calphas, Nchains, 'GENERAL-HLXPH', 2, 0, [], [], [], [])
        outpt = pret  #how to format pret
        output.append((outpt[1], err)) #how to format err

    print('\nextracting accommodation index data...:REPORT')


