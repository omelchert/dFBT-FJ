import numpy as np

def fetchRawData(inFileName):
        rArray   = []
        zArray   = []
        ArzArray = []
        ctr = 0
        with open(inFileName,'r') as f:
            for line in f:
                    c = line.split()
                    # FILTER FOR r,z MESH PROPERTIES --------------------------
                    if len(c)>1 and c[0]!='#':
                        ctr+=1
                        if ctr == 5: 
                            dz, dr = float(c[0]), float(c[1])
                        if ctr == 6:
                            Nz, Nr = int(c[0]), int(c[1])

                    # FILTER FOR INTERNAL PHOTON DISTRIBUTION -----------------
                    if len(c)==1 and c[0]=='A_rz':
                            c = f.next().split()
                            while(len(c)>=1):
                                ArzArray.append(c)
                                c = f.next().split()

        rAxis = np.linspace(0,dr*Nr,Nr,endpoint=False)
        zAxis = np.linspace(0,dz*Nz,Nz,endpoint=False)
        ArzArray = np.asarray(ArzArray,dtype=float).reshape((Nr,Nz))

        return rAxis,zAxis,ArzArray


def writeNpz(r,z,Wrz,fName='Wrz_dataROI.npz'):
        np.savez_compressed(fName, r=r, z=z, Wrz=Wrz)


def readNpz(fName):
        npzDict = np.load(fName)
        return npzDict['r'], npzDict['z'], npzDict['Wrz']


def writeGpl(r,z,Wrz,f=None):
        f = open(f,'w') if f else sys.stdout 
        f.write("# Nr Nz %d %d\n"%(r.size, z.size))
        f.write("%d "%(r.size))
        for ir in range(r.size):
            f.write("%lf "%(r[ir]))
        f.write("\n")
        for iz in range(z.size):
            f.write("%lf "%(z[iz]))
            for ir in range(r.size):
                f.write("%lf "%(Wrz[ir,iz]))
            f.write("\n")
        f.close()

# EOF: dataIO.py
