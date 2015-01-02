from math import *
import os
import shutil
#adjust multiple polygon files in a directory

inDir = r"D:\research\Tiling\prog\newspm\data\polygons\fx_size\green_stone_egg"
outDir = r"D:\research\Tiling\prog\newspm\data\polygons\temp\\"
#scale = sqrt(130000.0*0.75)
scale = sqrt(29.2479/3000/2.39253) #302: 100 tiles
os.chdir(inDir)

filenames = os.listdir(inDir)

for f in filenames:
    ext = os.path.splitext(f)
    if ( ext[1] == ".tply"):
        ifh = open(f, "r")
        ofh = open(outDir+f, "w")
        #ofh = open(ext[0]+ext[1], "w") #double original files
        nbVert = int(ifh.readline())
        if (nbVert == 2):
            print "invalid polygon indexed at"
            print f
        ofh.write(str(nbVert)+"\n")
        #coords = []
        for i in range(0, nbVert):
            aline = ifh.readline().split()
            #print aline
            x = float(aline[0])
            y = float(aline[1])
            ofh.write(str(scale*x)+" "+str(scale*y)+" "+aline[2]+" "+aline[3]+"\n")
            #coords.append(str(scale*x)+" "+str(scale*y)+" "+aline[2]+" "+aline[3]+"\n")
        #reverse its orientation
        #for i in range(len(coords)-1, -1, -1):
        #    ofh.write(coords[i])
        ifh.close()
        ofh.close()
    #else:
    #    shutil.copyfile(f, outDir+f)
        #shutil.copyfile(f, ext[0]+ext[1])
      
