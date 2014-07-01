from math import *
import os
import shutil
#adjust multiple polygon files in a directory

inDir = r"D:\research\Tiling\prog\newspm\data\polygons\kitten_cartoon"
outDir = r"D:\research\Tiling\prog\newspm\data\polygons\fx_size\cartoon_kitten\\"
#scale = sqrt(130000.0*0.75)
scale = sqrt(0.5) #302: 100 tiles
os.chdir(inDir)

filenames = os.listdir(inDir)

extend_len = 0

for f in filenames:
    ext = os.path.splitext(f)
    #idx = int(ext[0][2:5])
    if ( ext[1] == ".tply"):
        ifh = open(f, "r")
        ofh = open(outDir+ext[0]+"_small"+ext[1], "w")
        #ofh = open(outDir+str(idx+extend_len)+ext[1], "w")
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
    else:
        shutil.copyfile(f, outDir+ext[0]+"_small"+ext[1])
        #shutil.copyfile(f, outDir+str(idx+extend_len)+ext[1])
        #shutil.copyfile(f, ext[0]+ext[1])
      
