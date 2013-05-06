# normalize polygons in a directory to unit area

import math
import os
import shutil
from numpy import *

def polygonArea( pgnVertList ):
    c = array([0.0, 0.0, 0.0]) #polygon geometric center
    for i in pgnVertList:
        c += i
    c /= len(pgnVertList)
    area = 0.0
    for i in range(0, len(pgnVertList)):
        tri_edge_0 = pgnVertList[i] - c
        tri_edge_1 = pgnVertList[(i+1)%len(pgnVertList)] - c
        crss = cross(tri_edge_0, tri_edge_1)
        area += math.sqrt(dot(crss, crss))
    return 0.5*area

inDir = r"D:\research\Tiling\prog\newspm\data\polygons\bwg_stones"
outDir = r"D:\research\Tiling\prog\newspm\data\polygons\bwg_stones_kitten\\"
os.chdir(inDir)
filenames = os.listdir(inDir)

for f in filenames:
    ext = os.path.splitext(f)
    if ( ext[1] == ".tply"):
        ifh = open(f, "r")
        ofh = open(outDir+f, "w")
        nbVert = int(ifh.readline())
        pgn_verts = []
        text_coords = []
        for i in range(0, nbVert):
            aline = ifh.readline().split()
            x = float(aline[0])
            y = float(aline[1])
            tx = float(aline[2])
            ty = float(aline[3])
            pgn_verts.append(array([x, y, 0.0]))
            text_coords.append((tx, ty))
        area = polygonArea(pgn_verts)
        scale = math.sqrt(1.0/area)
        ofh.write(str(nbVert)+"\n")
        for j in range(0, nbVert):
            scaled_vert = scale*pgn_verts[j]
            ofh.write(str(scaled_vert[0]) + " " + str(scaled_vert[1]) + " ")
            ofh.write(str(text_coords[j][0]) + " " + str(text_coords[j][1]) + "\n")
        ifh.close()
        ofh.close()
    else:
        shutil.copyfile(f, outDir+f)
    
