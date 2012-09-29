from math import *

N = 8
R1 = 10.0 ;

for i in range(N):
   x =  R1 * sin(i * 2 * pi / N) 
   y = -R1 * cos(i * 2 * pi / N)
   print "v " + str(x) + " " + str(y) ;


for i in range(N):
   print "e " + str(i+1) + " " + str(((i + 1) % N)+1)


