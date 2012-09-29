from math import *

N = 200
R1 = 10.0 ;

for i in range(N):
   R = R1 + 3 * sin(5 * i * 2 * pi / N)
   x =  R * sin(i * 2 * pi / N) 
   y = -R * cos(i * 2 * pi / N)
   print "v " + str(x) + " " + str(y) 


for i in range(N):
   print "e " + str(i+1) + " " + str(((i + 1) % N)+1)


