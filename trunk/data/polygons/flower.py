from math import *

N = 100
R1 = 10.0 ;
R2 = 3.0 ;

for i in range(N):
   alpha = i * 2 * pi / N ;
   R = R1 + R2 * sin(5*alpha)
   x =  R * sin(alpha)
   y = -R * cos(alpha)
   print "v " + str(x) + " " + str(y) ;


for i in range(N):
   print "e " + str(i+1) + " " + str(((i + 1) % N)+1)


