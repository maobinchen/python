#! /usr/bin/env python
'this script use generate lines as described in Laster paper to build hyperboloid of one sheet'
from math import sqrt,asin,sin,tan,cos,atan,pi,floor
from scipy.optimize import fsolve,fmin
import argparse
import string
import sys
parser = argparse.ArgumentParser(description='Program')
parser.add_argument('n', action="store",type=int, help="Number of beta strands of beta-barrel")
parser.add_argument('S', action="store",type=int,help="Shear number of beta-barrel")
parser.add_argument('nres', action="store",type=int,help="Number of residues in each strand")
parser.add_argument('dtw', action="store",type=float,default=0,help="twist angle difference relative to the opitmal")
parser.add_argument('b', action="store",type=float,default=4.4,help="inter strand distance")
parser.add_argument('da', action="store",type=float,default=0,help="radius difference relative to the opitmal")
parser.add_argument('db', action="store",type=float,default=0,help="radius difference relative to the opitmal")
parser.add_argument('topology', action="store",default='1,2,5,4,3,6',help="toplogy")
parser.add_argument('--a', action="store_true",help="Anti-parralel")

a=3.3
#topology = [1,2,5,4,3,6] #Greek key
PI = pi
args = parser.parse_args()
topology_str = args.topology #Greek key
topology = topology_str.split(',')
topology = [int(c) for c in topology]
n=args.n
S=args.S
nres = int(args.nres)
b=args.b
da=args.da 
db=args.db
dtw=args.dtw

theta=atan(S*a/(n*b)) #the tilt angle of beta strand
def twistAndCoil(x):
    return (x-20*PI/180)**2+((2*PI*sin(theta)-n*x)/S)**2+((2*PI*cos(theta)-S*a*x/b)/n)**2
tw1 = fmin(twistAndCoil,0)[0] #twist angle between strand
nres += 4
tw1 *= dtw
C = sqrt((b**2*(1+cos(tw1))/(4*(1-cos(tw1)))))
r = sqrt(b**2*(1+cos(tw1))/(4*(cos(tw1)-cos(2*PI/n))))
A=r*da
B=r*db

print "A=%f,B=%f,C=%f,twist=%f,tilt=%f\n" %(A,B,C,tw1*180/PI,theta*180/PI)




def vecLength(vec):
    sum = 0
    for x in vec:
	sum += x**2
    return sqrt(sum) 

coor = []
ref = [] #reference positions(in the middle of the barrel)
sign0 = 1 #the sign of middle residue
shifts = [0,]
for i in range(1,n+1):
    omega = 2*PI*(i-1)/n 
    x0,y0,z0 = A*cos(omega),B*sin(omega),0
    #direction vector of two skewer lines family
    v = (-A*sin(omega),B*cos(omega),C)
    dt = a/vecLength(v) #since (xt,yt,zt) = (x0,y0,z0)+vt
    if i>1:
	def pos(t):
	    xt,yt,zt = x0+v[0]*t,y0+v[1]*t,z0+v[2]*t #coordinate of points on strand i
	    vec_refpt = (xt-ref[i-2][0],yt-ref[i-2][1],zt-ref[i-2][2]) #the vector of previous reference point to the current vector
	    out = vec_refpt[0]*v[0]+vec_refpt[1]*v[1]+vec_refpt[2]*v[2] #inner product of vector p(i,t0)-p(i-1,t0) and direction vector
	    return out
	ti = fsolve(pos,-0.5)[0] #get the value of parameter t which corresponding to the ca coordinate of strand i+1 in register with reference ca of strand i 
	k0 = 0
	for k in range(-10,10):
	    tk = ti+dt*k
	    if abs(tk)<abs(ti):
		k0 = k
		sign0 *= (-1)**k #the sign when t=0, j=(nres-1)/2
	ti = ti+dt*k0
	shifts.append(k0)
	x0,y0,z0 = x0+v[0]*ti,y0+v[1]*ti,z0+v[2]*ti
	
    ref.append([x0,y0,z0])
    coor_strand = []
    for j in range(nres):
	resi = int(floor((j-(nres-1)/2.0)))
	t = resi*dt  
	xt,yt,zt = x0+v[0]*t,y0+v[1]*t,z0+v[2]*t
	sign = sign0*((-1)**abs(resi))
	vn = 2*xt/(A**2),2*yt/(B**2),-2*zt/(C**2) #the normal vector of point at the surface
	dtn = sign*sqrt(1.9*1.9-1.65*1.65)/vecLength(vn) 
	xn,yn,zn = xt+vn[0]*dtn,yt+vn[1]*dtn,zt+vn[2]*dtn
	dd = xn-xt,yn-yt,zn-zt
        coor_strand.append([xn,yn,zn]) 
    coor.append(coor_strand)
    

outfile = 'test_ca.pdb'
chain_ids = string.uppercase
fo = open(outfile,'w')
resi = 1
atomi = 2
for i in range(len(coor)):
    coor_strand = coor[topology[i]-1]
    if args.a and i%2 == 1: coor_strand.reverse()
    for j in range(len(coor_strand)):
        (x,y,z) = coor_strand[j]
        fo.write("ATOM  %5d  CA  ALA %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (atomi,chain_ids[i],resi,x,y,z))
        resi += 1
        atomi += 4

fo.close()
###adding meaningless comments
