#postprocess.py
# takes a file and turns it into clusters and plots 
#usage :
#python postprocess.py [filename] [numberOfFrames] [startFrame]
#Jim Lakowski 7/8/2013
#the way this works...
#3d arrays arrays (array[a])
#a is either 0,1,2  
#the zeroth index is always x
#the first always is y
#the third is always z
#in the 6d arrays, 
#the first array[m][k]
#m is the dancer
#
#TODO



import scipy.cluster.hierarchy as hcluster
import numpy as np
import fileinput
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import os
import math
import time
import sys
#helper functions 
def movingaverage(a, n) :
    #using a moving average algorithm 
    #https://en.wikipedia.org/wiki/Moving_average
    ret = np.cumsum(a, dtype=float)
    return (ret[n - 1:] - ret[:1 - n]) / n

def matrixToAngle(mat):
    #using a formula
    #http://en.wikipedia.org/wiki/Axis-angle_representation
    a = math.acos((np.trace(mat)-1)/2)
    return a

def matrixToAxis(mat, angle):
    #same link as above
    axis = 1/(2*math.sin(angle)) * np.array([mat[2][1] - mat[1][2], mat[0][2] - mat[2][0], mat[1][0]-mat[0][1]])
    return axis

def findLargest(vector):
    #finds the largest number in a 1D array
    largest = 0
    for i in range(len(vector)):
        if(vector[i]>largest):
            largest = vector[i]

    return largest

def countInstances(vector):
#finds the number of points in each cluster
    quantities = np.array([])
    numclus = findlargest(vector)
    for k in range(numclus):
        quantities = np.append(quantities, [k,0])

    for i in range(len(vector)):
        for j in range(len(quantities)):
            if(quantities[j][0] == vector[i]):
                quantities[j][1] = quantities[j][1] + 1
    return quantities                           


#array initiallization
xdot = np.array([])
ydot = np.array([])
zdot = np.array([])
x1dot= np.array([])
x1= np.array([])
x2dot= np.array([])
y1dot= np.array([])
y2dot= np.array([])
z1dot= np.array([])
z2dot= np.array([])
speedt = np.array([])
comt = [0,0,0]#total center of mass
comtp = [0,0,0]#previous center of mass
dcom = [[0,0,0],[0,0,0]]
com = [[0,0,0],[0,0,0]]
J1 = np.array([[0,0,0],[0,0,0],[0,0,0]])#Inerita Tensors!
J2 = np.array([[0,0,0],[0,0,0],[0,0,0]])
jar = [J1,J2]
comx= np.array([])
comy= np.array([])
comz= np.array([])
spcom1 = np.array([])#for plotting
spcom2 = np.array([])
J1xdot = np.array([])
J1ydot = np.array([])
J1zdot = np.array([])
J2xdot = np.array([])
J2ydot = np.array([])
J2zdot = np.array([])
JCOM = np.array([[0,0,0],[0,0,0],[0,0,0]])
dthetat = np.array([])
R1 = np.array([[0,0,0],[0,0,0],[0,0,0]])
R2 = np.array([[0,0,0],[0,0,0],[0,0,0]])
rot = [R1,R2]
dR1 = np.array([[0,0,0],[0,0,0],[0,0,0]])
dR2 = np.array([[0,0,0],[0,0,0],[0,0,0]])
rotcom = np.array([[0,0,0],[0,0,0],[0,0,0]])
totalrot = np.array([0,0,0])
timespread = np.array([])
axisx = np.array([])
axisy = np.array([])
axisz = np.array([])
numclusters = np.array([])
c1c2 = 550
#user interface
for i in range(len(sys.argv)):
    if(sys.argv[i] == 'help' or sys.argv[i] == '--help'):
        print "to postprocess a file enter 'python clusfile.py filename numberofframes startingframe'"
        sys.exit()

startframe = 0 #default is zero

fname = sys.argv[1]
totlines =  sum(1 for line in open(fname))    
numfr = totlines-7

if len(sys.argv) > 2:
    startframe = int(sys.argv[3])
    numfr = int(sys.argv[2])

#userproofing

if(numfr + startframe > totlines):
    print "Invalid frame parameters"
    sys.exit()

#print description
print "Postprocessing %d frames of " %numfr + fname 
#setup progress bar
progwidth = 60
sys.stdout.write("[%s]"%(" " * progwidth))
sys.stdout.flush()
sys.stdout.write("\b" * (progwidth +1))

#start program 
tstart = time.time()
f=open(fname, 'r')
#skips to the start frame
for i in range(0,6+startframe):
    f.readline()

#this is the main loop
#NOTE each loop uses a different letter for indexing
for k in range(0,numfr):
    l = f.readline()
    s = l.split()

    fr = float(s[0])#frame
    time1 = float(s[1])
    npts = (len(s) - 2) / 3
    ndata =[[],[],[]]
    comtp = comt
    comt = [0,0,0]
    #this loop finds the center of mass and catalouges each point
    for j in range(0, npts):
        x = float(s[2 + 3*j])
        y = float(s[3 + 3*j])
        z = float(s[4 + 3*j])
        comt[0] += x/npts
        comt[1] += y/npts
        comt[2] += z/npts
        #x,y,z for the clustering algorithm
        ndata[0].append(x)
        ndata[1].append(y)
        ndata[2].append(z)
        
    
    #calculating the Interia Tensor for the whole array
    rotcomp = rotcom
    JCOM = np.array([[0,0,0],[0,0,0],[0,0,0]]) 
    for w in range(0, npts):
        x = float(s[2 + 3*j])
        y = float(s[3 + 3*j])
        z = float(s[4 + 3*j])
        r = [x-comt[0],y-comt[1],z-comt[2]]
        
        JCOM = JCOM + (r[0]*r[0] + r[1]*r[1] +r[2]*r[2])*np.identity(3) - np.array([[r[0]*r[0], r[0]*r[1], r[0]*r[2]],[r[1]*r[0], r[1]*r[1], r[1]*r[2]],[r[2]*r[0], r[2]*r[1], r[2]*r[2]]])
    
    #calculating angular velocity
    #change in orientation of interia tensor
    qrcom = np.linalg.qr(JCOM)
    rotcom = qrcom[0]

    deltat = np.dot(rotcom, np.transpose(rotcomp))
    thetacom = matrixToAngle(deltat)
    dthetat = np.append(dthetat, thetacom)
    
    #calculating the 'spread-outness'
    #by taking the 'R' part of the
    #QR decomposition
    spread = qrcom[1]
    sprx = np.dot(spread,[1,0,0])
    spry = np.dot(spread,[0,1,0])
    sprz = np.dot(spread,[0,0,1])
    
    sprxx = sprx[0] + spry[0] + sprz[0]
    spryy = sprx[1] + spry[1] + sprz[1]
    sprzz = sprx[2] + spry[2] + sprz[2]
    
    sprtot = math.sqrt(sprxx*sprxx + spryy*spryy + sprzz*sprzz)
    timespread = np.append(timespread, sprtot)
    
    
    #calculate the axis of rotation (unit vector)
    axis = matrixToAxis(deltat, thetacom)

    mag = math.sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2])

    axisx = np.append(axisx, axis[0])
    axisy = np.append(axisy, axis[1])
    axisz = np.append(axisz, axis[2])
    
        
    #calculating the COM velocity
    comx = np.append(comx,comt[0])
    comy = np.append(comy,comt[1])
    comz = np.append(comz,comt[2])
    
    dx = comt[0] - comtp[0]
    dy = comt[1] - comtp[1]
    dz = comt[2] - comtp[2]
    
    xdot = np.append(xdot, dx)
    ydot = np.append(ydot, dy)
    zdot = np.append(zdot, dz)
    speedt = np.append(speedt, math.sqrt(dx*dx + dy*dy + dz*dz))
    #cluster needs an numpy array instead of a regular one
    data = np.asarray(ndata) 
    #The two lines below are for clustering. 
    thresh = 660 #works best from guessing and checking was 550 . no dynamically changes based on the distance between center of mass 
    clusters = hcluster.fclusterdata(np.transpose(data), thresh, criterion="distance")
    nclusters = findLargest(clusters)
    numclusters = np.append(numclusters, nclusters)
    #calculate center of mass for each cluster (person)
    #lots of parellel arrays(for speed)
    #the 0th index corresponds to the 1st dancer
    #TODO generalize to n dancers 
    comp = com
    com =[[0,0,0],[0,0,0]]
    dcom = com
    numpts = [0,0]
    xt = [0,0]
    yt = [0,0]
    zt = [0,0]
        #this loop sepearates the the two clusters(dancers)
    #by means of parellel arrays
    for i in range(0, npts):
        if(clusters[i] == 1):
            tg = 0
        else:
            tg = 1

        numpts[tg] = numpts[tg] + 1
        xt[tg] = xt[tg] + ndata[0][i]
        yt[tg] = yt[tg] + ndata[1][i]
        zt[tg] = zt[tg] + ndata[2][i]

        
    #this loop calculates the center of mass and the 
    #change is the com for each cluster.
    for m in range(0,2):
        try:
            com[m] = [xt[m]/numpts[m], yt[m]/numpts[m], zt[m]/numpts[m]]
            #Change in the x,y,z components of each dancer-cluster.
            dcom[m] = [com[m][0] - comp[m][0], com[m][1] - comp[m][1], com[m][2] - comp[m][2]]
            
        
        except ZeroDivisionError:
            pass
        except ValueError:
            pass

        
    jarp = jar
    rotp = rot
#need to redeclare them up here so we aren't adding ad infinitum
    J1 = np.array([[0,0,0],[0,0,0],[0,0,0]])#Inerita Tensors
    J2 = np.array([[0,0,0],[0,0,0],[0,0,0]])
    jar = [J1,J2]
    
    R1 = np.array([[0,0,0],[0,0,0],[0,0,0]])
    R2 = np.array([[0,0,0],[0,0,0],[0,0,0]])
    rot = [R1,R2]


    #calculate the Inertia Tensor
    for q in range(0, npts):
        x = float(s[2 + 3*q])
        y = float(s[3 + 3*q])
        z = float(s[4 + 3*q])
        
        #clusters are either 1 or 2, not 0 or 1
        if (clusters[q] == 1):
            tg = 0
        else:
            tg = 1
        
        r = [x-com[tg][0], y-com[tg][1], z-com[tg][2]]

                #Armand's Inertia tensor algorithm
                
            
        jar[tg] = jar[tg] + (r[0]*r[0] + r[1]*r[1] +r[2]*r[2])*np.identity(3) - np.array([[r[0]*r[0], r[0]*r[1], r[0]*r[2]],[r[1]*r[0], r[1]*r[1], r[1]*r[2]],[r[2]*r[0], r[2]*r[1], r[2]*r[2]]])
        
        
    rot[0] = np.linalg.qr(jar[0])[0]
    rot[1] = np.linalg.qr(jar[1])[0]
    
    deltar1 = np.dot(rot[0], np.transpose(rotp[0]))
    deltar2 = np.dot(rot[1], np.transpose(rotp[1]))
    
    ang1 = matrixToAngle(deltar1)
    ang2 = matrixToAngle(deltar2)
    
    dR1 = append(dR1, ang1)
    dR2 = append(dR2, ang2)
    #we found the angle
    #
    #
    x1 = np.append(x1, com[0][0])
            #adding data to the arrays used for plotting
    x1dot = np.append(x1dot, dcom[0][0])
    y1dot = np.append(y1dot, dcom[0][1])
    z1dot = np.append(z1dot, dcom[0][2])

    x2dot = np.append(x2dot, dcom[1][0])
    y2dot = np.append(y2dot, dcom[1][1])
    z2dot = np.append(z2dot, dcom[1][2])
    
    
    #c1c2 = math.fabs(np.linalg.norm(com[0]) - np.linalg.norm(com[1]))
    #updating progressbar
    if k % (numfr/progwidth) == 0:
        sys.stdout.write("-")
        sys.stdout.flush()
        #print "Frame ",k
        #print mag
                
    #Display qunatities of interest...
    """
    print "Center of mass for dancer 1:"
    print com[0]
    print "Center of mass for dancer 2:"
    print com[1]
    print "Center of mass for system:"
    print comt
    print "X Velocity for dancer 1: %f", com[0][0]
    """

######Plotting###########
#plot the xy
#subplot(kmn) k-> total numner m->1, n-> level


ymax = 150
"""
plt.figure("Raw Velocity Components for whole COM")

plt.subplot(311)
plt.plot(xdot)

plt.subplot(312)
plt.plot(ydot)

plt.subplot(313)
plt.plot(zdot)
"""
xdotf = movingaverage(xdot,400)
ydotf = movingaverage(ydot,400)
zdotf = movingaverage(zdot,400)

plt.figure("Filtered Velocity Components for whole COM")

plt.subplot(311)
plt.plot(xdotf)

plt.subplot(312)
plt.plot(ydotf)

plt.subplot(313)
plt.plot(zdotf)

#plt.subplot(414)
#plt.plot(speedt)
"""
plt.figure("Raw Position Components for COM")

plt.subplot(311)
plt.plot(comx)

plt.subplot(312)
plt.plot(comy)

plt.subplot(313)
plt.plot(comz)
"""
comxf = movingaverage(comx,30)
comyf = movingaverage(comy,30)
comzf = movingaverage(comz,30)

plt.figure("Filtered Position Coordinates for COM")
plt.subplot(311)
plt.plot(comxf)

plt.subplot(312)
plt.plot(comyf)

plt.subplot(313)
plt.plot(comzf)
"""
x1dotf = movingaverage(x1dot,150)
y1dotf = movingaverage(y1dot,150)
z1dotf = movingaverage(z1dot,150)
x1f = movingaverage(x1,100)

plt.figure("Position Components for Dancer 1")
plt.plot(x1f)

plt.figure("Velocity Components for Dancer 1")
plt.plot(x1dotf)

plt.subplot(311)
plt.plot(x1dotf)

plt.subplot(312)
plt.plot(y1dotf)

plt.subplot(313)
plt.plot(z1dotf)


plt.figure("Angular Velcotiy Components for Dancer 1")
plt.subplot(311)
plt.plot(J1xdot)

plt.subplot(312)
plt.plot(J1ydot)

plt.subplot(313)
plt.plot(J1zdot)

plt.figure("Angular Velcotiy Components for Dancer 2")
plt.subplot(311)
plt.plot(J2xdot)

plt.subplot(312)
plt.plot(J2ydot)

plt.subplot(313)
plt.plot(J2zdot)
"""
"""
drf = movingaverage(dR1,250)
plt.figure("Total Angular Velocity for Dancer1")
plt.plot(drf)
"""
#QR Angular Velocity Plot
dtcf = movingaverage(dthetat, 75)
plt.figure("QR Angular Velocity for ALL Points")
plt.plot(dtcf)

#Spread Plot
tsprf = movingaverage(timespread, 75)
plt.figure("QR Scaling For ALL Points")
plt.plot(tsprf)
#Axis Plot
axxf = movingaverage(axisx, 75)
ayyf = movingaverage(axisy, 75)
azzf = movingaverage(axisz, 75)
plt.figure("Axis Of Rotation For ALL Points")
plt.subplot(311)
plt.plot(axxf)
plt.subplot(312)
plt.plot(ayyf)
plt.subplot(313)
plt.plot(azzf)
#plot number of clusters (should be 2)
plt.figure("Number of Clusters")
plt.plot(numclusters)

#fps stats
tst = time.time()

ttotal = tst - tstart
fps = numfr / ttotal

print "Processed %d frames in %f seconds" %(numfr, ttotal)
print "rate: %f frames per second" %fps

plt.show()
