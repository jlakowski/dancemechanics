#mocablib
#a python module with helpful functions 
#for interfacing with and analyzing data
#from the ISL's mocap array
#
# James Lakowski December 23, 2013
#requirements
#python 2.6
#numpy
#scipy

#TODO: Write documentation
import math
import numpy as np
import scipy.cluster.hierarchy as hcluster
import re
import copy

# helper method used by the other methods in this
#library 
def matrixToAngle(mat):
    #using a formula
    #http://en.wikipedia.org/wiki/Axis-angle_representation
    a = math.acos((np.trace(mat)-1)/2)
    return a
#a helper method used by other methods
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

#From Armand's dancer.py (Spring 2013)
def is_number(s):
    # A utility function to make certain that an argument is a float
    try:
        float(s)
        return True
    except ValueError:
        return False

#A method which finds the Center of mass for a points vector
#points is the raw list of x y z coordinates from the mocap system
def findCOM(points):
    comt = [0,0,0]
    npts = (len(points) -2 )/3
    for j in range(0, npts):
        x = float(points[2 + 3*j])
        y = float(points[3 + 3*j])
        z = float(points[4 + 3*j])
        comt[0] += x/npts
        comt[1] += y/npts
        comt[2] += z/npts
    
    return comt
#finds the center of mass for a cluster of points
# a cluster is an array of [x,y,z] points
# 
def clusCOM(cluster):
    com = [0,0,0]
    for i in range(len(cluster)):
        com[0] += cluster[i][0]/len(cluster)
        com[1] += cluster[i][1]/len(cluster)
        com[2] += cluster[i][2]/len(cluster)
    return com
#A method which  calculates the inertia tensor for a points vector
#com is a 3 dimensional array that holds the x,y,z center of mass
#Returns a NUMPY array with the Interia Tensor of the points
def findInertiaTensor(points, com):
    npts = (len(points)-2)/3
    JCOM = np.array([[0,0,0],[0,0,0],[0,0,0]]) 
    for w in range(0, npts):
        x = float(points[2 + 3*w])
        y = float(points[3 + 3*w])
        z = float(points[4 + 3*w])
        r = [x-com[0],y-com[1],z-com[2]]
        
        JCOM = JCOM + (r[0]*r[0] + r[1]*r[1] +r[2]*r[2])*np.identity(3) - np.array([[r[0]*r[0], r[0]*r[1], r[0]*r[2]],[r[1]*r[0], r[1]*r[1], r[1]*r[2]],[r[2]*r[0], r[2]*r[1], r[2]*r[2]]])
    
    return JCOM

#uses a QR Factorization to calculate the change in orientation
#of the Inertia Tensor which we use as angular momentum
#previousIT is the Inertia Tensor from the last frame
#currentIT is the Inertia Tensor from the current frame
#returns a vector in the following form
# [change in orientation of the two ITs,
#  unit vector normal to that rotation,
#  the 'spreadoutness' of the CURRENT point cloud]
def findAngularVelocity(previousIT, currentIT):
    
    qrp = np.linalg.qr(previousIT)
    rotp = qrp[0]
    qrc = np.linalg.qr(currentIT)
    rotc = qrcp[0]
    deltat = np.dot(rotcom, np.transpose(rotcomp))
    dtheta = matrixToAngle(deltat)
    axis = matrixToAxis(deltat, dtheta)
    #Using the R portion of the qr factorization
    #the Amount of space the point cloud takes up
    #can be estimated
    spread = qrc[1]
    
    sprx = np.dot(spread,[1,0,0])
    spry = np.dot(spread,[0,1,0])
    sprz = np.dot(spread,[0,0,1])
    
    sprxx = sprx[0] + spry[0] + sprz[0]
    spryy = sprx[1] + spry[1] + sprz[1]
    sprzz = sprx[2] + spry[2] + sprz[2]
    
    sprtot = math.sqrt(sprxx*sprxx + spryy*spryy + sprzz*sprzz)
    
    return[dtheta,axis,sprtot]

# A method which takes in a point cloud and returns back
# an array of point cluster arrays 
# 
# generalizes the clustering to n Dancers
#output in the following form
# [ points vector , list of numbers which correspond to the points vector, 
# array of clustered points]
def cluster(points, thresh):
    #the x,y,z points must first be separated out
    ndata = [[],[],[]]
    npts = (len(points)-2)/3
    for j in range(0,npts):
        x = float(points[2 + 3*j])
        y = float(points[3 + 3*j])
        z = float(points[4 + 3*j])
        ndata[0].append(x)
        ndata[1].append(y)
        ndata[2].append(z)
    data = np.asarray(ndata)

    clusterlist = hcluster.fclusterdata(np.transpose(data), thresh, criterion="distance")
    
    nclusters = findLargest(clusterlist)
    
    #initializes an array to the right size
    #http://stackoverflow.com/questions/7745562/appending-to-2d-lists-in-python
    clusters = [[] for i in range(nclusters)] 
    #assingns points to the correct cluster
    for i in range(0, npts):
        #print clusters[clusterlist[i]-1]
        
        clusters[clusterlist[i]-1].append([ndata[0][i],ndata[1][i],ndata[2][i]])
    return [data, clusterlist, clusters]


#takes two a 3-D Vectors
#And calculates the difference
def findVelocity(prev,curr):
    return [curr[0]-prev[0],curr[1]-prev[1],prev[2]-prev[2]]

# method that used TCP/IP and the
#'start-record' protocol we've been using to get a
# frame from the mocap server
# pretty much copied and pasted
# from dancer.py of Spring 2013
def parseFrame(rawIn):
    
     # If we can't keep up, there may be multiple records.  Just take the first one.
    records = rawIn.split("\n")
        
    for record in records:
        if re.search('start-of-record',record): # =='start-of-record':
            ss = record.split(",")
            break

        # Remove the first, "start-of-record" string
    ss.pop(0)
    fs = []
    
    for pt in ss:
        f2b = pt.replace('[','').replace(']','')
        if is_number( f2b):
            fs.append( float( f2b ) )    # pt.replace('[','').replace(']','') ) )
        else:
            print 'bad float', f2b
            print record
        
    return fs
