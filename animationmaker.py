#animationmaker.py
#Jim Lakowski 7/15/2013-8/13/2013
#turns a list of points into an animation!
#and twitter?
#8/13 works well. renders 20000 plus frames,
#a progress bar would be really sick

import scipy.cluster.hierarchy as hcluster
import numpy
import fileinput
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import os
import twitter
import time
import gc
import scitools
#twitter stuff
consumer_key = "8URMpPDxJ2YAfGRVWVUGtQ"
consumer_secret = "7Y8d2ehwTGKYtDnocRZUSRP6k6GBFEW1ThOa2167vq0"
access_key = "398477484-gYw8bpydJKe6FtsYrPnCLasju89AR0gj7nyWDcmb"
access_secret = "lOGX2mM3WGtxXC86SubcBTc67LAfN2Lw05UGpXEsw"

api = twitter.Api()
api = twitter.Api(consumer_key=consumer_key, consumer_secret=consumer_secret, access_token_key=access_key, access_token_secret=access_secret)

#ion()
blocknum = 0
fig_hndl = figure()
frames= []
f=open('simonkirstie-130328-14531-Unnamed.trc', 'r')
startframe = 6 #start of frame data
numframes = 22000 #originally 2210

tstart = time.time()

for i in range(0,startframe):
    f.readline()
    

filename = "animation%d"%blocknum 
for k in range(0,numframes):
    l = f.readline()
    s = l.split()

    fr = float(s[0])
    time = float(s[1])
    npts = (len(s) - 2) / 3
    ndata =[[],[],[]]
    comt = [0,0,0]
    for i in range(0, npts):
        x = float(s[2 + 3*i])
        y = float(s[3 + 3*i])
        z = float(s[4 + 3*i])
        comt[0] += x/npts
        comt[1] += y/npts
        comt[2] += z/npts

        ndata[0].append(x)
        ndata[1].append(y)
        ndata[2].append(z)
    data = numpy.asarray(ndata)

    thresh = 550
    clusters = hcluster.fclusterdata(numpy.transpose(data), thresh, criterion="distance")
#calculate center of mass for each cluster (person)
    com =[[0,0,0],[0,0,0]]
    numpts = [0,0]
    xt = [0,0]
    yt = [0,0]
    zt = [0,0]

    for i in range(0, npts):
        if(clusters[i] == 1):
            tg = 0
        else:
            tg = 1

        numpts[tg] = numpts[tg] + 1
        xt[tg] = xt[tg] + ndata[0][i]
        yt[tg] = yt[tg] + ndata[1][i]
        zt[tg] = zt[tg] + ndata[2][i]

    for i in range(0,2):
        try:
            com[i]=[xt[i]/numpts[i], yt[i]/numpts[i], zt[i]/numpts[i]]
        except ZeroDivisionError:
            pass
    #print "Center of mass for dancer 1:"
    #print com[0]
    #print "Center of mass for dancer 2:"
    #print com[1]
    #print "Center of mass for system:"
    #print comt

    # plotting
    fig = plt.figure(dpi=1000) #hi-res

    #ax = fig.add_subplot(111,projection='3d')
    ax = Axes3D(fig)
    ax.scatter(*data, c=clusters)
    ax.scatter(com[0][0], com[0][1], com[0][2], c='orange') #com1
    ax.scatter(com[1][0], com[1][1], com[1][2], c='orange') #com2
    ax.scatter(comt[0], comt[1], comt[2], c = 'purple')#total COM
    ax.set_autoscale_on(False)
    ax.set_xlim(left = -2000, right = 2000)
    ax.set_ylim(bottom = -2000, top = 2000)
    ax.set_zlim(bottom =-500, top = 1400)
    title = "threshold: %f, number of clusters: %d, frame %d" % (thresh, len(set(clusters)), k + startframe)
    plt.title(title)
    #gives the files the right name to be rendered in order...
    if k <10:
        fname = 'renders/_frame0000%d.png'%k
    elif k < 100 and k > 9:
        fname = 'renders/_frame000%d.png'%k
    elif k < 1000 and k > 99:
        fname = 'renders/_frame00%d.png'%k
    elif k < 10000 and k> 999:
        fname = 'renders/_frame0%d.png'%k
    elif k < 100000 and k >9999:
        fname = 'renders/_frame%d.png'%k
    plt.savefig(fname , bbox_inches=0, dpi=300)
    if k % 100 ==0:
        print "rendered frame %d"%k
    frames.append(fname)
    plt.close()#THIS LINE IS KEY!! 
    #plt.show()
   
#timestop = time.time()
#animation
os.system("mencoder 'mf://renders/_frame*.png' -mf type=png:fps=100 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation13.avi")

for fname in frames: 
    os.remove(fname)


#timetotal = timestop - tstart
#fps = numframes / timetotal
#print "Rendered at a rate of %f frames per second" %fps        

