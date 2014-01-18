#Realtime clustering and plotting
#Gets mocapdata from a TCP/IP socket
#And plots it in real time
# with color coded clusters!
#J Lakowski 1/13/2014

import mocaplib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import socket

fig = plt.figure() #hi-res dpi=1000
plt.ion()
plt.show()

TCP_IP = "zx81.isl.uiuc.edu" #car8.isl.uiuc.edu:4712 for realtime
TCP_PORT = 4710

s = socket.socket() 
s.connect((TCP_IP, TCP_PORT))

while True:
    data = s.recv(8192)
    fs = mocaplib.parseFrame(data)
    clust = mocaplib.cluster(fs, 660)
    nclus = len(clust[2]) #get the number of clusters
    print nclus
    coms = [[]]*nclus #initialize array of the correct size
    print clust[2][0]
    for i in range(0,nclus):
        coms[i] = mocaplib.clusCOM(clust[2][i])
    
    ax = Axes3D(fig)
    ax.scatter(*clust[0], c=clust[1])
    #display an orange dot at the centers of mass
    for j in range(0,nclus):
        ax.scatter(coms[j][0], coms[j][1], coms[j][2], c='orange')

    plt.draw()
