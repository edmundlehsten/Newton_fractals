import newtfraclib as newtonFractal
import numpy as np
from matplotlib import pyplot as pt
from matplotlib.colors import LinearSegmentedColormap 

start_xdim = 0.2
start_ydim = 0.2
focus = (0,0.4)
frames = 150
end_xdim = 0.001
end_ydim = 0.001
resolution = 200
directory = "frames/fin"
xm = (end_xdim - start_xdim) / frames 
ym = (end_ydim - start_ydim) / frames 


def g(x,y):
    return np.array([(x**3)-(3*x)*(y**2)-(2*x)-2,3*(x**2)*y-(y**3)-2*y])
fractal=newtonFractal.fractal2D(g)
for i  in range(frames):
    #get A and B plot
    x1 = focus[0]-(xm*i+start_xdim)/2
    x2 = focus[0]+(xm*i+start_xdim)/2
    y1 = focus[1]-(ym*i+start_ydim)/2
    y2 = focus[1]+(ym*i+start_ydim)/2

    Y, X=np.meshgrid(np.linspace(x1,x2,resolution),np.linspace(y1,y2,resolution),indexing='ij')
    A,B = fractal.plot(resolution,x1,x2,y1,y2,0,False)
    #plot and save plot
    pt.clf() #clear current figure
    B = ((B-1)%10)/10
    color_map = pt.cm.get_cmap('Set1')
    color_map.set_under('0')
    pt.pcolormesh(X,Y,A,cmap=pt.cm.Set1,edgecolor = 'face', linewidth=0,rasterized=True, vmin=0,vmax = A.max())
    #make custom color map
    custom = [(0,(0, 0, 0, 0)), (1,(0, 0, 0, 0.5))]
    my_cmap = LinearSegmentedColormap.from_list('something',custom)
    pt.pcolormesh(X,Y,B,edgecolor = 'face',cmap=my_cmap,linewidth=0,rasterized=True)
    pt.savefig(directory+str(i))
    print(i+1,"  finished")
print("Done")
