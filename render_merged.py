import Main as fractal
import numpy as np
from matplotlib import pyplot as pt
from matplotlib.colors import LinearSegmentedColormap 

#define peramiters
x1,x2 = -1,1
y1,y2 = -1,1
resolution = 100
method = 0
frames = 10
save_name = "frames/merge_test"

#define functions
def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])

def g(x,y):
    return np.array([(x**3)-(3*x)*(y**2)-(2*x)-2,3*(x**2)*y-(y**3)-2*y])
    
def get_function(f,g,l):
    def ret(x,y):
        return np.array(f(x,y))+l*np.array(g(x,y))
    return ret

#generate frames:
for i in range(frames):
    func = get_function(f,g,i/frames)
    frac = fractal.fractal2D(func)
    A,B  = frac.plot(resolution, x1, x2, y1, y2, method, False)
    Y, X=np.meshgrid(np.linspace(x1,x2,resolution),np.linspace(y1,y2,resolution),indexing='ij')
    #plot and save plot
    pt.clf() #clear current figure
    B = ((B-1)%10)/10
    color_map = pt.cm.get_cmap('Set1')
    color_map.set_under('0')
    pt.pcolormesh(X,Y,A,cmap=pt.cm.Set1,edgecolor = 'face', linewidth=0,rasterized=True, vmin=0,vmax = 9)
    #make custom color map
    custom = [(0,(0, 0, 0, 0)), (1,(0, 0, 0, 0.5))]
    my_cmap = LinearSegmentedColormap.from_list('something',custom)
    pt.pcolormesh(X,Y,B,edgecolor = 'face',cmap=my_cmap,linewidth=0,rasterized=True)
    pt.savefig(save_name+str(i))
    print(i+1,"  finished")
