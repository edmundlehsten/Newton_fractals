# Newton_fractals
NUMA01 Project 4
import time 
"""need to import this"""

"""
for function k
"""
nk=[]
sk=[]
ok=[]
N=[]
"""
epty lists to collect data
"""
for i in range (50,550, 50):
    N.append(i)
    start = time.time()
    k.plot(i,-1,1,-1,1,0)
    end = time.time()
    nk.append(end-start)
    
    start = time.time()
    k.plot(i,-1,1,-1,1,1)
    end = time.time()
    sk.append(end-start)
    
    start = time.time()
    k.plot(i,-1,1,-1,1,2)
    end = time.time()
    ok.append(end-start)
    
    
pt.plot(N,nk,'r', label='Newton method')
pt.plot(N,sk,'b',label='Simplified Newton method')
pt.plot(N,ok,'g',label='Optimized Newton method')
pt.xlabel('Resolution')
pt.ylabel('time (s)')
pt.title('Execution time for different implementations')
pt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)

"""
for function l
"""

nl=[]
sl=[]
ol=[]
N=[]

for i in range (50,550, 50):
    N.append(i)
    start = time.time()
    l.plot(i,-1,1,-1,1,0)
    end = time.time()
    nl.append(end-start)
    
    start = time.time()
    l.plot(i,-1,1,-1,1,1)
    end = time.time()
    sl.append(end-start)
    
    start = time.time()
    l.plot(i,-1,1,-1,1,2)
    end = time.time()
    ol.append(end-start)
    
    
pt.plot(N,nl,'r', label='Newton method')
pt.plot(N,sl,'b',label='Simplified Newton method')
pt.plot(N,ol,'g',label='Optimized Newton method')
pt.xlabel('Resolution')
pt.ylabel('time (s)')
pt.title('Execution time for different implementations')
pt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)

"""
for function m
"""
nm=[]
sm=[]
om=[]
N=[]

for i in range (50,550, 50):
    N.append(i)
    start = time.time()
    m.plot(i,-1,1,-1,1,0)
    end = time.time()
    nm.append(end-start)
    
    start = time.time()
    m.plot(i,-1,1,-1,1,1)
    end = time.time()
    sk.append(end-start)
    
    start = time.time()
    m.plot(i,-1,1,-1,1,2)
    end = time.time()
    om.append(end-start)
    
    
pt.plot(N,nm,'r', label='Newton method')
pt.plot(N,sm,'b',label='Simplified Newton method')
pt.plot(N,om,'g',label='Optimized Newton method')
pt.xlabel('Resolution')
pt.ylabel('time (s)')
pt.title('Execution time for different implementations')
pt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
