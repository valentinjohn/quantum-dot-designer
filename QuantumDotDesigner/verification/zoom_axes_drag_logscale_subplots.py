import matplotlib.pyplot as plt
import numpy as np

class ZoomPan:
    """
    combine two examples from  https://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel
    in order to (1) zoom x-axis, y-axis, and both axes with mouse scroll (2) drag with mouse click  
    (3) zoom and drag a plot that is log scale
    """
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None
        self.zoomaxs = []
        self.panaxs = []
        self.panax = None

    def zoom_factory(self, axs=[], base_scale = 1.4, axiszoomrange = 0.4):
        if isinstance(axs,np.ndarray): axs = list(axs.flatten())
        elif not isinstance(axs,list): axs = [axs]       
        if len(axs)==0: return None        
        self.zoomaxs = axs
        def zoom(event):
            x = event.x 
            y = event.y        
            
            # cur_xlim = ax.get_xlim()
            # cur_ylim = ax.get_ylim()
            # xdata = event.xdata # get event x location
            # ydata = event.ydata # get event y location
            
            if event.inaxes is not None:
                if event.inaxes in self.zoomaxs:
                    ax = event.inaxes
                else:
                    return
            else:
                pos_rel = []
                for i, ax_iter in enumerate(self.zoomaxs):
                    tranP2A = ax_iter.transAxes.inverted().transform
                    xa,ya = tranP2A((x,y))
                    # xylist += [ [xa,ya] ]
                    if 0 <= xa <= 1 and -axiszoomrange <= ya <= 0:
                        pos_rel += [   ['xaxis', ya, ax_iter]  ]
                    if 0 <= ya <= 1 and -axiszoomrange <= xa <= 0:
                        pos_rel += [   ['yaxis', xa, ax_iter]  ]
                if len(pos_rel) == 0:
                    return
                else:
                    i_pos_rel = np.argmax(np.array([ item[1] for item in pos_rel]))
                    ax = pos_rel[i_pos_rel][2]

            #convert pixels to axes
            tranP2A = ax.transAxes.inverted().transform
            #convert axes to data limits
            tranA2D= ax.transLimits.inverted().transform
            #convert the scale (for log plots)
            tranSclA2D = ax.transScale.inverted().transform
            
            if event.button == 'down':
                # deal with zoom in
                scale_factor = base_scale
            elif event.button == 'up':
                # deal with zoom out
                scale_factor = 1/base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1
            #get my axes position to know where I am with respect to them
            xa,ya = tranP2A((x,y))
            zoomx = False
            zoomy = False 
            if(ya < 0):
                if(xa >= 0 and xa <= 1):
                    zoomx = True
                    zoomy = False
            elif(ya <= 1):
                if(xa <0): 
                    zoomx = False
                    zoomy = True
                elif(xa <= 1):
                    zoomx = True
                    zoomy = True
                else:
                    zoomx = False
                    zoomy = True
            else:
                if(xa >=0 and xa <= 1):
                    zoomx = True
                    zoomy = False
            # print(ax, event.inaxes)
            new_alimx = (0,1)
            new_alimy = (0,1)
            
 
            if(zoomx):
                new_width =  scale_factor    
                relx = 1-xa
                new_alimx = [xa - new_width * (1-relx), xa + new_width * (relx)]
            if(zoomy):
                new_width =  scale_factor    
                rely = 1-ya
                new_alimy = [ya - new_width * (1-rely), ya + new_width * (rely)]
            #now convert axes to data
            new_xlim0,new_ylim0 = tranSclA2D(tranA2D((new_alimx[0],new_alimy[0])))
            new_xlim1,new_ylim1 = tranSclA2D(tranA2D((new_alimx[1],new_alimy[1])))

            ax.figure.canvas.toolbar.push_current()
            ax.set_xlim([new_xlim0,new_xlim1])
            ax.set_ylim([new_ylim0,new_ylim1])
            ax.figure.canvas.draw()            
        # fig = ax.get_figure()   
        fig = axs[0].get_figure() # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom

    def pan_factory(self, axs=[]):
        if isinstance(axs,np.ndarray): axs = list(axs.flatten())
        elif not isinstance(axs,list): axs = [axs]
        self.panaxs = axs
        def onPress(event):
            # if event.inaxes != ax: return
            i = 0
            while self.panaxs[i] != event.inaxes:
                i += 1
                if i == len(self.panaxs): return

            self.panax = event.inaxes            
            self.press = self.x0, self.y0, event.x, event.y
            self.xpress,self.ypress = event.x, event.y
        def onRelease(event):
            if self.press is None: return
            self.press = None
            if event.inaxes != self.panax: return
            event.inaxes.figure.canvas.draw()
            
        def onMotion(event):
            if self.press is None: return
            if event.inaxes != self.panax: return
            ax = self.panax
            #convert pixels to axes
            tranP2A = ax.transAxes.inverted().transform
            #convert axes to data limits
            tranA2D= ax.transLimits.inverted().transform
            #convert the scale (for log plots)
            tranSclA2D = ax.transScale.inverted().transform
            
            xa,ya = tranP2A((event.x, event.y))
            selfxa,selfya = tranP2A((self.xpress,self.ypress))
            dxa = xa - selfxa
            dya = ya - selfya
            
            new_xlim0,new_ylim0 = tranSclA2D(tranA2D((-dxa,-dya)))
            new_xlim1,new_ylim1 = tranSclA2D(tranA2D((-dxa+1,-dya+1)))            
            

            self.xpress, self.ypress = event.x, event.y
            ax.figure.canvas.toolbar.push_current()
            ax.set_xlim(new_xlim0, new_xlim1)
            ax.set_ylim(new_ylim0, new_ylim1)
            ax.figure.canvas.draw()
        fig = axs[0].get_figure()

        # attach the call back
        fig.canvas.mpl_connect('button_press_event',onPress)
        fig.canvas.mpl_connect('button_release_event',onRelease)
        fig.canvas.mpl_connect('motion_notify_event',onMotion)

        #return the function
        return onMotion
    
'''
#%%

fig, axs = plt.subplots(2,3,figsize=(10,5))
ax = axs[0,0]
x = np.linspace(0,10,101)
y = np.linspace(0,10,101)
ax.pcolormesh(x,y, np.sin(x[:,np.newaxis]*y))

ax = axs[0,1]
x = np.linspace(0,10)
ax.plot(x, np.exp(x))
ax.set_yscale('log')

ax = axs[1,0]
x = np.linspace(0,10,101)
y = np.linspace(0,10,101)
ax.contour(x,y, np.sin(x[:,np.newaxis]*y))

ax = axs[1,1]
x = np.linspace(0,10)
ax.plot(x, np.exp(x))
# ax.set_yscale('log')

ax = axs[1,2]
x,y,s,c = np.random.rand(4,50)
s *= 200
ax.scatter(x,y,s,c)


zp = ZoomPan()
# figZoom = zp.zoom_factory(list(axs.flatten()), base_scale = 1.2)
# figPan = zp.pan_factory(   list(axs.flatten())   )
figZoom = zp.zoom_factory( axs, base_scale = 1.2)
figPan = zp.pan_factory(   axs   )
plt.show()



#%%
axs = []

plt.figure()
plt.subplot(121)
x = np.linspace(0,10,101)
y = np.linspace(0,10,101)
plt.pcolormesh(x,y, np.sin(x[:,np.newaxis]*y))
# axs += [plt.gca()]

plt.subplot(122)
x = np.linspace(0,10)
plt.plot(x, np.exp(x))
plt.yscale('log')
axs += [plt.gca()]

zp = ZoomPan()
figZoom = zp.zoom_factory(axs, base_scale = 1.2)
figPan = zp.pan_factory(  axs  )
plt.show()



#%%
fig, ax = plt.subplots()
x = np.linspace(0,10)
ax.plot(x, np.exp(x))
# ax.set_yscale('log')
zp = ZoomPan()
figZoom = zp.zoom_factory(ax)
figPan = zp.pan_factory(ax)
plt.show()
    
#%%

fig = plt.figure()

ax = fig.add_subplot(111, xlim=(0,1), ylim=(0,1), autoscale_on=False)

ax.set_title('Click to zoom')
x,y,s,c = np.random.rand(4,200)
s *= 200

ax.scatter(x,y,s,c)
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = 1.5)
figPan = zp.pan_factory(ax)
plt.show()



#%%
fig, ax = plt.subplots()
x = np.linspace(0,10,101)
y = np.linspace(0,10,101)
ax.contour(x,y, np.sin(x[:,np.newaxis]*y))
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = 1.4)
figPan = zp.pan_factory(ax)
plt.show()

'''
