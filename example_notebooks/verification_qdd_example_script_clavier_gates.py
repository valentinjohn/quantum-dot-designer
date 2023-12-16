# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 11:43:37 2023

@author: cawang
"""
import os
import sys
# sys.path.append('../')
path_notebook=str(os.getcwd()).replace('\\example_notebooks','')
sys.path.insert(1, path_notebook)



#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#%%  finding disconnected fanout lines
from QuantumDotDesigner.verification.verification import cellconnection

layer_name_idx = {}
for item in qda.elements.values():
    # print(item['layer'].name)
    layer_name_idx[item['layer'].name] = set()
    layer_name_idx[item['layer'].name].add( item['layer'].fine )
    layer_name_idx[item['layer'].name].add( item['layer'].coarse )
print('layers overview')
print(layer_name_idx , '\n')
layer_mapping = {
    (6,5): 5,
    (22,21): 21,
    }

# gives the gdstk.Cell object and the layer mapping layer_mapping
# union the polygons in the same layer
# after union, calculate the number of polygons and their area
# return the information as a dictionary 
# By looking at the numbers and area of polygons, this could be helpful in 
# finding disconnected fanout lines.

# printflag=True for printing the result in the console
results = cellconnection(qda.main_cell, layer_mapping, printflag=True)

#%% the shortest distance from a group of points to a line segment
from QuantumDotDesigner.verification.verification import linesection_point_distance

# define a line segment 
# left end  (x,y)=(-1,0)
# right end (x,y)=(1,0)
linepoint1 = [-1,0]
linepoint2 = [1,0]

# define three points
# x coordinate = 0.1, 0.1, 1.1
# y coordinate =   1,   2,   1
pointxs = [0.1, 0.1, 1.1]
pointys = [  1,   2,  1]
dis, xy, dis_new, xy_new, shortestpoint_idx = linesection_point_distance(linepoint1, linepoint2, pointxs, pointys )

# For these three points, the shortest distance to the infinite line going through linepoint1 - linepoint2
print(dis)

# For these three points, the closest points on the infinite line
print(xy)

# For these three points, the shortest distance to the line segment linepoint1 - linepoint2
print(dis_new )

# For these three points, the closest points on the line segment
print(xy_new )

# calculate the distance manually, to verify that dis_new is the shortest distance for these three points
dx = xy_new[:,0]-np.array(pointxs)
dy = xy_new[:,1]-np.array(pointys)
print( np.sqrt(dx**2+dy**2) )

# Find out if xy_new is on the line or on the end points
# 0 means on the line, 1 means linepoint1, 2 means linepoint2
print(shortestpoint_idx)





#%%  Visualize the shortest distance of a layer in a given gdstk.Cell
from QuantumDotDesigner.verification.verification import layermap, uniongdstk, polys_distance





layer_mapping = {
    (6,5): 5,
    (22,21): 22,
    (32,31): 31,
    }

# make a list of polygons which are in qda.main_cell 
#  the polygons of layer5 and layer6  becomes layer5
# layer 22 and 21 becomes 22
# layer 32 and 31 becomes 31
polyslist = layermap(qda.main_cell, layer_mapping)

# polyslist[1] is a list of polygons (gdstk.Polygon)
polysU, Npolys, NpolysU  = uniongdstk(polyslist[1]) 
# polysU is the polygons after we union all the polygons in polyslist[1]  
# Npolys is the number of polygons in polyslist[1]
# NpolysU is the number of polygons in polysU 

    
    
fig, ax = plt.subplots(1,1, figsize=(6, 6 ))

# plot all the polygons in polysU
for i in range(len(polysU)):
    ax.add_patch(matplotlib.patches.Polygon(polysU[i].points, facecolor='C0', edgecolor='k'))


# compute the shortest distance between polygons in polysU
vals_f, xy_f, report_new = polys_distance(polysU, Npolys, NpolysU)

# the shortest distance of polygon_0 to other polygons
print(vals_f[0])

# the shortest distance of polygon_1 to other polygons
print(vals_f[1])

# the coordinates of a line segment reporesents the shortest distance of polygon_0 to other polygons
print(xy_f[0] )

# visualize the line segments of the shortest distance
for i in range(NpolysU):
    ax.plot( xy_f[i][:,0], xy_f[i][:,1], 'r', linewidth=1)    


ax.set_xlim(-0.4,0.4)
ax.set_ylim(-0.4,0.4)


# Enalbe zoom and drag functions by giving the axes to the following functions 
from QuantumDotDesigner.verification.zoom_axes_drag_logscale_subplots import ZoomPan
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = 1.4)
figPan = zp.pan_factory(ax)
plt.show()







#%% check intersecting edges of a single polygon

from QuantumDotDesigner.verification.verification import intersect_check


# np2darr represents a rectangle with width and height of 1
np2darr = np.array([[0, 0], [0, 1], [1,1], [1, 0]])
# Check the polygon defined by np2darr, to see if it has two edges intersect
intersect_flag, report = intersect_check(np2darr)
print(intersect_flag)
print(report)


# Here np2darr represents a polygon contains two intersections
np2darr = np.array([[0, 0], [1, 1], [2,0], [2,1], [1,0], [0, 1]])
intersect_flag, report = intersect_check(np2darr)
print(intersect_flag)
# the report contains the coordinates of four points which make the intersections
# the 0th intersection 
print(report['intersect'][0])
# the 1st intersection 
print(report['intersect'][1])

# also contains the positions in the numpy array defining the polygon 
print(report['intersectidx'])

#%% find the normal vectors pointing ourward of a polygon 
from QuantumDotDesigner.verification.verification import outward_normvec

# np2darr represents a rectangle with width and height of 1
np2darr = np.array([[0, 0], [0, 1], [1,1], [1, 0]])
# Check the polygon defined by np2darr, to see if the vertices are ordered clockwise or anticlockwise
clockwiseness, vec_out = outward_normvec(np2darr)
print(clockwiseness)
# compute the normal vectors at every vertex on the polygon, 
# which are normalized, perpendicular to the edges, and point outward
print(vec_out)


# np2darr represents a polygon contains intersections
np2darr = np.array([[0, 0], [1,1], [0,1], [1, 0]])
clockwiseness, vec_out = outward_normvec(np2darr)
# It returns these values if the polygon contains intersections
print(clockwiseness)
print(vec_out)


#%% check overlapping edges of a single polygon
from QuantumDotDesigner.verification.verification import edgeoverlap_check

#  Check if there are overlapping edges in a given polygon

np2darr = np.array([[0, 0], [0, 4], [4,4], [4, 0], [0, 0], [0.4, 0.8], [1, 1], [2, 1], [2, 2], [1, 2], [1,1], [0.4, 0.8]])
overlap3pts_flag, overlap4pts_flag, report = edgeoverlap_check(np2darr, precision=1e-14)

# plot the polygon
fig, ax = plt.subplots(1,1, figsize=(6, 6 ))
ax.add_patch(matplotlib.patches.Polygon(np2darr, facecolor='C0', edgecolor='none'))
ax.plot( np2darr[:,0], np2darr[:,1], 'r', linewidth=1, marker='o')
ax.plot( (np2darr[-1,0],np2darr[0,0]) , (np2darr[-1,1],np2darr[0,1]), 'r', linewidth=1)
ax.set_xlim(-1,5)
ax.set_ylim(-1,5)

print(overlap3pts_flag)
print(overlap4pts_flag)
# the report contains the coordinates of two overlapping edges
# the coordinates of the 0th overlapping edges 
print(report['overlap4pts'][0])
# the position in the numpy array np2darr
print(report['overlap4ptsidx'][0])

# the 1st overlapping edges 
print(report['overlap4pts'][1])
# the position in the numpy array np2darr
print(report['overlap4ptsidx'][1])




#%% check overlapping edges of a single polygon
from QuantumDotDesigner.verification.verification import edgeoverlap_check

#  Check if there are overlapping edges in a given polygon

np2darr = np.array([[-1, 0], [1, 0], [1,2], [1, 1], [-1, 1], ])
overlap3pts_flag, overlap4pts_flag, report = edgeoverlap_check(np2darr, precision=1e-14)

# plot the polygon
fig, ax = plt.subplots(1,1, figsize=(6, 6 ))
ax.add_patch(matplotlib.patches.Polygon(np2darr, facecolor='C0', edgecolor='none'))
ax.plot( np2darr[:,0], np2darr[:,1], 'r', linewidth=1, marker='o')
ax.plot( (np2darr[-1,0],np2darr[0,0]) , (np2darr[-1,1],np2darr[0,1]), 'r', linewidth=1)
ax.set_xlim(-1.5,1.5)
ax.set_ylim(-0.5,2.5)

print(overlap3pts_flag)
# the report contains the coordinates of two overlapping edges
# the coordinates of the 0th overlapping edges 
print(report['overlap3pts'])
# the position in the numpy array np2darr
print(report['overlap3ptsidx'])



#%% expand or shrink a given polygon
from QuantumDotDesigner.verification.verification import polygonbuffer

# Buffer operation. Want to do something similar to shapely.buffer()
# It expand or shrink a given polygon with a given distance 

np2darr = np.array([[0, 0], [0, 3], [3,3], [3, 0], [0, 0], [0.4, 0.8], [1, 1], [2, 1], [2, 2], [1, 2], [1,1], [0.4, 0.8]])

fig, ax = plt.subplots(1,1, figsize=(6, 6 ))
ax.add_patch(matplotlib.patches.Polygon(np2darr, facecolor='C0', edgecolor='none'))
ax.plot( np2darr[:,0], np2darr[:,1], 'k', linewidth=1, marker='o')
ax.plot( (np2darr[-1,0],np2darr[0,0]) , (np2darr[-1,1],np2darr[0,1]), 'k', linewidth=1)

# expand the polygon by 0.2 um
# round the corner with an angle parameter, resulting 3 line segments for a 90 degree corner
buffer_distance = 0.2
corner_arc_dtheta = 31*np.pi/180
polynew = polygonbuffer(np2darr, buffer_distance, corner_arc_dtheta)
# polynew is a list of 2d numpy arrays representing polygons

# plot the result
for i in range(len(polynew)):
    xyplot = np.vstack((polynew[i], polynew[i][0:1,:]))
    patch = ax.plot( xyplot[:,0], xyplot[:,1], 'r', linewidth=1)

# shrink the polygon by 0.2 um
# round the corner with a angle parameter, resulting 2 line segments for a 90 degree corner
buffer_distance = -0.2
corner_arc_dtheta = 46*np.pi/180
polynew = polygonbuffer(np2darr, buffer_distance, corner_arc_dtheta)

for i in range(len(polynew)):
    xyplot = np.vstack((polynew[i], polynew[i][0:1,:]))
    patch = ax.plot( xyplot[:,0], xyplot[:,1], 'C1', linewidth=1)

ax.set_xlim( min(np2darr[:,0])-0.5, max(np2darr[:,0])+0.5 )
ax.set_ylim( min(np2darr[:,1])-0.5, max(np2darr[:,1])+0.5 )

#%% expand or shrink a given polygon: debug mode
from QuantumDotDesigner.verification.verification import polygonbuffer

np2darr = np.array([[0, 0], [0, 3], [3,3], [3, 0], [0, 0], [0.4, 0.8], [1, 1], [2, 1], [2, 2], [1, 2], [1,1], [0.4, 0.8]])

fig, axs = plt.subplots(1,3, figsize=(15, 5 ))
for ax in axs:
    ax.add_patch(matplotlib.patches.Polygon(np2darr, facecolor='C0', edgecolor='none'))
    ax.plot( np2darr[:,0], np2darr[:,1], 'k', linewidth=1, marker='o')
    ax.plot( (np2darr[-1,0],np2darr[0,0]) , (np2darr[-1,1],np2darr[0,1]), 'k', linewidth=1)


# Buffer operation with extra options:
# debugflag=True enables debug mode, returns two lists: rec_debug and arc_debug
# rec_debug is a list of 2d numpy arrays. Each of them should be rectangles forming the expanding region
# arc_debug is a list of 2d numpy arrays. Each of them should be fan-shaped filling the corner of the expanded polygons.
# if buffer_distance>0, polynew is the union (gdstk.boolean('or') ) of np2darr, rec_debug and arc_debug
# if buffer_distance<0, polynew is the subtraction (gdstk.boolean('not') ) of np2darr by rec_debug and arc_debug
# boolean_precision is the precision for the functions gdstk.boolean()

buffer_distance = 0.3
corner_arc_dtheta = 31*np.pi/180
polynew, rec_debug, arc_debug = polygonbuffer(np2darr, buffer_distance, corner_arc_dtheta, 
                                              debugflag=True, boolean_precision=1e-6)


for i in range(len(polynew)):
    xyplot = np.vstack((polynew[i], polynew[i][0:1,:]))
    axs[0].plot( xyplot[:,0], xyplot[:,1], 'C1', linewidth=1)
    
for i in range(len(rec_debug)):
    xyplot = np.vstack((rec_debug[i].points, rec_debug[i].points[0:1,:]))
    axs[1].plot( xyplot[:,0], xyplot[:,1], 'C1', linewidth=1)
    axs[1].set_title('rec_debug')

for i in range(len(arc_debug)):
    xyplot = np.vstack((arc_debug[i].points, arc_debug[i].points[0:1,:]))
    axs[2].plot( xyplot[:,0], xyplot[:,1], 'C1', linewidth=1)
    axs[2].set_title('arc_debug')
    
for ax in axs:    
    ax.set_xlim( min(np2darr[:,0])-0.5, max(np2darr[:,0])+0.5 )
    ax.set_ylim( min(np2darr[:,1])-0.5, max(np2darr[:,1])+0.5 )


# Enalbe zoom and drag functions by giving the axes 'axs' to the following functions 
from QuantumDotDesigner.verification.zoom_axes_drag_logscale_subplots import ZoomPan
zp = ZoomPan()
figZoom = zp.zoom_factory(axs, base_scale = 1.4)
figPan = zp.pan_factory(axs)
plt.show()


#%% expand or shrink polygons in a layer of gdstk.main_cell

from QuantumDotDesigner.verification.verification import layermap, uniongdstk, polygonbuffer

buffer_distance = 0.02
corner_arc_dtheta = 31*np.pi/180

layer_mapping = {
    (6,5): 5,
    (22,21): 21,
    (32,31): 31,
    }

polyslist = layermap(qda.main_cell, layer_mapping)

# Pick polygons belonging to one of the layers in layer_mapping
polys = polyslist[1]
polysU, Npolys, NpolysU  = uniongdstk(polys) 

for j in range(len(polysU)):
    
    # convert gdstk.Polygon to a numpy array representing a single polygon
    poly_origin = np.array(polysU[j].points)
    
    # Buffer a polygon
    polynew = polygonbuffer(poly_origin, buffer_distance, corner_arc_dtheta)




    # Plot original polygons
    fig, ax = plt.subplots(1,1, figsize=(6, 6 ))
    patch = ax.add_patch(matplotlib.patches.Polygon(poly_origin, facecolor='C0', edgecolor='none'))
    
    # Plot polygons after buffer
    for k in range(len(polynew)):
        xyplot = np.vstack((polynew[k], polynew[k][0:1,:]))
        patch = ax.plot( xyplot[:,0], xyplot[:,1], 'r', linewidth=1)


    # Enalbe zoom and drag functions by giving the axes 'ax' to the following functions 
    from QuantumDotDesigner.verification.zoom_axes_drag_logscale_subplots import ZoomPan
    zp = ZoomPan()
    figZoom = zp.zoom_factory(ax, base_scale = 1.4)
    figPan = zp.pan_factory(ax)
    plt.show()


#%% expand or shrink polygons in multiple layers of gdstk.main_cell

from QuantumDotDesigner.verification.verification import layermap, uniongdstk, polygonbuffer

buffer_distance = 0.02
corner_arc_dtheta = 31*np.pi/180

layer_mapping = {
    (6,5): 5,
    (22,21): 21,
    (32,31): 31,
    }

polyslist = layermap(qda.main_cell, layer_mapping)


# Pick polygons belonging to multiple layers in layer_mapping

fig, ax = plt.subplots(1,1, figsize=(6, 6 ))

for i in range(len(polyslist)):
    polys = polyslist[i]
    polysU, Npolys, NpolysU  = uniongdstk(polys) 
    for j in range(len(polysU)):
        
        # convert gdstk.Polygon to a numpy array representing a single polygon
        poly_origin = np.array(polysU[j].points)
        
        # Buffer a polygon
        polynew = polygonbuffer(poly_origin, buffer_distance, corner_arc_dtheta)
    
        # Plot original polygons
        patch = ax.add_patch(matplotlib.patches.Polygon(poly_origin, facecolor=f'C{i}', edgecolor='none'))
        
        # Plot polygons after buffer
        for k in range(len(polynew)):
            xyplot = np.vstack((polynew[k], polynew[k][0:1,:]))
            patch = ax.plot( xyplot[:,0], xyplot[:,1], f'C{i}', linewidth=1)


# Enalbe zoom and drag functions by giving the axes 'ax' to the following functions 
from QuantumDotDesigner.verification.zoom_axes_drag_logscale_subplots import ZoomPan
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = 1.4)
figPan = zp.pan_factory(ax)
plt.show()
    

#%% minimal example python interactive plot
from QuantumDotDesigner.verification.verification import ClickSimplePolygon

# allows mouse click: left or right buttons

# np2darr = qda.main_cell.get_polygons()[0].points
np2darr = np.array([[1,-0.1], [1, -1], [-1, -1], 
                    [-0.5, -0.5], [0, -0.5],  [0, 0], [-0.5, 0],  [-0.5, -0.5],
                    [-1, -1], [-1,1], [1, 1], [1,0.1], [2,0.1], [2,1], [3,1], [3,-1],[1.2,-1],[1.2,-0.8],[2,-0.8],[2,-0.1],])


fig, ax = plt.subplots(figsize=(6,6))
poly = matplotlib.patches.Polygon(np2darr, facecolor='C0', edgecolor='none')
ax.add_patch(poly)

patch_minimal = ClickSimplePolygon(poly, ax, bd_list=np.linspace(-0.1,0.1,21) )
patch_minimal.connect()
plt.show()

ax.set_xlim( min(np2darr[:,0])-0.5, max(np2darr[:,0])+0.5 )
ax.set_ylim( min(np2darr[:,1])-0.5, max(np2darr[:,1])+0.5 )

# a single instance of ClickSimplePolygon allows the following methods:
patch_minimal.clear()  # remove the buffered polygon lines from the plot
patch_minimal.compute_cache()  # compute all the buffered polygons with buffer_distance in the list patch.bufferdis_list
patch_minimal.set_buffer( np.linspace(-0.2,0.3,21) )  # change the list patch.bufferdis_list

# the methods below require information of shortest distance between polygons, not useful here
patch_minimal.set_formatxy( '{:.4f}')  
patch_minimal.get_formatxy()           


#%% plot polygons of multiple layers and use ClickSimplePolygon
from QuantumDotDesigner.verification.verification import layermap, uniongdstk, polygonbuffer, ClickSimplePolygon

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)


layer_mapping = {
    (3,4): 3,
    (6,5): 5,
    (22,21): 21,
    (32,31): 31,
    }

polyslist = layermap(qda.main_cell, layer_mapping)


# collect all the instances in a list 
patchlist = []
for i in range(len(polyslist)): 
    polys = polyslist[i]
    polysU, Npolys, NpolysU  = uniongdstk(polys) 
    
    for j in range(len(polysU)):
        p = np.array(polysU[j].points)
    
    
        poly = matplotlib.patches.Polygon(p, facecolor=f'C{i}', edgecolor='none')
        ax.add_patch(poly)
        # give the patchlist and layer id when creating a instance
        patch = ClickSimplePolygon(poly, ax, patchlist,  layerid=polyslist[i][0].layer )
        patch.connect()
        patchlist.append(patch)

ax.set_xlim( -0.5, 0.5 )
ax.set_ylim( -0.5, 0.5 )

# Enalbe zoom and drag functions by giving the axes 'ax' to the following functions 
from QuantumDotDesigner.verification.zoom_axes_drag_logscale_subplots import ZoomPan
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = 1.4)
figPan = zp.pan_factory(ax)
plt.show()
    

# ClickSimplePolygon has an attribute self.patchlist, a list of all the ClickSimplePolygon instances.
# This allows us to access other instances via one instance, in this case 'patchlist[0] '.
# The following methods are achieved in this way:
patchlist[0].print_bypatchid()  # print patch id vs layer id 
patchlist[0].print_bylayerid()  # print layer id vs patch id and numbers of patches in each layer
patchlist[0].printarea_bylayerid( sortarea=True, multiline=False)  # print layer id vs area of patches in each layer
patchlist[0].show_patchid( visible=True, globalscope=True, Nskip=10)  # on the plot mark the patch id
patchlist[0].all_clear()  # for all the patches run .clear(). this remove all buffered polygon lines
patchlist[0].all_compute_cache()  # for all the patches run .compute_cache(). compute all buffered polygons
patchlist[0].all_set_buffer( [0.005,0.015])  # for all the patches run .set_buffer( [0.005,0.015])
patchlist[0].plotall_buffer()  # plot all the buffered polygons
patchlist[0].plot_onelayer_buffer( 31)  # plot all the buffered polygons in layer 31

# the methods below require information of shortest distance between polygons, will gives warning in this cell
patchlist[0].setall_formatxy( '{:.3g}')  
patchlist[0].plotall_shortestdistance( print_flag=False)  
patchlist[0].plot_onelayer_shortestdistance( 21, print_flag=False)  
patchlist[0].printdistance_bylayerid()   


#%% a full example of ClickSimplePolygon, with information of shortest distance between polygons
from QuantumDotDesigner.verification.verification import layermap, uniongdstk, polys_distance, polygonbuffer, ClickSimplePolygon

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
plt.subplots_adjust(right=0.8,top=0.8)
ax_color = fig.add_axes([0.85, 0.1, 0.1, 0.6])


layer_mapping = {
    (3,4): 3,
    (6,5): 5,
    (22,21): 21,
    (32,31): 31,
    }

polyslist = layermap(qda.main_cell, layer_mapping)

colors = plt.cm.jet(np.linspace(0, 1, len(polyslist)))


patchlist = []
for i in range(len(polyslist)): 
    polys = polyslist[i]
    polysU, Npolys, NpolysU  = uniongdstk(polys) 
    
    # compute the shortest distance of a layer
    val_f, xy_f, _ = polys_distance(polysU, Npolys, NpolysU)

    for j in range(len(polysU)):
        p = np.array(polysU[j].points)


        poly = matplotlib.patches.Polygon(p, facecolor=colors[i], edgecolor='none')
        ax.add_patch(poly)
        # give shortest distance information when creating the instance
        patch = ClickSimplePolygon(poly, ax, patchlist, sd_info=(val_f[j], xy_f[j]), layerid=polyslist[i][0].layer )
        patch.connect()
        patchlist.append(patch)





L0 = 1 / (3*len(colors) )
W0 = 0.5
x_legend = 0.1
for i, c in enumerate(colors):
    y_legend = 0.1 +  i* L0*2.5
    ax_color.text(x_legend, y_legend - 0.5*L0, f'layer #{polyslist[i][0].layer }')
    ax_color.add_patch(matplotlib.patches.Rectangle( (x_legend, y_legend), W0, L0,  facecolor=colors[i], edgecolor='none'))

ax_color.set_xticks ([])
ax_color.set_yticks ([])
ax_color.set_xticklabels ([])
ax_color.set_yticklabels ([])
ax_color.set_xlim(0,1)
ax_color.set_ylim(0,1)
ax_color.axis('off')


# bbox = np.array(cell_gdstk.bounding_box())
# bboxX = max(bbox[:,0]) - min(bbox[:,0])
# bboxY = max(bbox[:,1]) - min(bbox[:,1])
# plt.xlim( min(bbox[:,0]), max(bbox[:,0]))
# plt.ylim( min(bbox[:,1]), max(bbox[:,1]))

ax.set_xlim( -0.5, 0.5 )
ax.set_ylim( -0.5, 0.5 )


# Enalbe zoom and drag functions by giving the axes 'ax' to the following functions 
from QuantumDotDesigner.verification.zoom_axes_drag_logscale_subplots import ZoomPan
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = 1.4)
figPan = zp.pan_factory(ax)
plt.show()
    




# pick one of the ClickSimplePolygon instances
   

# A ClickSimplePolygon allow the following methods:
patchlist[0] .clear()
patchlist[0] .compute_cache()
patchlist[0] .set_buffer( [0.005,0.015])
patchlist[0] .set_formatxy( '{:.4f}')  # change the format for printing shortest distance points coordinates
patchlist[0] .get_formatxy()  # get the current format for printing shortest distance points coordinates


# The following methods still works like the previous cell:
patchlist[0].print_bypatchid()  # print patch id vs layer id 
patchlist[0].print_bylayerid()  # print layer id vs patch id and numbers of patches in each layer
patchlist[0].printarea_bylayerid( sortarea=True, multiline=False)  # print layer id vs area of patches in each layer
patchlist[0].show_patchid( visible=True, globalscope=True, Nskip=10)  # on the plot mark the patch id
patchlist[0].all_clear()  # for all the patches run .clear(). this remove all buffered polygon lines
patchlist[0].all_compute_cache()  # for all the patches run .compute_cache(). compute all buffered polygons
patchlist[0].all_set_buffer( [0.005,0.015])  # for all the patches run .set_buffer( [0.005,0.015])
patchlist[0].plotall_buffer()  # plot all the buffered polygons
patchlist[0].plot_onelayer_buffer( 31)  # plot all the buffered polygons in layer 31

# the methods below require information of shortest distance between polygons, will gives warning in this cell
patchlist[0].setall_formatxy( '{:.3g}')   # for all the patches run .set_formatxy( '{:.3g}')
patchlist[0].plotall_shortestdistance( print_flag=False)   # plot shortest distance between polygons in all the layers, and print the information
patchlist[0].plot_onelayer_shortestdistance( 21, print_flag=False)  # plot shortest distance between polygons in layer 21, and print the information
patchlist[0].printdistance_bylayerid()   # print layer id vs shortest distance in each layer

# There are additional options when creating a ClickSimplePolygon instance 

#%%




 