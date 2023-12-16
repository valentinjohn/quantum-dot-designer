# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 16:33:45 2023

@author: cawang
"""

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import gdstk

def layermap(cell_gdstk, layer_mapping):
    '''
    input
        cell_gdstk: gdstk.Cell object
        layer_mapping: dictionary
    return 
        polyslist: nested list. list of lists of polygons
        
    example:
        layer_mapping = {
            (6,5): 5,
            (22,21): 21,
            (32,31): 31,
            }  
        This dictionary means the polygons in 
        cell_gdstk, layer = 5 and 6 polygons 
        will become layer = 5 polygons in polyslist
        Similarly, layer 22 and 21 become 21, 
        layer 32 and 31 become 31
    '''
    polyslist = []
    for key, value in layer_mapping.items():
        polys = []
        layers = []
        for item in cell_gdstk.get_polygons():
            if item.layer in key:
                item_new = item.copy()
                item_new.layer = value
                polys += [ item_new  ] 
        polyslist += [ polys ]
    return polyslist

def uniongdstk(polys, precision=1e-12):
    # Notice: the return polygon is an gdstk.Polygon, with gdstk.Polygon().layer = 0
    # which is likely not the same as the layer number of the input
    # precision is the precision for the union operation gdstk.boolen('or')
    polysU = polys.copy()
    Npolys = len(polys)
    NpolysU = len(polysU)
    NpolysU_before = NpolysU
    
    count_loop = 0
    while(count_loop==0 or NpolysU < NpolysU_before):
        count_loop += 1
        NpolysU = len(polysU)
        NpolysU_before = NpolysU
       
        for i in range(Npolys):
            for j in range(i+1,Npolys):
    
                if j>=NpolysU:
                    break
                polyunion = gdstk.boolean(polysU[i], polysU[j], 'or', precision=precision )
                
                if len(polyunion) == 1:
                    polysU[i] = polyunion[0]
                    polysU.pop(j)
                    NpolysU = len(polysU)
            if i>=NpolysU:
                break
    return polysU, Npolys, NpolysU

def singlelayerconnection(polys):
    # polys: list of gdstk.Polygon
    if isinstance(polys,list):
        for i in range(len(polys)):
            if isinstance(polys[i], gdstk.Polygon):
                if polys[i].layer != polys[0].layer:
                    ValueError('layers should be the same')
            else:
                ValueError('all elements in the list should be gdstk.Polygon object')
    else:
        ValueError('input should be a list')
    polysU, Npolys, NpolysU  = uniongdstk(polys)              

    # print(len(polysU))
    polyarea = [p.area() for p in polys]
    polyUarea = [p.area() for p in polysU]
    
    
    return polysU, Npolys, NpolysU, polyarea, polyUarea





def cellconnection(cell_gdstk, layer_mapping, printflag=True):
    layer_mapping_inv = {}
    for key, value in layer_mapping.items():
        layer_mapping_inv[value] = key
    
    polyslist = layermap(cell_gdstk, layer_mapping)
    # polyslist = []
    # for key, value in layer_mapping.items():
    #     polys = []
    #     layers = []
    #     for item in cell_gdstk.get_polygons():
    #         if item.layer in key:
    #             item_new = item.copy()
    #             item_new.layer = value
    #             polys += [ item_new  ] 
    
    #     polyslist += [ polys ]
    
    results = {}
    for i in range(len(polyslist)):
        polys = polyslist[i]
        polysU, Npolys, NpolysU, polyarea, polyUarea = singlelayerconnection(polys)
        for p in polysU:
            p.layer = polys[0].layer
            
        polyarea_sorted = sorted(polyarea)
        polyUarea_sorted = sorted(polyUarea)
        if printflag:
            print('number of polygons in layer_'+ str(polys[0].layer) + '  (includes layer ' + str(layer_mapping_inv[polys[0].layer])+ ')' + ' \n'\
                          + '* Before union operations: ' + str(Npolys) + '\n'  \
                          +  f'area minimum = {min(polyarea):.3g},  maximum = {max(polyarea):.3g}' + ' \n'\
                          + '* After union operations: ' +  str(NpolysU) + '\n' \
                          +  f'area minimum = {min(polyUarea):.3g},  maximum = {max(polyUarea):.3g}' + '\n' \
                          +  f'sorted area = ' + ''.join(['{:.2g}, ']*len(polyUarea_sorted)).format(*polyUarea_sorted) + '\n' )
                
        results['layer_'+ str(polys[0].layer)] = dict(polysU=polysU, Npolys=Npolys, NpolysU=NpolysU, polyarea=polyarea, polyUarea=polyUarea)
    return results



def linesection_point_distance(linepoint1, linepoint2, pointxs, pointys ):
    # x0 = np.array(allx)
    # y0 = np.array(ally)
    for arg_in in [linepoint1, linepoint2]:
        if not (isinstance(arg_in, np.ndarray) or isinstance(arg_in, list) or isinstance(arg_in, tuple)):
            raise ValueError('The inputs \'linepoint1\' and \'linepoint2\', each of them should be either numpy array, list or tuple')
        elif np.array(arg_in).shape != (2,):
            raise ValueError('The inputs \'linepoint1\' and \'linepoint2\', each of them should have shape = (2,)')
                

    if isinstance(pointxs, np.ndarray):
        if len(np.array(pointxs).shape) != 1:
            raise ValueError('The input numpy array \'pointxs\'should be 1d')
    elif isinstance(pointxs, list) or isinstance(pointxs, tuple):
        pointxs = np.array(pointxs)
        if len(np.array(pointxs).shape) != 1:
            raise ValueError('The input numpy array \'pointxs\' should be 1d')
    elif np.isscalar(pointxs):
        pointxs = np.array([pointxs])            
    else:
        raise ValueError('The input \'pointxs\' should be either numpy array, list or tuple')
        
    if isinstance(pointys, np.ndarray):
        if len(np.array(pointys).shape) != 1:
            raise ValueError('The input numpy array \'pointys\'should be 1d')
    elif isinstance(pointys, list) or isinstance(pointys, tuple):
        pointys = np.array(pointys)
        if len(np.array(pointys).shape) != 1:
            raise ValueError('The input numpy array \'pointys\' should be 1d')
    elif np.isscalar(pointys):
        pointys = np.array([pointys])      
    else:
        raise ValueError('The input \'pointys\' should be either numpy array, list or tuple')        
    x1, y1 = linepoint1
    x2, y2 = linepoint2

  
    x0 = pointxs
    y0 = pointys
    
    a = y2-y1
    b = -x2+x1
    c = -y2*x1 + x2*y1
    
    
    x = (b*(b*x0 - a*y0) - a*c )/(a**2 + b**2)
    y = (a*(-b*x0 + a*y0) - b*c)/(a**2 + b**2)
    xy = np.vstack((x, y))
    
    dis = np.abs(a*x0 + b*y0 + c)/np.sqrt(a**2 + b**2)

    
    # dis12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    dis01 = np.sqrt((x0-x1)**2 + (y0-y1)**2)
    dis02 = np.sqrt((x0-x2)**2 + (y0-y2)**2)
    
    maskx = np.logical_or(np.logical_and(x<x1, x<x2 ), np.logical_and(x>x1, x>x2 ) )
    masky = np.logical_or(np.logical_and(y<y1, y<y2 ), np.logical_and(y>y1, y>y2 ) )
    
    mask = np.logical_or(maskx, masky)
    dis_endpts = np.min( np.vstack((dis01,dis02))  , axis = 0 )
    dis_new = np.array(dis)
    dis_new[mask] = dis_endpts[mask]
    # print(mask)
    C1 = np.logical_and(mask, dis01 < dis02)
    C2 = np.logical_and(mask, dis01 > dis02)
    x_new = np.array(x)
    # print(linepoint1)
    x_new[C1] = x1
    x_new[C2] = x2
    
    y_new = np.array(y)
    y_new[C1] = y1
    y_new[C2] = y2   
    xy_new = np.vstack((x_new, y_new))

    shortestpoint_idx = np.zeros((len(mask),), dtype=int)  #  # shortest point is on the line
    shortestpoint_idx[C1] = 1  # shortest point is the end point: linepoint1
    shortestpoint_idx[C2] = 2  # shortest point is the end point: linepoint2
    return dis, xy.T, dis_new, xy_new.T, shortestpoint_idx


def polys_distance(polysU, Npolys, NpolysU):
   
    report = dict(vals=np.zeros((NpolysU,NpolysU)), 
                  xs0=np.zeros((NpolysU,NpolysU)), 
                  ys0=np.zeros((NpolysU,NpolysU)), 
                  xs1=np.zeros((NpolysU,NpolysU)), 
                  ys1=np.zeros((NpolysU,NpolysU)), )
    for pA in range(NpolysU):
        for pB in range(NpolysU):
            if pB != pA:
                polyA = polysU[pA]
                polyB = polysU[pB]
                vals = np.zeros((polyA.points.shape[0],))
                xs0 = np.zeros((polyA.points.shape[0],))
                ys0 = np.zeros((polyA.points.shape[0],))
                xs1 = np.zeros((polyA.points.shape[0],))
                ys1 = np.zeros((polyA.points.shape[0],))
                for ipA in range(-1,polyA.points.shape[0]-1):
                # for ipA in [0]:
                    x1, y1 = polyA.points[ipA]
                    x2, y2 = polyA.points[ipA+1]
            
                    
                    x0 = polyB.points[:,0]
                    y0 = polyB.points[:,1]
                    _, xy, dis_new, xy_new, _ = linesection_point_distance( (x1, y1),(x2, y2), x0, y0 )
     
                    
                    idx_min = np.argmin(dis_new)
                    vals[ipA] = dis_new[idx_min] 
                    xs0[ipA] = x0[idx_min]
                    ys0[ipA] = y0[idx_min]
                    xs1[ipA] = xy_new[idx_min,0]
                    ys1[ipA] = xy_new[idx_min,1]
                    # print(dis_new)
                    # print(x0)
                    
                idx_min = np.argmin(vals)
                report['vals'][pA,pB] =  vals[idx_min] 
                report['xs0'][pA,pB] =  xs0[idx_min]
                report['ys0'][pA,pB] =  ys0[idx_min]
                report['xs1'][pA,pB] =  xs1[idx_min]
                report['ys1'][pA,pB] =  ys1[idx_min]
                # print(report['vals'])
                
                
    report_new = dict(vals=np.zeros((NpolysU,)), 
                  xs0=np.zeros((NpolysU,)), 
                  ys0=np.zeros((NpolysU,)), 
                  xs1=np.zeros((NpolysU,)), 
                  ys1=np.zeros((NpolysU,)), )
    
    for i in range(NpolysU):
        X = []
        Y = []
        for j in range(NpolysU):
            if i!= j:
                Y += [ report['vals'][i,j] , report['vals'][j,i] ]
                X += [ [i,j], [j,i] ]
        idx = np.argmin(Y)
        ii, jj = X[idx]
        report_new['vals'][i] = Y[idx]
        report_new['xs0'][i] = report['xs0'][ii, jj]
        report_new['ys0'][i] = report['ys0'][ii, jj]
        report_new['xs1'][i] = report['xs1'][ii, jj]
        report_new['ys1'][i] = report['ys1'][ii, jj]
        
        vals_f = []
        xy_f = []
    for i in range(NpolysU):
        vals_f  += [  report_new['vals'][i]  ]
        xy_f += [ np.array([ [report_new['xs1'][i], report_new['ys1'][i]],
                             [report_new['xs0'][i], report_new['ys0'][i]],  ]) ] 
    return vals_f, xy_f, report_new



def intersect_check(p, precision=1e-15):
    # p is 2d numpy array

    # p = np.array([[0, 0], [0, 4], [4,4], [4, 0], [1, 1], [2, 1], [2, 2], [1, 2]])
    # p = np.array([[0, 0], [0, 4], [4,4], [4, 0], [1, 1], [1, 2], [2, 2], [2, 1]])
    # print(intersect_check(p))

    if not isinstance(p, np.ndarray):
        raise ValueError('The input should be a numpy array')
    elif len(p.shape) != 2:
        raise ValueError('The input numpy array should be 2d')
    elif p.shape[1] != 2:
        raise ValueError('The input numpy array should be n by 2')
        
    Lp = p.shape[0]
    report = dict(intersect=[], intersectidx=[])
    for i in [-1] + list(range(Lp-1)):
        x1, y1 = p[i,:]
        x2, y2 = p[i+1,:]
        a = y2-y1
        b = -x2+x1
        c = -y2*x1 + x2*y1
        for j in list(range(i+2, Lp-1)):
            if j+1 != Lp+i:
                x3, y3 = p[j,:]
                x4, y4 = p[j+1,:]
                d = y4-y3
                e = -x4+x3
                f = -y4*x3 + x4*y3            
                
                C1 = d*x1 + e*y1 + f
                C2 = d*x2 + e*y2 + f
                C3 = a*x3 + b*y3 + c
                C4 = a*x4 + b*y4 + c
                
                # if C1*C2 < 0 and C3*C4 <0:
                if C1*C2 < 0 and C3*C4 <0 and abs(C1*C2)>precision and abs(C3*C4)>precision:
                    # print(C1, C2, C3, C4)
                    report['intersect'] += [ [ [p[i,:], p[i+1,:]], 
                                               [p[j,:], p[j+1,:]]  ] ]
                    report['intersectidx'] += [ [ [i, i+1], 
                                                  [j, j+1]  ] ]                    
    if  len(report['intersect'])==0:
        intersect_flag = False
    else:
        intersect_flag = True
    return intersect_flag, report

                   
                
def outward_normvec(p):
    # p = np.array([[1,-0.1], [1, -1], [-1, -1], [-1,1], [1, 1], [1,0.1], [3,0.1], [3,1], [4,1], [4,-1],[3,-1],[3,-0.1],]) 
    if not isinstance(p, np.ndarray):
        raise ValueError('The input should be a numpy array')
    elif len(p.shape) != 2:
        raise ValueError('The input numpy array should be 2d')
    elif p.shape[1] != 2:
        raise ValueError('The input numpy array should be n by 2')
        
    e = np.diff( p, axis=0, append=p[0:1,:])
    en = e/np.sqrt(e[:,0:1]**2+e[:,1:2]**2)
    
    ex = en[:,0]
    ey = en[:,1]
    sin = ex*np.roll(ey,-1) - np.roll(ex,-1)*ey
    
    eclockwise = np.sum(np.arcsin(sin)) < 0
    
    vec_out = np.roll(en,1,axis=1)
    if intersect_check(p)[0]:
        clockwiseness = 'edges intersect'
        vec_out = None
    elif eclockwise:
        clockwiseness = 'clockwise'
        vec_out[:,0] = -vec_out[:,0]
    else:
        clockwiseness = 'counter clockwise'
        vec_out[:,1] = -vec_out[:,1]
        
    return clockwiseness, vec_out


def edgeoverlap_check(p, precision=1e-14):
    # p is 2d numpy array    
    # p = np.array([[0, 0], [0, 4], [4,4], [4, 0], [1, 1], [1, 2], [2, 2], [2, 1]])
    # p = np.array([[0, 0], [0, 4], [4,4], [4, 0], [0, 0], [1, 1], [2, 1], [2, 2], [1, 2], [1,1]])
    if not isinstance(p, np.ndarray):
        raise ValueError('The input should be a numpy array')
    elif len(p.shape) != 2:
        raise ValueError('The input numpy array should be 2d')
    elif p.shape[1] != 2:
        raise ValueError('The input numpy array should be n by 2')  
        
    if intersect_check(p)[0]:
        raise ValueError('The input polygon should not have intersect edges. It should be a simple planar polygon.')
            
    Lp = p.shape[0] 
    report = dict(overlap3pts=[], overlap3ptsidx=[], overlap4pts=[], overlap4ptsidx=[])
    for i in [-1] + list(range(Lp-1)):
        x0, y0 = p[i-1,:]
        x1, y1 = p[i,:]
        x2, y2 = p[i+1,:]        
        # linesection_point_distance( (x0, y0), (x1, y1), x2, y2)
        _, _, dis_new, _, _ = linesection_point_distance( (x0, y0), (x1, y1), x2, y2)
        if dis_new < precision:
            report['overlap3pts'] += [ [ [x0, y0], 
                                         [x1, y1],
                                         [x2, y2]] ]
            report['overlap3ptsidx'] += [ [ i-1, i, i+1, ] ]            
        
       

    for i in [-1] + list(range(Lp-1)):
        for j in list(range(i+2, Lp-1)):
            if j+1 != Lp+i:
                
                if np.all(p[i,:] == p[j+1,:]) and np.all(p[i+1,:] == p[j,:]):
                    report['overlap4pts'] += [ [ [p[i,:], p[i+1,:]], 
                                                 [p[j,:], p[j+1,:]]  ] ]
                    report['overlap4ptsidx'] += [ [ [i, i+1], 
                                                   [j, j+1]  ] ] 
                    
    overlap3pts_flag = False if  len(report['overlap3pts'])==0 else True
    overlap4pts_flag = False if  len(report['overlap4pts'])==0 else True

    return overlap3pts_flag, overlap4pts_flag, report



def polygonbuffer(p, distance, corner_seg_dtheta, debugflag=False, boolean_precision=1e-5):
    '''
    expand or shrink a given polygon
    example:        
        p = np.array([[0, 0], [0, 3], [3,3], [3, 0], [0, 0], [0.4, 0.8], [1, 1], [2, 1], [2, 2], [1, 2], [1,1], [0.4, 0.8]])
        polynew = polygonbuffer(p, 0.02, 170*np.pi/180)
    '''
    if not isinstance(p, np.ndarray):
        raise ValueError('The input should be a numpy array')
    elif len(p.shape) != 2:
        raise ValueError('The input numpy array should be 2d')
    elif p.shape[1] != 2:
        raise ValueError('The input numpy array should be n by 2')


    intersect_check_result = intersect_check(p, precision=1e-15)
    if intersect_check_result[0]:
        print('The input polygon contains intersect edges:')
        for i in range(len(intersect_check_result[1]['intersect'])):
            print(f'#{i} idx = ' + str(intersect_check_result[1]['intersectidx'][i]))
            print(f'#{i} coordinates = ' + str(intersect_check_result[1]['intersect'][i]))
            # print('\n')
        raise ValueError('The input polygon should not have intersect edges. It should be a simple planar polygon.')
  
    overlap3pts_flag, overlap4pts_flag, reportoverlap = edgeoverlap_check(p, precision=1e-14)
    if overlap3pts_flag:
        # row_id = []
        # print('The input polygon contains trivial inner points, the script is removing them:')
        # for i in range(len(reportoverlap['overlap3pts'])):
        #     print(f'#{i} idx = ' + str(reportoverlap['overlap3ptsidx'][i]) )
        #     print(f'#{i} coordinates = ' + str(reportoverlap['overlap3pts'][i])  )
        #     row_id += [ reportoverlap['overlap3ptsidx'][i][1] ]
        #     ddd = p[reportoverlap['overlap3ptsidx'][i][0],:] - p[reportoverlap['overlap3ptsidx'][i][2],:]
        #     if np.sqrt(np.sum(ddd**2)) < 1e-12:
        #         row_id += [ reportoverlap['overlap3ptsidx'][i][2] ]
        # p = np.delete(p, row_id, axis=0)
        print('The input polygon contains trivial inner points (one vertex very close to an edge).')
        
        


    Lp = p.shape[0]
    clockwiseness, vec_out = outward_normvec(p)
        
    overlap4ptsidx = []
    if overlap4pts_flag:
        for i in range(len(reportoverlap['overlap4ptsidx'])):
            overlap4ptsidx += [ reportoverlap['overlap4ptsidx'][i][0][0]    ] 
            overlap4ptsidx += [ reportoverlap['overlap4ptsidx'][i][1][0]    ] 


    # print(reportoverlap['overlap4ptsidx'])  
    # print(overlap4ptsidx)
    # print(Lp)

    addrecs1 = []
    for i in list(range(Lp-1)) + [-1]:
        if i in overlap4ptsidx:
            continue
        x0, y0 = p[i,:]
        x1, y1 = p[i+1,:]
        x2, y2 = p[i+1,:] + distance*vec_out[i]
        x3, y3 = p[i,:]  + distance*vec_out[i]  
        addrec1 = np.array([[x0, y0],
                            [x1, y1],
                            [x2, y2],
                            [x3, y3],]) 
        addrecs1 += [addrec1]
        

    # print(vec_out)
    addarcs1 = []
    # print(reportoverlap['overlap4ptsidx'])
    # print(overlap4ptsidx)
    for i in list(range(Lp-1)) + [-1]:
        if (i in overlap4ptsidx) or \
           (i in np.array(overlap4ptsidx)+1) or \
           (i in np.array(overlap4ptsidx)+1-Lp) :
            # print(i)
            continue        
        
        vec_out1 = vec_out[i,:]
        vec_out2 = vec_out[i-1,:]
        x0, y0 = p[i,:]
        x1, y1 = p[i,:] + distance*vec_out1
        x2, y2 = p[i,:] + distance*vec_out2
        
        Mrot = np.array([[ vec_out1[0], vec_out1[1] ],
                         [-vec_out1[1], vec_out1[0],]])
        
        vec_out2_rot = Mrot @ vec_out2

        theta1 = np.arctan2(vec_out1[1], vec_out1[0])
        theta2 = theta1 + np.arctan2(vec_out2_rot[1], vec_out2_rot[0])

        n_arc = int(np.floor(abs((theta2-theta1)/corner_seg_dtheta)))
        thetas = np.linspace(theta1, theta2, 2+n_arc)[1:-1]
        # print(theta1, theta2)
        # print(thetas)
        p_arc = [ [x0, y0], [x1, y1] ]
        for j in range(len(thetas)):
            x_arc = x0 + distance*np.cos(thetas[j])
            y_arc = y0 + distance*np.sin(thetas[j])
            p_arc += [ [x_arc,y_arc] ]
        p_arc += [ [x2,y2] ]
        addarcs1 += [np.array(p_arc)]
        
    # print(addarcs1)
    # print(overlap4ptsidx)
    for i in range(len(reportoverlap['overlap4ptsidx'])):
        idx = reportoverlap['overlap4ptsidx'][i]
        # print(idx)
        # idx0 = idx[0][0]-1
        # idx1 = idx[0][0]
        # idx2 = idx[1][1] 
        for idx0, idx1, idx2 in [ (idx[0][0]-1,idx[0][0],idx[1][1]), 
                                  (idx[1][0]-1,idx[1][0],idx[0][1]) ]:
            # print(idx0)
            # print(idx0 % Lp, overlap4ptsidx)
            if idx0%Lp  in np.array(overlap4ptsidx)%Lp:
                continue
            # print('*** ')
            vec_out1 = vec_out[idx0,:]
            vec_out2 = vec_out[idx2,:]
            x0, y0 = p[idx1,:]
            x1, y1 = p[idx1,:] + distance*vec_out1
            x2, y2 = p[idx1,:] + distance*vec_out2
            
            Mrot = np.array([[ vec_out1[0], vec_out1[1] ],
                             [-vec_out1[1], vec_out1[0],]])
            
            vec_out2_rot = Mrot @ vec_out2
    
            theta1 = np.arctan2(vec_out1[1], vec_out1[0])
            theta2 = theta1 + np.arctan2(vec_out2_rot[1], vec_out2_rot[0])
    
            n_arc = int(np.floor(abs((theta2-theta1)/corner_seg_dtheta)))
            thetas = np.linspace(theta1, theta2, 2+n_arc)[1:-1]
            # print(theta1, theta2)
            # print(thetas)
            p_arc = [ [x0, y0], [x1, y1] ]
            for j in range(len(thetas)):
                x_arc = x0 + distance*np.cos(thetas[j])
                y_arc = y0 + distance*np.sin(thetas[j])
                p_arc += [ [x_arc,y_arc] ]
            p_arc += [ [x2,y2] ]
            addarcs1 += [np.array(p_arc)]  
            # print(p_arc)
        # x2, y2 = p[idx,:]
    # print(addarcs1)
    # print('2'*10)
    polynew = [gdstk.Polygon(p)]
    rec_debug = []
    
    # for i in range(Lp):
    for i in range(len(addrecs1)):
        if distance > 0:
            polynew_1 = gdstk.boolean(polynew, gdstk.Polygon(addrecs1[i]), "or", precision=boolean_precision)
            polynew = polynew_1 
            rec_debug += [gdstk.Polygon(addrecs1[i])]
        elif distance < 0:
            polynew_1 = gdstk.boolean(polynew, gdstk.Polygon(addrecs1[i]), "not", precision=boolean_precision)       
            polynew = polynew_1 
            rec_debug += [gdstk.Polygon(addrecs1[i])]  
    
    arc_debug = []
    # for i in range(Lp):
    for i in range(len(addarcs1)):
        # break
        if distance > 0:
            polynew_1 = gdstk.boolean(polynew, gdstk.Polygon(addarcs1[i]), "or", precision=boolean_precision)
            polynew = polynew_1 
            arc_debug += [gdstk.Polygon(addarcs1[i])]
        elif distance < 0:
            polynew_1 = gdstk.boolean(polynew, gdstk.Polygon(addarcs1[i]), "not", precision=boolean_precision)
            polynew = polynew_1 
            arc_debug += [gdstk.Polygon(addarcs1[i])]
    pointslist = [x.points for x in polynew]
    if debugflag:
        return pointslist, rec_debug, arc_debug
    else:
        return pointslist




class ClickSimplePolygon:
    '''
    This is a class that can plot buffered polygons upon mouse left/right click  
    on matplotlib.patches.Polygon instances in a plot represented by self.ax.
    
    minimal example:
        np2darr = np.array([[0, 0], [0, 1], [1,1], [1, 0]])
        fig, ax = plt.subplots(figsize=(6,6))
        poly = matplotlib.patches.Polygon(np2darr, facecolor=colors[i], edgecolor='none')
        ax.add_patch(poly)
        patch = ClickSimplePolygon(poly, ax, bd_list=np.linspace(-0.2,0.3,21) )
        patch.connect()
        plt.show()
        
    There are additional options when creating a ClickSimplePolygon instance: 
        titlefont : size of the plot title
        plotlines_arg : dictionary. for plotting the buffered polygons
        formatxy : a string. the format of printing the coordinates of shortest distance points     
    important attributes:
        self.poly0 : matplotlib.patches.Polygon instance. The buffered polygons are computed based on the coordinates stored in this instance
        self.press : a variable responsible for clicking events
        self.ax : axes created at beginning
        self.patchlist : list of ClickSimplePolygon instances
        self.patchid : an integer labelling this instance among all the other instances. The rule currently using is the position in the list self.patchlist
        
        self.polyid : this variable is currently unused
        self.lines : a list storing the plotted matplotlib.lines.Line2D instances
        self.bufferdis_list : a list of float numbers, buffer distance
        self.bufferdis_list_origin : the self.bufferdis_list when the instance is created
        self.bufferdis_id : an integer indicates the buffer distance in self.bufferdis_list
        self.buffercache : a list storing polygon coordinates that have been plotted, to speed up the display afterward
        self.shortestdistance_info = sd_info
        self.layerid : an integer represents layer id. It can be assigned when creating the instance. 
        self.button1_count : a variable responsible for plotting when after mouse lift button click
        
        self.titlefont = titlefont
        self.plotlines_arg : a float number describing the corners of buffered polygons
        
        self.text_patchid : a list of matplotlib.pyplot.Text instances after running self.show_patchid()
        self.Nskip : an integer describing the number of vertices should skip when running self.show_patchid()
        
        self.formatxy : a string. the current format for printing shortest distance points coordinates
    '''

    
    def __init__(self, poly0, ax, patchlist=[], bd_list=[0.01,0.02,0.03], corner_arc_angle=10*np.pi/180, 
                 sd_info=None, layerid=None,
                 titlefont=15, plotlines_arg={'color':'r','linewidth':1}, formatxy='{:.3f}'):
        self.poly0 = poly0
        self.press = None
        self.ax = ax
        self.patchlist = patchlist
        self.patchid = len(self.patchlist)
        
        self.polyid = len(ax.patches) - 1  # this variable is currently unused
        self.lines = []
        self.bufferdis_list = bd_list.copy()
        self.bufferdis_list_origin = bd_list.copy()
        self.bufferdis_id = -1
        self.buffercache = [None for _ in self.bufferdis_list]
        self.corner_arc_angle = corner_arc_angle  
        self.shortestdistance_info = sd_info
        self.layerid = layerid
        self.button1_count = 0 if sd_info is None else 2
        
        self.titlefont = titlefont
        self.plotlines_arg = plotlines_arg
        
        self.text_patchid = []
        self.Nskip = None
        
        self.formatxy = formatxy
    def connect(self):
        #'connect to all the events we need'
        self.cidpress = self.poly0.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)

    def compute_polygonbuffer(self, p):
        return polygonbuffer(np.array(p), self.bufferdis_list[self.bufferdis_id], self.corner_arc_angle, boolean_precision=1e-6)
    
    def plotlines(self, xs, ys):
        # return self.ax.plot( xs, ys, 'r', linewidth=1)[0]
        return self.ax.plot( xs, ys, **self.plotlines_arg)[0]
    
    def clear(self):
        if len( self.lines) != 0:
            for i in range(len(self.lines)):
                self.ax.lines.remove(self.lines[i])
            self.lines = []   
        self.bufferdis_id = -1
        self.button1_count = 0 if self.shortestdistance_info is None else 2
        

    
    def check_patchlist(self):      
        if len(self.patchlist) == 0:
            print('The argument \'patchlist\' is not assigned. Give a list containing the same object to itself during initialization.')
            return  True  
        else:
            return False
        
    def compute_cache(self):
        bufferdis_id = self.bufferdis_id
        for i in range(len(self.buffercache)):
            if self.buffercache[i] is None:
                self.bufferdis_id = i
                p = self.poly0.xy[:-1,:]
                polynew = self.compute_polygonbuffer(p)
                self.buffercache[i] = polynew.copy()        
        self.bufferdis_id = bufferdis_id
        
    def print_bypatchid(self):
        if self.check_patchlist(): return
        for i, patch in enumerate(self.patchlist):
            str1, str2 = '', ''
            if patch.layerid is not None:
                str1 = f'in layer #{patch.layerid}, '
            if patch.shortestdistance_info is not None:
                str2 = f'shortest distance = {patch.shortestdistance_info[0]:.3g}'
            print(f'patch #{patch.patchid} ' + str1 + str2  )
    
    def set_formatxy(self, formatxy):
        self.formatxy = formatxy
    def get_formatxy(self):
        return self.formatxy
            
    def prtxy(self,x,y):
        return (self.formatxy + ', ' + self.formatxy).format(x,y)

    def setall_formatxy(self, formatxy):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            patch.set_formatxy(formatxy)
    
    def print_bylayerid(self):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            if patch.layerid is None:
                print('At least one polygon is not specified layer id.')
                return   
            

        info = {}
        for i, patch in enumerate(self.patchlist):          
            if patch.layerid not in info.keys():
                info[patch.layerid] = []
            info[patch.layerid] += [ patch.patchid ]

        for key, value in info.items():
            strvalue = str(value).replace('[','').replace(']','')
            print(f'Layer #{key} contains {len(value)} patches: patchid = #' + strvalue)


    def printdistance_bylayerid(self):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            if patch.layerid is None:
                print('At least one polygon is not specified layer id.')
    
        info_dis = {}
        info_xy = {}
        if True not in [patch.shortestdistance_info is None for patch in self.patchlist]:
            
            for i, patch in enumerate(self.patchlist):
                
                x0, y0 = patch.shortestdistance_info[1][0,:]
                x1, y1 = patch.shortestdistance_info[1][1,:]              
                if patch.layerid not in info_dis.keys():
                    info_dis[patch.layerid] = []
                    info_xy[patch.layerid] = []
                
                info_dis[patch.layerid] += [ patch.shortestdistance_info[0] ] 
                info_xy[patch.layerid] += [ [x0, y0, x1, y1, patch] ]
            
            for key, value in info_dis.items():
                idxmin = np.argmin(value)
                x0, y0, x1, y1, patch = info_xy[key][idxmin]
                strvalue = str(value).replace('[','').replace(']','')
                print(f'Layer #{key} shortest distance = {value[idxmin]:.3g} \n' \
                      + f'--> from patch #{patch.patchid} coordinates (x0,y0)=(' + patch.prtxy(x0,y0) + ') to (x1,y1)=(' + patch.prtxy(x0,y0) +  ')' )                
    
    def printarea_bylayerid(self, sortarea=True, multiline=False):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            if patch.layerid is None:
                print('At least one polygon is not specified layer id.')
                return   
            
        info_area = {}
        info_id = {}     
        for i, patch in enumerate(self.patchlist):
            if patch.layerid not in info_area.keys():
                info_area[patch.layerid] = []
                info_id[patch.layerid] = []
            x = patch.poly0.xy[:-1,0]
            y = patch.poly0.xy[:-1,1]
            info_area[patch.layerid] += [ 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1))) ] 
            #  Shoelace formula  https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
            info_id[patch.layerid] += [ patch.patchid ]


        for key, areas in info_area.items():
            idxs = np.argsort(areas) if sortarea else range(len(areas))
            areasort = [ info_area[key][iii] for iii in idxs ]
            patchid = [ info_id[key][iii] for iii in idxs ]
            print(f'Layer #{key}, (patchid, area) = ')
            prtstr = ''
            for i in range(len(areasort)):
                prtstr += f'({patchid[i]},{areasort[i]:.2g}) '
                if multiline: prtstr += '\n'
            print(prtstr + '\n' )


    def show_patchid(self, visible=True, globalscope=True, Nskip=10):
        if Nskip is None: return
        if self.check_patchlist(): return

        for i, patch in enumerate(self.patchlist):
            if len(patch.text_patchid) == 0:
                arr = patch.poly0.xy[:-1,:]
                for i in range(0,arr.shape[0], Nskip+1):
                    x, y = arr[i,:]
                    patch.text_patchid += [  self.ax.text(x, y, str(patch.patchid) , clip_on=True) ]
                if len(patch.text_patchid) > 0:
                    patch.Nskip = Nskip
            elif Nskip != patch.Nskip:
                for text in patch.text_patchid:
                    text.remove()
                patch.text_patchid = []
                arr = patch.poly0.xy[:-1,:]
                for i in range(0,arr.shape[0], Nskip+1):
                    x, y = arr[i,:]
                    patch.text_patchid += [  self.ax.text(x, y, str(patch.patchid) , clip_on=True) ]                
                patch.Nskip = Nskip
                    
        if globalscope:
            for i, patch in enumerate(self.patchlist):
                for text in patch.text_patchid:
                    text.set_visible(visible)   
        else:
            for i, patch in enumerate(self.patchlist):
                xL, xR = self.ax.get_xlim()
                yL, yR = self.ax.get_ylim()            
                for text in patch.text_patchid:
                    x, y = text.get_position()
                    if xL<x<xR and yL<y<yR:
                        text.set_visible(visible)
         
                    
            #         text.set_visible(visible)
    def all_clear(self):
        if self.check_patchlist(): return
        for i in range(len(self.patchlist)):
            self.patchlist[i].clear()

    def all_compute_cache(self):
        if self.check_patchlist(): return
        for i in range(len(self.patchlist)):
            self.patchlist[i].compute_cache()    

    def plotall_shortestdistance(self, print_flag=False):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            if patch.shortestdistance_info is None or patch.layerid is None:
                print('At least one polygon is not specified layer id or shortest distance information.')
                return        
        dislist = []
        info = []
        for i in range(len(self.patchlist)):
            self.patchlist[i].button1_count = 2
            self.patchlist[i].press_event_button_1() 
            x0, y0 = self.patchlist[i].shortestdistance_info[1][0,:]
            x1, y1 = self.patchlist[i].shortestdistance_info[1][1,:]            
            
            dislist += [ self.patchlist[i].shortestdistance_info[0] ]
            info += [ dict(patchid=self.patchlist[i].patchid, patch=patch, xy=(x0, y0, x1, y1)  )  ]
            if print_flag:
                patch = self.patchlist[i]
                print(f'Patch #{self.patchlist[i].patchid} in patchlist')
                print(f'distance = {self.patchlist[i].shortestdistance_info[0]:.3g}' )
                print(f'coordinates (x0,y0)=(' + patch.prtxy(x0,y0) + ') to (x1,y1)=(' + patch.prtxy(x0,y0) +  ') \n'  )   
                
        x0, y0, x1, y1 = info[np.argmin(dislist)]['xy']
        minpatch = info[np.argmin(dislist)]['patchid']
        patch = info[np.argmin(dislist)]['patch']
        titlestr = 'Plot all shortest distance in the same layer \n' \
            + f' shortest distance = {min(dislist):.3g} from patch #{minpatch} \n' \
                + f'(x0,y0)=(' + patch.prtxy(x0,y0) + ') to (x1,y1)=(' + patch.prtxy(x0,y0) +  ') \n'
        self.ax.set_title(titlestr, fontsize=self.titlefont)  

    def plot_onelayer_shortestdistance(self, layerid, print_flag=False):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            if patch.layerid is None or patch.shortestdistance_info is None:
                print('At least one polygon is not specified layer id or shortest distance information.')
                return
        dislist = []
        info = []
        for patch in self.patchlist:
            if layerid == patch.layerid:
                patch.button1_count = 2
                patch.press_event_button_1()
                x0, y0 = patch.shortestdistance_info[1][0,:]
                x1, y1 = patch.shortestdistance_info[1][1,:]            
                
                dislist += [ patch.shortestdistance_info[0] ]
                info += [ dict(patchid=patch.patchid, patch=patch, xy=(x0, y0, x1, y1)  )  ]
                if print_flag:
                    print(f'Patch #{patch.patchid} in patchlist')
                    print(f'distance = {patch.shortestdistance_info[0]:.3g}' )
                    print(f'coordinates (x0,y0)=(' + patch.prtxy(x0,y0) + ') to (x1,y1)=(' + patch.prtxy(x0,y0) +  ') \n'  )   
            else:
                patch.clear()
               
        if len(dislist) != 0:
            x0, y0, x1, y1 = info[np.argmin(dislist)]['xy']
            minpatch = info[np.argmin(dislist)]['patchid']
            patch = info[np.argmin(dislist)]['patch']
            titlestr = f'Plot shortest distance in the layer #{layerid} \n' \
                + f' shortest distance = {min(dislist):.3g} from patch #{minpatch} \n' \
                    + f'(x0,y0)=(' + patch.prtxy(x0,y0) + ') to (x1,y1)=(' + patch.prtxy(x0,y0) +  ') \n'             
            self.ax.set_title(titlestr, fontsize=self.titlefont)  
        else:
            self.ax.set_title(f'The given layerid #{layerid} not found', fontsize=self.titlefont)   
            
    def set_buffer(self, new_bd_list):
        if isinstance(new_bd_list, list):
            self.bufferdis_list = new_bd_list
            self.buffercache = [None for _ in self.bufferdis_list]            
    def all_set_buffer(self, new_bd_list):
        if self.check_patchlist(): return
        for i in range(len(self.patchlist)):
            self.patchlist[i].set_buffer(new_bd_list) 

    def plotall_buffer(self):
        if self.check_patchlist(): return

        bufferdis_list_all = self.patchlist[0].bufferdis_list_origin.copy()
        bufferdis_id_all = (self.patchlist[0].bufferdis_id + 1) % len(bufferdis_list_all)

        for patch in self.patchlist:
            # tempvar = patch.bufferdis_list.copy()
            # patch.bufferdis_list = bufferdis_list_all
            if len(patch.bufferdis_list)==len(bufferdis_list_all) and np.max(np.abs(np.array(patch.bufferdis_list)-np.array(bufferdis_list_all))) < 1e-12:
                pass
            else:
                patch.set_buffer(bufferdis_list_all)
            patch.bufferdis_id = bufferdis_id_all
            patch.button1_count = 0
            patch.press_event_button_1()   
            # patch.bufferdis_list = tempvar.copy()
        self.ax.set_title(f'Plot all buffer distance = {self.bufferdis_list[bufferdis_id_all]:.3g}', fontsize=self.patchlist[0].titlefont)
    


    def plot_onelayer_buffer(self, layerid):
        if self.check_patchlist(): return
        for patch in self.patchlist:
            if patch.layerid is None:
                print('At least one polygon is not specified layer id.')
                return
        findlayer = False
        for patch in self.patchlist:
            if layerid == patch.layerid:
                if findlayer is False:
                    bufferdis_list_this = patch.bufferdis_list_origin.copy()
                    bufferdis_id_this = (patch.bufferdis_id + 1) % len(bufferdis_list_this)                    
                # patch.bufferdis_list = bufferdis_list_this
                if len(patch.bufferdis_list)==len(bufferdis_list_this) and np.max(np.abs(np.array(patch.bufferdis_list)-np.array(bufferdis_list_this))) < 1e-12:
                    pass
                else:
                    patch.set_buffer(bufferdis_list_this)
                patch.bufferdis_id = bufferdis_id_this
                patch.button1_count = 0
                patch.press_event_button_1()         
                findlayer = True
            else:
                patch.clear()
               
        if findlayer:           
            self.ax.set_title(f'Plot layer #{layerid} buffer distance = {bufferdis_list_this[bufferdis_id_this]:.3g}', fontsize=self.titlefont)  
        else:
            self.ax.set_title(f'The given layerid #{layerid} not found', fontsize=self.titlefont)
        

        
    def verticesmean(self,pp):
        return np.mean(pp[:-1,:],axis=0)
    
    def press_event_button_1(self):
        if len( self.lines) != 0:
            for i in range(len(self.lines)):
                self.ax.lines.remove(self.lines[i])
            self.lines = []   

        if self.button1_count == 0:
            if self.buffercache[self.bufferdis_id] is None:
                p = self.poly0.xy[:-1,:]
                polynew = self.compute_polygonbuffer(p)
                self.buffercache[self.bufferdis_id] = polynew.copy()
            else:
                polynew = self.buffercache[self.bufferdis_id]
                
            for i in range(len(polynew)):
                xyplot = np.vstack((polynew[i], polynew[i][0:1,:]))
                self.lines  += [ self.plotlines(xyplot[:,0], xyplot[:,1]) ]    
            
            titlestr_begin = '' if len(self.patchlist)==0 else f'#{self.patchid} in patchlist \n'
            titlestr =  titlestr_begin + f'buffer distance = {self.bufferdis_list[self.bufferdis_id]:.3g}'
            self.ax.set_title(titlestr, fontsize=self.titlefont)                    
        elif self.button1_count == 1:
            self.ax.set_title('')
        elif self.button1_count == 2:
            if self.shortestdistance_info is not None and self.layerid is not None:
                self.lines  += [  self.plotlines( self.shortestdistance_info[1][:,0],  self.shortestdistance_info[1][:,1]) ] 
                titlestr = f'#{self.patchid} in patchlist' + '\n' + f'shortest distance in the same layer #{self.layerid} = {self.shortestdistance_info[0]:.3g}'
                self.ax.set_title(titlestr, fontsize=self.titlefont)                
            else:
                print('Should provide these two values when initialize the object: layerid and shortestdistance_info ')
        self.poly0.figure.canvas.draw()
        Nmode = 2 if self.shortestdistance_info is None else 3
        self.button1_count = (self.button1_count + 1) % Nmode 

    def press_event_button_3(self):
        self.button1_count = 1
        
        self.bufferdis_id = (self.bufferdis_id+1) % len(self.bufferdis_list)
        if len( self.lines) != 0:
            for i in range(len(self.lines)):
                self.ax.lines.remove(self.lines[i])
            self.lines = []
            
        if self.buffercache[self.bufferdis_id] is None:
            p = self.poly0.xy[:-1,:]
            polynew = self.compute_polygonbuffer(p)
            self.buffercache[self.bufferdis_id] = polynew.copy()
        else:
            polynew = self.buffercache[self.bufferdis_id]
        
        for i in range(len(polynew)):
            xyplot = np.vstack((polynew[i], polynew[i][0:1,:]))                
            self.lines  += [ self.plotlines(xyplot[:,0], xyplot[:,1])  ]              
        titlestr_begin = '' if len(self.patchlist)==0 else f'#{self.patchid} in patchlist \n'
        titlestr =  titlestr_begin + f'buffer distance = {self.bufferdis_list[self.bufferdis_id]:.3g}'
        self.ax.set_title(titlestr, fontsize=self.titlefont)
        self.poly0.figure.canvas.draw()
        
        
    def on_press(self, event):
        #'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.poly0.axes: return
        
        contains, attrd = self.poly0.contains(event)
        if not contains: return
        # print('event contains', self.verticesmean(self.poly0.xy)  )
       
        if event.button==1: 
            self.press_event_button_1()
            
        elif event.button==3:
            self.press_event_button_3()

    def disconnect(self):
        self.poly0.figure.canvas.mpl_disconnect(self.cidpress)

