import numpy as np

def fresh_volume(conc, pars):
    
    fresh_list=[]

    for conc_t in conc:
        fresh = 0
        for cell in pars.aquifer_cells:
            if (pars.Lx-pars.L)/pars.dx < cell[-1] and conc_t[cell[0], cell[1], cell[2]] <= 0.35:
                fresh += 1
    
        fresh_list.append(fresh*pars.dx*pars.dz)
    return fresh_list

def fresh_volume_total(conc, pars):
    # fresh volume including aquitard and onshore portion of the aquifer
    fresh_list=[]

    for conc_t in conc:
        fresh = 0
        for cell in pars.aquifer_cells + pars.aquitard_cells:
            if conc_t[cell[0], cell[1], cell[2]] <= 0.35 and cell[2]>pars.x_b/pars.dx:
                fresh += 1
    
        fresh_list.append(fresh*pars.dx*pars.dz)
    return fresh_list

def interface(conc, pars):
    toe_list=[]
    tip_list=[]

    for conc_t in conc:
        for cell in pars.bottom_cells:
            if conc_t[cell[0], cell[1], cell[2]] > 0.35:
                toe_list.append(cell[2]*pars.dx-(pars.Lx-pars.L))
                break
        for cell in pars.top_cells:
            if conc_t[cell[0], cell[1], cell[2]] > 0.35:
                tip_list.append(cell[2]*pars.dx-(pars.Lx-pars.L))
                break

    return toe_list, tip_list

def mixing(conc, pars):
    centroid_list = []
    volume_list = []

    for conc_t in conc:
        mix = 0
        centroid = 0
        for cell in pars.aquifer_cells + pars.aquitard_cells:
            if 31.5 > conc_t[cell[0], cell[1], cell[2]] > 0.35:
                mix += 1
                centroid += cell[2]*pars.dx-(pars.Lx-pars.L)

        volume_list.append(mix*pars.dx*pars.dz)
        try:   
            centroid_list.append(centroid/mix)
        except:
            centroid_list.append(0.0)

    return volume_list, centroid_list

def sgd(conc, qz, pars):
    fresh_flux_list = []
    fresh_centroid_list = []
    sal_flux_list = []
    sal_centroid_list = []

    for conc_t, qz_t in zip(conc, qz):
        f_flux = 0
        f_centroid = 0
        s_flux = 0
        s_centroid = 0
        f_count = 0
        s_count = 0

        for cell in pars.top_cells:
            if qz_t[cell[0], cell[1], cell[2]] > 0:
                if 0 <= conc_t[cell[0], cell[1], cell[2]] < 0.35:
                    f_flux += abs(qz_t[cell[0], cell[1], cell[2]])
                    f_centroid += cell[2]*pars.dx-(pars.Lx-pars.L)
                    f_count += 1
                elif 0.35 <= conc_t[cell[0], cell[1], cell[2]] < 40:
                    s_flux += abs(qz_t[cell[0], cell[1], cell[2]])
                    s_centroid += cell[2]*pars.dx-(pars.Lx-pars.L)
                    s_count += 1

        fresh_flux_list.append(f_flux)
        if f_count != 0:
            fresh_centroid_list.append(f_centroid/f_count)
        else:
            fresh_centroid_list.append(np.nan)
        sal_flux_list.append(s_flux)
        if s_count != 0:
            sal_centroid_list.append(s_centroid/s_count)
        else:
            sal_centroid_list.append(np.nan)

    return fresh_flux_list, fresh_centroid_list, sal_flux_list, sal_centroid_list

def abstracted_flux(conc, qx, pars):
    flux_list = []
    for conc_t, qx_t in zip(conc, qx):
        f_flux = 0

        for cell in pars.aquifer_cells:
            if cell[-1] == 0: 
                if conc_t[cell[0], cell[1], cell[2]] <= 0.35:
                    f_flux += qx_t[cell[0], cell[1], cell[2]]
        flux_list.append(f_flux)
    return flux_list

def paleo_volume(conc, qx, pars, set_x_paleo=None):
    paleo_volumes = [] 
    if set_x_paleo == None:
        for cell in pars.bottom_cells:
            if qx[0][cell[0], cell[1], cell[2]] < 0:
                x_paleo = cell[2]
                break
    else:
        x_paleo = set_x_paleo

    for conc_t in conc:
        if x_paleo == 0:
            paleo_volumes.append(0)
        else:   
            pv = 0           
            for cell in pars.aquifer_cells+pars.aquitard_cells:
                if cell[2]>=x_paleo and conc_t[cell[0], cell[1], cell[2]] <= 0.35:
                    pv += 1
            pv = pv*pars.dx*pars.dz

    paleo_volumes.append(pv)
    if set_x_paleo==None:
        return paleo_volumes, x_paleo
    else:
        return paleo_volumes

def time_constants(concs, head, pars, extent0, zextent0):

    shoreline = pars.bottom_cells[int(pars.x_b/pars.dx)]

    dls = [extent0[2]*pars.dx-pars.x_b]
    for i in range(len(concs)):
        dls.append(extent(concs[i], pars))
    dhh = head[shoreline[0]][shoreline[1]][shoreline[2]] - head[extent0[0]][extent0[1]][extent0[2]]
    dds = [pars.D-((zextent0[0])*pars.dz-(pars.Lz-pars.z0))]
    for i in range(len(concs)):
        zextent = get_z_extent_cell(concs[i], pars)
        dds.append(pars.D-((zextent[0])*pars.dz-(pars.Lz-pars.z0)))
    dhv = head[int(shoreline[0]-pars.H/pars.dz), shoreline[1], shoreline[2]] - head[zextent0[0], zextent0[1], zextent0[2]]

    if dls[-1] <= 0:
        T=500/365*np.argmax([np.array(dls) <= 0])
    else: 
        T=np.min([100/(np.abs((dls[0]-dls[-1]))/dls[-1]), 100/((dds[0]-dds[-1])/dds[-1])])

    return dls[0], dhh, dds[0], dhv, T

def hf(conc, head, pars):
    # calculate hf using method outlined by simmons and post in the book variable density flow
    hf = np.zeros_like(conc)
    for col in range(conc.shape[1]):
        for row in range(conc.shape[0]):
            row = pars.nlay-row-1
            zi = (pars.z_b-pars.dz/2-pars.dz*row) 
            # print(row)
            hf[row, col] = zi + (head[row, col]-zi)*(1000+0.7143*conc[row,col])/1000
            # print((1000+0.7143*conc[row,col])/1000)
    return hf

def extent(conc, pars, return_cell=False):

    x_cell = np.inf
    for cell in pars.aquifer_cells:
        try:
            if conc[cell[0]][cell[1]][cell[2]] >= 0.35 and cell[2] < x_cell:
                x_cell = cell[2]
                extent_cell = cell
        except: 
            pass
    
    if not return_cell:
        return pars.dx*x_cell-pars.x_b
    else:
        return pars.dx*x_cell-pars.x_b, extent_cell
    
def get_z_extent_cell(conc, pars):

    z = np.argmax(np.fliplr(np.array([(36>conc[:, 0, int(pars.x_b/pars.dx)+1]) & (conc[:, 0, int(pars.x_b/pars.dx)+1]>0.35)])))
    if z <= 0:
        z = int((pars.D+pars.z_b)/pars.dz)
    
    return [pars.nlay-z, 0, int(pars.x_b/pars.dx)+1]