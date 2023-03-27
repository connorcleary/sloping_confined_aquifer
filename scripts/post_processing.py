import numpy as np

def fresh_volume(conc, pars):
    
    fresh_list=[]

    for conc_t in conc[0::5]:
        fresh = 0
        for cell in pars.aquifer_cells:
            if (pars.Lx-pars.L)/pars.dx < cell[-1] and conc_t[cell[0], cell[1], cell[2]] <= 0.35:
                fresh += 1
    
        fresh_list.append(fresh*pars.dx*pars.dz)
    return fresh_list

def interface(conc, pars):
    toe_list=[]
    tip_list=[]

    for conc_t in conc[0::5]:
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

    for conc_t in conc[0::5]:
        mix = 0
        centroid = 0
        for cell in pars.aquifer_cells:
            if 31.5 > conc_t[cell[0], cell[1], cell[2]] > 0.35:
                mix += 1
                centroid += cell[2]*pars.dx-(pars.Lx-pars.L)

        volume_list.append(mix*pars.dx*pars.dz)
        centroid_list.append(centroid/mix)

    return volume_list, centroid_list

def sgd(conc, qz, pars):
    fresh_flux_list = []
    fresh_centroid_list = []
    sal_flux_list = []
    sal_centroid_list = []

    for conc_t, qz_t in zip(conc[0::5], qz[0::5]):
        f_flux = 0
        f_centroid = 0
        s_flux = 0
        s_centroid = 0
        f_count = 0
        s_count = 0

        for cell in pars.top_cells:
            if qz_t[cell[0], cell[1], cell[2]] < 0:
                if conc_t[cell[0], cell[1], cell[2]] < 0.35:
                    f_flux += abs(qz_t[cell[0], cell[1], cell[2]])
                    f_centroid += cell[2]*pars.dx-(pars.Lx-pars.L)
                    f_count += 1
                else:
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