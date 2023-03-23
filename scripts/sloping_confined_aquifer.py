from re import L
import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
from pars import ModelParameters, load_parameters
import matplotlib.pyplot as plt


def build_steady_model(pars):
    '''
        A function to build a coastal aquifer model.

        Input: 
            pars: parameters object

        Outputs:
            swt: groundwater flow and transport model
    '''
    
    model_ws = f"./model_files/{pars.name}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    modelname = pars.name
    swt = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=pars.exe_path)

    delc = 10
    delv= pars.dz
    delr=pars.dx

    top = pars.Lz-pars.z0
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)

    ipakcb = 53

    sea_levels = []
    times = []

    if not pars.steady:
        if pars.T_init==0:
            nT = np.min([int(20000*365/pars.T), int((pars.z0-pars.H-pars.D)/pars.v_slr/(pars.T)+1)])
            
        else:
            nT = int(-pars.sea_level/pars.v_slr/pars.T)+1
        perlen = [pars.T_init]*(pars.T_init!=0) + [pars.T for i in range(nT)]
        nT += 1*(pars.T_init!=0)
        nstp = [int(T/pars.dt) for T in perlen]
        tsmult = [1.05 for i in nstp]
    else:
        nT = int(1)
        perlen = [pars.T_init]
        nstp = [int(pars.T_init/pars.dt)]
    
    if pars.T_modern != 0:
        nT += 1
        perlen.append(pars.T_modern) 
        nstp.append(int(1000))
        tsmult.append(1.05)

    dis = flopy.modflow.ModflowDis(
            swt,
            pars.nlay,
            pars.nrow,
            pars.ncol,
            nper=nT,
            itmuni=4,
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=perlen,# perlen, 
            nstp=nstp,
            tsmult=tsmult,
            steady=[False for i in range(nT)]
        )

    # Variables for the BAS package
    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)

    for col in range(pars.ncol):
        ibound[0:int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)),:,col] = 0
        ibound[int(pars.D/pars.dz+pars.H/pars.dz+np.ceil(np.tan(pars.beta)*col*pars.dx/pars.dz)):pars.nlay,:,col] = 0

    strt = top*np.ones((pars.nlay, pars.nrow, pars.ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    # laytyp=np.zeros(pars.nlay)
    # laywet=np.zeros(pars.nlay)
    laytyp=np.ones(pars.nlay)
    laywet=np.ones(pars.nlay)
    
    aquitard_cells = []
    aquifer_cells = []
    top_cells = []
    bottom_cells = []

    for col in range(pars.ncol):
        first = True
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz)): 
            if first:
                first=False
                top_cells.append([lay, 0, col])
            aquitard_cells.append([lay, 0, col])
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+int(pars.D/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz+pars.H/pars.dz)):           
            if lay<pars.nlay:
                aquifer_cells.append([lay, 0, col])
        bottom_cells.append([lay, 0, col])

    hk = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    for cell in aquifer_cells: hk[cell[0], cell[1], cell[2]] = pars.K_aquifer
    for cell in aquitard_cells: hk[cell[0], cell[1], cell[2]] = pars.K_aquitard

    lpf = flopy.modflow.ModflowLpf(swt, hk=hk, vka=pars.anis_aquifer, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet, layvka=1, ss=1e-6)
    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, rclose=1e-5, npcond=1, mxiter=500, iter1=30)

    oc_spd = {} 
    for i in range(len(nstp)):
        for j in range(nstp[i]):
            oc_spd[(i, j)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    offshore_boundary_cells_end = []
    for k in range(0, int(np.floor(pars.D/pars.dz+pars.H/pars.dz))):
            onshore_boundary_cells.append([k, 0, 0])
            offshore_boundary_cells_end.append([k+np.trunc((pars.Lz-pars.D-pars.H)/pars.dz), 0, pars.ncol-1])

    for col in range(int(pars.x_b/pars.dx),pars.ncol):
        offshore_boundary_cells.append([np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz), 0, col])
        if np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz) != np.floor(np.tan(pars.beta)*(col+1)*pars.dx/pars.dz):
           pass #offshore_boundary_cells.append([np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+1, 0, col]) 

    # set up constant head stress period data
    chd_data = {}
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ghb_data = {}
    rch_data = {}
    wel_data = {}

    # sea level at start 
    if not pars.steady:
        if pars.T_init == 0:
            sl0 = -pars.z0+pars.D+pars.H
            sea_levels+=[sl0]
            times+=[-np.min([int(20000*365), np.ceil((pars.z0-pars.H-pars.D)/pars.v_slr/perlen[1])*perlen[1]-1*pars.T])]
        else: 
            sl0 = pars.sea_level
            sea_levels += int(pars.T_init/pars.dt)*[pars.sea_level]
            times+=int(pars.T_init/pars.dt)*[-np.min([int(20000*365), int((-pars.sea_level)/pars.v_slr)-1*pars.T])]
    else:
        sl0 = pars.sea_level

    for i in range(nT-bool(pars.T_init)-2):
        sea_levels += [sl0+pars.v_slr*(i+1)*perlen[i]]
        times += [times[-1]+perlen[i]]

    # f, ax = plt.subplots()
    # ax.pcolormesh(np.flipud(ibound[:,0,:]))
    # for cell in onshore_boundary_cells:
    #     ax.scatter(cell[2], 76-cell[0], color="red")
    # for cell in offshore_boundary_cells:
    #     ax.scatter(cell[2], 76-cell[0], color="blue")
    # plt.show()

    
    for i in range(nT-bool(pars.T_init)-2):
        sea_cells = []
        chd_temp=[]
        ghb_temp = []
        ssm_temp=[]
        rch_temp=[]
        wel_temp=[]
        last_cell = [-1, -1, -1]
        #sea_levels += [sl0+pars.v_slr*i*perlen[i]]
        #times += [times[-1]+perlen[i]]
        for cell in offshore_boundary_cells:
            if pars.Lx-cell[2]*pars.dx<=(pars.Lx-((pars.x_b*np.tan(pars.beta))-sl0)/np.tan(pars.beta))+pars.v_slr*i*perlen[i]/np.tan(pars.beta):
            #if -(pars.dz*(cell[0])+pars.Lz-pars.z0) <= sl0+pars.v_slr*i*perlen[i]:
                #if not 
                if cell[2] != last_cell[2]:
                    h_f = sl0+pars.v_slr*i*perlen[i]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*i*perlen[i])*0.025/1
                else:
                    h_f = sl0+pars.v_slr*i*perlen[i]+((pars.dx*last_cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*i*perlen[i])*0.025/1
                ssm_temp.append([cell[0], cell[1], cell[2], 35.7, 5])
                chd_temp.append([cell[0], cell[1], cell[2], h_f, h_f])
                # ghb_temp.append([cell[0], cell[1], cell[2], h_f, 10000*pars.K_aquitard*pars.anis_aquitard/pars.dz*pars.dx])
                sea_cells.append(cell)
                last_cell = cell
            else:
                pass
               # rch_temp.append(
        
        for cell in offshore_boundary_cells_end[1:-1]:
            h_f = sl0+pars.v_slr*i*perlen[i]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*i*perlen[i])*0.025/1
            ssm_temp.append([cell[0], cell[1], cell[2], 35.7, itype["CHD"]])
            chd_temp.append([cell[0], cell[1], cell[2], h_f, h_f])    
            # ghb_temp.append([cell[0], cell[1], cell[2], h_f, pars.K_aquifer*pars.dx*pars.dz])

     
        for cell in onshore_boundary_cells:
            #if cell[0] >= pars.D/pars.dz:
                #wel_temp.append([cell[0], cell[1], cell[2], 0.0003])
            ssm_temp.append([cell[0], cell[1], cell[2], 0, itype["CHD"]])
            chd_temp.append([cell[0], cell[1], cell[2], pars.Lz-pars.z0, pars.Lz-pars.z0])
            # rch_temp.append([cell[0], cell[1], cell[2],

        ssm_data[i] = ssm_temp
        chd_data[i] = chd_temp
        # ghb_data[i] = ghb_temp
        wel_data[i] = wel_temp
        # rch_data[i] = rch_temp
        
    # rch = flopy.modflow.ModflowRch(
        # model=swt,
        # rech = 0.00001,
        # ipakcb=ipakcb
    # )
    

    if pars.T_modern != 0:
        chd_temp=[]
        chd_temp2=[]
        ghb_temp = []
        ssm_temp=[]
        rch_temp=[]
        wel_temp=[]
        i+=1 
        last_cell = [-1, -1, -1]
        sea_levels += int(pars.T_modern/pars.dt)*[sea_levels[-1]]
        times += [times[-1]+pars.dt_modern*(i+1) for i in range(nstp[-1])]
        for cell in offshore_boundary_cells:
            if pars.Lx-cell[2]*pars.dx<=(pars.Lx-((pars.x_b*np.tan(pars.beta))-sl0)/np.tan(pars.beta))+pars.v_slr*(i-1)*perlen[i-1]/np.tan(pars.beta):
            #if -(pars.dz*(cell[0])+pars.Lz-pars.z0) <= sl0+pars.v_slr*i*perlen[i]:
                #if not 
                if cell[2] != last_cell[2]:
                    h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*(i-1)*perlen[i-1])*0.025/1
                else:
                    h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]+((pars.dx*last_cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*(i-1)*perlen[i-1])*0.025/1
                ssm_temp.append([cell[0], cell[1], cell[2], 35.7, itype["CHD"]])
                chd_temp.append([cell[0], cell[1], cell[2], top, h_f])
                chd_temp2.append([cell[0], cell[1], cell[2], h_f, h_f])
                # ghb_temp.append([cell[0], cell[1], cell[2], h_f, pars.K_aquitard*pars.anis_aquitard/pars.dz*pars.dx])
                sea_cells.append(cell)
                last_cell = cell
            else:
                pass
               # rch_temp.append(
        
        for cell in offshore_boundary_cells_end[:-1]:
            h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*(i-1)*perlen[i-1])*0.025/1
            ssm_temp.append([cell[0], cell[1], cell[2], 35.7, itype["CHD"]])
            chd_temp.append([cell[0], cell[1], cell[2], h_f, h_f])    

            # ghb_temp.append([cell[0], cell[1], cell[2], h_f, pars.K_aquifer*pars.dx*pars.dz])

     
        for cell in onshore_boundary_cells:
            if cell[0]+1 > (pars.Lz-pars.z0-pars.h_modern)/pars.dz:
            # if cell[0] >= (pars.Lz-pars.z0)*pars.dz:
                chd_temp.append([cell[0], cell[1], cell[2], pars.h_modern, pars.h_modern])
                ssm_temp.append([cell[0], cell[1], cell[2], 0.0, itype["CHD"]])

        ssm_data[i] = ssm_temp
        chd_data[i] = chd_temp
        ghb_data[i] = ghb_temp

        #wel_data[i] = wel_temp


    chd = flopy.modflow.ModflowChd(swt, stress_period_data=chd_data, ipakcb = ipakcb)
    # wel  = flopy.modflow.ModflowWel(swt, stress_period_data=wel_data, ipakcb = ipakcb)
    # ghb = flopy.modflow.ModflowGhb(swt, stress_period_data=ghb_data, ipakcb = ipakcb)

    sconc = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    sconc[:, :, -1] = 35.7

    timprs=[1000*(x+1) for x in range(nT)]
    for i in range(int(365*200/100)):
        timprs.append(timprs[-1]+100)

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nstp=nstp,
        tsmult=tsmult,
        nprs=-1,
        # timprs=timprs,
        prsity=pars.porosity, # can change this
        sconc=sconc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=1000
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5
        )
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=pars.alpha_L, trpt=pars.alpha_anisT, trpv=pars.alpha_anisV, dmcoef=pars.diff)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-6)
    mxss = 10*int(len(offshore_boundary_cells+offshore_boundary_cells_end+onshore_boundary_cells))*nT
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=mxss)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=pars.dt,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()
    pars.sea_levels = sea_levels
    pars.times = times
    pars.aquifer_cells = aquifer_cells
    pars.aquitard_cells = aquitard_cells
    pars.top_cells = top_cells
    pars.bottom_cells = bottom_cells
    pars.save_parameters()
    return swt 

def run_model(swt):
    """
        A function to run the seawat model

        Inputs: 
            swt: model object
        Outputs:
            None
    """
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")

def build_modern_model(name):

    model_ws_old = f"./model_files/{name}"
    ucnobj = bf.UcnFile(os.path.join(model_ws_old, "MT3D001.UCN"))
    cbbobj = bf.CellBudgetFile(os.path.join(model_ws_old, f'{name}.cbc'))
    headobj = bf.HeadFile(os.path.join(model_ws_old, f'{name}.hds'))
 
    # get head and concentration data
    concentration = ucnobj.get_alldata()[-2]
    head = headobj.get_alldata()[-2]

    pars = load_parameters(name)
    pars.name = pars.name + "_modern"
    model_ws = f"./model_files/{pars.name}"

    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    modelname = pars.name
    swt = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=pars.exe_path)

    delc = 10
    delv= pars.dz
    delr=pars.dx

    top = pars.Lz-pars.z0
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)

    ipakcb = 53

    nT  = 1
    perlen = 365*200
    nstp = 365*2
    tsmult = 1

    dis = flopy.modflow.ModflowDis(
            swt,
            pars.nlay,
            pars.nrow,
            pars.ncol,
            nper=nT,
            itmuni=4,
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=perlen,# perlen, 
            nstp=nstp,
            tsmult=tsmult,
            steady=[False for i in range(nT)]
        )

    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)

    for col in range(pars.ncol):
        ibound[0:int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)),:,col] = 0
        ibound[int(pars.D/pars.dz+pars.H/pars.dz+np.ceil(np.tan(pars.beta)*col*pars.dx/pars.dz)):pars.nlay,:,col] = 0

    strt = head[:, 0, :]
    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.ones(pars.nlay)
    laywet=np.ones(pars.nlay)
    
    aquitard_cells = []
    aquifer_cells = []
    top_cells = []
    bottom_cells = []

    for col in range(pars.ncol):
        first = True
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz)): 
            if first:
                first=False
                top_cells.append([lay, 0, col])
            aquitard_cells.append([lay, 0, col])
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+int(pars.D/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz+pars.H/pars.dz)):           
            if lay<pars.nlay:
                aquifer_cells.append([lay, 0, col])
        bottom_cells.append([lay, 0, col])

    hk = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    for cell in aquifer_cells: hk[cell[0], cell[1], cell[2]] = pars.K_aquifer
    for cell in aquitard_cells: hk[cell[0], cell[1], cell[2]] = pars.K_aquitard

    lpf = flopy.modflow.ModflowLpf(swt, hk=hk, vka=pars.anis_aquifer, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet, layvka=1, ss=1e-6)
    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, rclose=1e-5, npcond=1, mxiter=500, iter1=30)

    oc_spd = {} 
    i = 0 
    for j in range(int(365*2)):
            oc_spd[(i, j)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    offshore_boundary_cells_end = []
    for k in range(0, int(np.floor(pars.D/pars.dz+pars.H/pars.dz))):
            onshore_boundary_cells.append([k, 0, 0])
            offshore_boundary_cells_end.append([k+np.trunc((pars.Lz-pars.D-pars.H)/pars.dz), 0, pars.ncol-1])

    for col in range(int(pars.x_b/pars.dx),pars.ncol):
        offshore_boundary_cells.append([np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz), 0, col])
       
    chd_data = {}
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ghb_data = {}
    rch_data = {}
    wel_data = {}
    
    chd_temp=[]
    chd_temp2=[]
    ghb_temp = []
    ssm_temp=[]
    rch_temp=[]
    wel_temp=[]
    i+=1 
    last_cell = [-1, -1, -1]
    for cell in offshore_boundary_cells:
        if pars.Lx-cell[2]*pars.dx<=pars.Lx-((pars.x_b*np.tan(pars.beta))):
            if cell[2] != last_cell[2]:
                h_f = (pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))*0.025/1
            else:
                h_f = (pars.dx*last_cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))*0.025/1
            ssm_temp.append([cell[0], cell[1], cell[2], 35.7, itype["CHD"]])
            chd_temp.append([cell[0], cell[1], cell[2], h_f, h_f])
            chd_temp2.append([cell[0], cell[1], cell[2], h_f, h_f])
            last_cell = cell
    
    for cell in offshore_boundary_cells_end[:-1]:
        h_f = (pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))*0.025/1
        ssm_temp.append([cell[0], cell[1], cell[2], 35.7, itype["CHD"]])
        chd_temp.append([cell[0], cell[1], cell[2], h_f, h_f])    
    
    for cell in onshore_boundary_cells:
        if cell[0]+1 > (pars.Lz-pars.z0-pars.h_modern)/pars.dz:
            chd_temp.append([cell[0], cell[1], cell[2], pars.h_modern, pars.h_modern])
            ssm_temp.append([cell[0], cell[1], cell[2], 0.0, itype["CHD"]])

    ssm_data[i] = ssm_temp
    chd_data[i] = chd_temp

    chd = flopy.modflow.ModflowChd(swt, stress_period_data=chd_data, ipakcb = ipakcb)

    sconc = np.abs(concentration[:, 0, :])

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nstp=365*2,
        nprs=-1,
        prsity=pars.porosity, # can change this
        sconc=sconc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=0
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5
        )
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=pars.alpha_L, trpt=pars.alpha_anisT, trpv=pars.alpha_anisV, dmcoef=pars.diff)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-6)
    mxss = 10*int(len(offshore_boundary_cells+offshore_boundary_cells_end+onshore_boundary_cells))*nT
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=mxss)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=pars.dt,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()
    pars.save_parameters()
    return swt 

def extract_results(name):
    """
        Open model results from binary files

        Inputs:
            name: name of model/realization/scenario
        Outputs:
            head: head matrix [nstp, nlay, nrow, ncol]
            qx: longitudinal flux matrix [nstp, nlay, nrow, ncol]
            qy: transverse flux matrix matrix [nstp, nlay, nrow, ncol]
            qz: vertical flux matrix matrix [nstp, nlay, nrow, ncol]
            concentration: concentration matrix [nstp, nlay, nrow, ncol]
    """
    pars = load_parameters(name)
    name = pars.name
    model_ws = f"./model_files/{name}"
    #nstp = pars.perlen/pars.dt

    # open binary files
    ucnobj = bf.UcnFile(os.path.join(model_ws, "MT3D001.UCN"))
    cbbobj = bf.CellBudgetFile(os.path.join(model_ws, f'{name}.cbc'))
    headobj = bf.HeadFile(os.path.join(model_ws, f'{name}.hds'))
 
    # get head and concentration data
    concentration = ucnobj.get_alldata()[:]
    head = headobj.get_alldata()[:]
    
    # select every n items
    times = ucnobj.get_times()
    concentration = concentration

    qx = np.zeros_like(concentration)
    qy = np.zeros_like(concentration)
    qz = np.zeros_like(concentration)

    # get fluxes
    for t in range(qx.shape[0]):
            qx[t] = cbbobj.get_data(text="flow right face")[0]
            if pars.nrow > 1:
                qy[t] = cbbobj.get_data(text="flow front face")[0]
            qz[t] = cbbobj.get_data(text="flow lower face")[0]
    # qx[0] = cbbobj.get_data(text="flow right face", totim=times[-1])[0]
    # #qy[0] = cbbobj.get_data(text="flow front face", totim=times[-1])[0]
    # qz[0] = cbbobj.get_data(text="flow lower face", totim=times[-1])[0]

    save_results(name, concentration, head, qx, qy, qz)
    return concentration, head, qx, qy, qz


def save_results(name, concentration, head, qx, qy, qz):
    """
        Save extracted results to a .npy file

        Inputs:
            name: model name
            concentration, head etc. : numpy arrays of model outputs
        Outputs:
            None
    """
    ws = os.path.join(f'./results/{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)

    with open(os.path.join(ws, f"qx.npy"), 'wb') as f: np.save(f, np.array(qx))
    with open(os.path.join(ws, f"qy.npy"), 'wb') as f: np.save(f, np.array(qy))
    with open(os.path.join(ws, f"qz.npy"), 'wb') as f: np.save(f, np.array(qz))
    with open(os.path.join(ws, f"head.npy"), 'wb') as f: np.save(f, np.array(head))
    with open(os.path.join(ws, f"concentration.npy"), 'wb') as f: np.save(f, np.array(concentration))


def load_results(name):
    """
        Load extracted results from .npy files

        Inputs:
            name: name of the model
        Outputs:
            concentration, head... : numpy matrices of results
    """
    ws = os.path.join(f'./results/{name}')

    with open(os.path.join(ws, f"qx.npy"), 'rb') as f: qx = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qy.npy"), 'rb') as f: qy = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qz.npy"), 'rb') as f: qz = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"head.npy"), 'rb') as f: head = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"concentration.npy"), 'rb') as f: concentration = np.load(f, allow_pickle=True)

    return concentration, head, qx, qy, qz, 
