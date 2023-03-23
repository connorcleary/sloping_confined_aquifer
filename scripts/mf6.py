import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
from pars import ModelParameters, load_parameters
import matplotlib.pyplot as plt

def build_steady_model(pars):

    mf6exe = pars.exe_path
    length_units = "m"
    time_units = "days"

    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97

    modelname = pars.name 

    ws = f"./model_files/{pars.name}"
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation(
        sim_name=pars.name, sim_ws=ws, exe_name=mf6exe
    )
    time_spd = []
    if not pars.steady:
        if pars.T_init==0:
            nper = np.min([int(20000*365/pars.T), int((pars.z0-pars.H-pars.D)/pars.v_slr/(pars.T)+1)])
            
        else:
            nper = int(-pars.sea_level/pars.v_slr/pars.T)+1
        perlen = [pars.T_init]*(pars.T_init!=0) + [pars.T for i in range(nper)]
        for i in range(nper):
            time_spd.append((pars.T, int(perlen[i]/pars.dt), 1.01))
        nper += 1*(pars.T_init!=0)
        nstp = [int(T/pars.dt) for T in perlen]
        tsmult = [1.01 for i in nstp]
    else:
        nper = int(1)
        perlen = [pars.T_init]
        nstp = [int(pars.T_init/pars.dt)]
    
    if pars.T_modern != 0:
        nper += 1
        perlen.append(pars.T_modern) 
        nstp.append(int(1000))
        tsmult.append(1.05)
        time_spd.append((pars.T_modern, 1000, 1.01))

    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=pars.name, save_flows=True)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )

    sim.register_ims_package(ims, [gwf.name])

    delc = 10
    delv= pars.dz
    delr=pars.dx
    top = pars.Lz-pars.z0
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)

    idomain = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)

    for col in range(pars.ncol):
        idomain[0:int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)),:,col] = 0
        idomain[int(pars.D/pars.dz+pars.H/pars.dz+np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)):pars.nlay,:,col] = 0

    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=pars.nlay,
        nrow=pars.nrow,
        ncol=pars.ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain
    )

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


    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=1,
        k=hk,
        k33=0.01,
    )
    strt = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    strt[:, :, -1] = 35.7
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    onshore_boundary_cells = []
    offshore_boundary_cells = []
    offshore_boundary_cells_end = []
    for k in range(0, int(np.floor(pars.D/pars.dz+pars.H/pars.dz))):
            onshore_boundary_cells.append([int(k), 0, 0])
            offshore_boundary_cells_end.append([int(k+np.trunc((pars.Lz-pars.D-pars.H)/pars.dz)), 0, pars.ncol-1])

    for col in range(int(pars.x_b/pars.dx),pars.ncol):
        offshore_boundary_cells.append([int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)), 0, col])
 
    sea_levels = []
    times=[]
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

    for i in range(nper-bool(pars.T_init)-2):
        sea_levels += [sl0+pars.v_slr*(i+1)*perlen[i]]
        times += [times[-1]+perlen[i]]

    ghb_data = {}
    ghbcond_v = pars.K_aquitard/pars.anis_aquitard * delr * delc / (0.5 * delv)
    ghbcond_ha = pars.K_aquifer * delv * delc / (0.5 * delr)
    ghbcond_hb = pars.K_aquitard * delv * delc / (0.5 * delr)

    for i in range(nper-bool(pars.T_init)-2):
        sea_cells = []
        ghb_temp = []
        last_cell = [-1, -1, -1]

        for cell in offshore_boundary_cells:
            if pars.Lx-cell[2]*pars.dx<=(pars.Lx-((pars.x_b*np.tan(pars.beta))-sl0)/np.tan(pars.beta))+pars.v_slr*i*perlen[i]/np.tan(pars.beta):
                if cell[2] != last_cell[2]:
                    h_f = sl0+pars.v_slr*i*perlen[i]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*i*perlen[i])*0.025/1
                    # h_f = sl0+pars.v_slr*i*perlen[i]
                else:
                    h_f = sl0+pars.v_slr*i*perlen[i]+((pars.dx*last_cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*i*perlen[i])*0.025/1
                    # h_f = sl0+pars.v_slr*i*perlen[i]
                
                ghb_temp.append([(cell[0], cell[1], cell[2]), h_f, ghbcond_v, 35.0])
                sea_cells.append(cell)
                last_cell = cell

        
        for cell in offshore_boundary_cells_end[1:-1]:
            h_f = sl0+pars.v_slr*i*perlen[i]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*i*perlen[i])*0.025/1
            # h_f = sl0+pars.v_slr*i*perlen[i]
            if cell[0] <= pars.Lz-pars.H:
                ghbcond=ghbcond_hb
            else:
                ghbcond=ghbcond_ha
            ghb_temp.append([(cell[0], cell[1], cell[2]), h_f, ghbcond, 35.0])
     
        for cell in onshore_boundary_cells:
            if cell[0] <= pars.Lz-pars.H:
                ghbcond=ghbcond_hb
            else:
                ghbcond=ghbcond_ha
            ghb_temp.append([(cell[0], cell[1], cell[2]), pars.Lz-pars.z0, ghbcond, 0.0])

        ghb_data[i] = ghb_temp

    if pars.T_modern != 0:
        ghb_temp = []
        i+=1 
        last_cell = [-1, -1, -1]
        sea_levels += int(pars.T_modern/pars.dt)*[sea_levels[-1]]
        times += [times[-1]+pars.dt_modern*(i+1) for i in range(nstp[-1])]
        for cell in offshore_boundary_cells:
            if pars.Lx-cell[2]*pars.dx<=(pars.Lx-((pars.x_b*np.tan(pars.beta))-sl0)/np.tan(pars.beta))+pars.v_slr*(i-1)*perlen[i-1]/np.tan(pars.beta):
                if cell[2] != last_cell[2]:
                    h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*(i-1)*perlen[i-1])*0.025/1
                    #h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]
                else:
                    h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]+((pars.dx*last_cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*(i-1)*perlen[i-1])*0.025/1
                    #h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]
                ghb_temp.append([(cell[0], cell[1], cell[2]), h_f, ghbcond_v, 35.0])
                sea_cells.append(cell)
                last_cell = cell

        for cell in offshore_boundary_cells_end[:-1]:
            h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]+((pars.dx*cell[2]*np.tan(pars.beta)-(pars.Lz-pars.z0))-sl0-pars.v_slr*(i-1)*perlen[i-1])*0.025/1
            # h_f = sl0+pars.v_slr*(i-1)*perlen[i-1]  
            if cell[0] <= pars.Lz-pars.H:
                ghbcond=ghbcond_hb
            else:
                ghbcond=ghbcond_ha
            ghb_temp.append([(cell[0], cell[1], cell[2]), h_f, ghbcond, 35.0])

     
        for cell in onshore_boundary_cells:
            if cell[0]+1 > (pars.Lz-pars.z0-pars.h_modern)/pars.dz:
                if cell[0] <= pars.Lz-pars.H:
                    ghbcond=ghbcond_hb
                else:
                    ghbcond=ghbcond_ha
                ghb_temp.append([(cell[0], cell[1], cell[2]), pars.h_modern, ghbcond, 0.0])

        ghb_data[i] = ghb_temp


    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_data,
        pname="GHB-1",
        auxiliary="CONCENTRATION",
    )

    head_filerecord = "{}.hds".format(pars.name)
    budget_filerecord = "{}.bud".format(pars.name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )

    sim.register_ims_package(imsgwt, [gwt.name])

    flopy.mf6.ModflowGwtdis(
        gwt,
        idomain=idomain,
        length_units=length_units,
        nlay=pars.nlay,
        nrow=pars.nrow,
        ncol=pars.ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    flopy.mf6.ModflowGwtmst(gwt, porosity=pars.porosity)

    flopy.mf6.ModflowGwtic(gwt, strt=strt)

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=pars.diff, alh=pars.alpha_L, ath1=pars.alpha_L*pars.alpha_anisT, alv=pars.alpha_L*pars.alpha_anisV)

    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
    ]

    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord="{}.cbc".format(gwt.name),
        concentration_filerecord="{}.ucn".format(gwt.name),
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    pars.sea_levels = sea_levels
    pars.times = times
    pars.aquifer_cells = aquifer_cells
    pars.aquitard_cells = aquitard_cells
    pars.top_cells = top_cells
    pars.bottom_cells = bottom_cells
    pars.save_parameters()

    return sim

def run_model(swt):
    swt.write_simulation()
    success, buff = swt.run_simulation(silent=True)
    if not success:
            print(buff)

    return success

def get_results(name):
    pars = load_parameters(name)

    ws = f"./model_files/{name}"
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation.load(
        sim_ws=ws,
        exe_name=pars.exe_path,
        verbosity_level=0,
    )
    gwf = sim.get_model(name)
    gwt = sim.get_model("trans")
    
    conc = gwt.output.concentration().get_alldata()
    head = gwf.output.head().get_alldata()
    bud = gwf.output.budget()
    
    times = bud.get_times()
    qx = np.zeros_like(conc)
    qy = np.zeros_like(conc)
    qz = np.zeros_like(conc)

    for i in range(qx.shape[0]):
        qx[i], qy[i], qz[i] = flopy.utils.postprocessing.get_specific_discharge(
            bud.get_data(text="DATA-SPDIS", totim=times[i])[0],
            gwf,
        )

    return conc, head, qx, qz, times-(times[-1]-365*200)*np.ones_like(times)