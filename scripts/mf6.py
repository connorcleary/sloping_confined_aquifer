import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
from pars import ModelParameters, load_parameters
import matplotlib.pyplot as plt
import pdb
import shutil
import scipy

def build_steady_model(pars, just_for_pars=False):

    spl_time_to_sl, spl_sl_to_time = get_time_series_spline()
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

    # transgression parameters
    sl0 = np.min([pars.z0-pars.D-pars.H,120])
    t0 = np.min([np.abs(spl_sl_to_time(-(pars.z0-pars.D-pars.H))),24])
    #pdb.set_trace()
    nper_trans = int(t0*365000/(pars.T))
    # time parameters
    nper = nper_trans+2
    perlen = [pars.T_init] + nper_trans*[pars.T] + [pars.T_modern] 
    
    for p in perlen:
    	time_spd.append((p, p/pars.dt, 1))

    time_spd[0] = (pars.T_init, 5*pars.T_init/pars.dt, 1)
    time_spd[-1] = (p, p/10, 1)
    	
    
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

    delc = 1
    delv= pars.dz
    delr=pars.dx
    top = pars.Lz-pars.z0
    #pdb.set_trace()
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)

    idomain = np.zeros((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)

    aquitard_cells = []
    aquifer_cells = []
    top_cells = []
    bottom_cells = []
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    offshore_boundary_cells_end = []

    for col in range(pars.ncol):
        first = True
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz)): 
            if first:
                first=False
                if col > int(pars.x_b/pars.dx):
                    offshore_boundary_cells.append([lay, 0, col])
            aquitard_cells.append([lay, 0, col])
            idomain[lay, 0, col]=1
            if col == 0:
                onshore_boundary_cells.append([lay, 0, col])
            elif col == pars.ncol-1 and not first:
                offshore_boundary_cells_end.append([lay, 0, col])
        first = True
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+int(pars.D/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz+pars.H/pars.dz)):           
            if first:
                first=False
                top_cells.append([lay, 0, col])
            if lay<pars.nlay:
                aquifer_cells.append([lay, 0, col])
            idomain[lay, 0, col]=1
            if col == 0:
                onshore_boundary_cells.append([lay, 0, col])
            elif col == pars.ncol-1 and not first:
                offshore_boundary_cells_end.append([lay, 0, col])
        bottom_cells.append([lay, 0, col])
    if offshore_boundary_cells_end[-1] != [pars.nlay-1, 0, pars.ncol-1]:
        aquifer_cells.append([pars.nlay-1, 0, pars.ncol-1])
        idomain[pars.nlay-1, 0, pars.ncol-1]=1
        offshore_boundary_cells_end.append([pars.nlay-1, 0, pars.ncol-1])
            
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
    hk = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    k33 = np.zeros_like(hk)

    for cell in aquifer_cells: 
        hk[cell[0], cell[1], cell[2]] = pars.K_aquifer
        k33[cell[0], cell[1], cell[2]] = pars.K_aquifer/pars.anis_aquifer

    for cell in aquitard_cells: 
        hk[cell[0], cell[1], cell[2]] = pars.K_aquitard
        k33[cell[0], cell[1], cell[2]] = pars.K_aquitard/pars.anis_aquitard

    # pdb.set_trace()
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        k33=k33,
        rewet_record="REWET WETFCT 1.0 IWETIT 1 IHDWET 0",
        wetdry=1
    )
    strt = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    # strt[:, :, -1] = 35.7
    flopy.mf6.ModflowGwfic(gwf, strt=8.5)

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)
    
    sea_levels = [-sl0] + [spl_time_to_sl(t0-(pars.T/365000)*(i+1)) for i in range(nper_trans)] + [0]
    #pdb.set_trace()
    ghb_data = {}
    chd_data= {}
    cnc_data= {}
    ghbcond_v = 1e5# pars.K_aquitard/pars.anis_aquitard * delr * delc / (delv)
    ghbcond_ha = 1e5 # pars.K_aquifer * delv * delc / (delr)
    ghbcond_hb = 1e5# pars.K_aquitard * delv * delc / (delr)

    for i in range(nper):
    
        ghb_temp = []
        chd_temp = []
        cnc_temp = []

        if i >= 1:
            for cell in offshore_boundary_cells:  
            	# if pars.Lx-cell[2]*pars.dx<=(pars.Lx-((pars.x_b*np.tan(pars.beta))-sea_levels[0])/np.tan(pars.beta))+(sl0 + spl_time_to_sl(t0-(pars.T/(365000))*(i+1)))/np.tan(pars.beta):
                if cell[2] >=(pars.z_b-sea_levels[i])/np.tan(pars.beta)/pars.dx:
                    chd_temp.append([(cell[0], cell[1], cell[2]), sea_levels[i], 35.0])
                    cnc_temp.append([(cell[0], cell[1], cell[2]), 35.0])
        
        for cell in offshore_boundary_cells_end[1:]:
            if -(sea_levels[i]) <= cell[0]*pars.dz-(pars.Lz-pars.z0):
                if cell[0] <= pars.Lz-pars.H:
                    ghbcond=ghbcond_hb
                else:
                    ghbcond=ghbcond_ha
                chd_temp.append([(cell[0], cell[1], cell[2]), sea_levels[i], 35.0])
                cnc_temp.append([(cell[0], cell[1], cell[2]), 35.0])
        
        for cell in onshore_boundary_cells:
            
            if cell[0] <= pars.Lz-pars.H:
                ghbcond=ghbcond_hb
            else:
                ghbcond=ghbcond_ha
            if i == nper-1:
                h_onshore = pars.h_modern
            else:   
                h_onshore = pars.gradient*(pars.Lz-pars.z0)/pars.x_b*(pars.Lx-(pars.z0-pars.H-pars.D+sea_levels[i])/np.tan(pars.beta))+sea_levels[i]
            if cell[0]*pars.dz-pars.Lz+pars.z0 >= -h_onshore and cell[0]*pars.dz>=pars.D:
            	chd_temp.append([(cell[0], cell[1], cell[2]), h_onshore, 0.0])
            	
        chd_data[i] = chd_temp
        ghb_data[i] = ghb_temp
        cnc_data[i] = cnc_temp

        print(h_onshore)
    # pdb.set_trace()
    # 
    # flopy.mf6.ModflowGwfghb(
    #     gwf,
    #     stress_period_data=ghb_data,
    #     pname="GHB-1",
    #     auxiliary=["CONCENTRATION","DENSITY"],
    # )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )
    
    # flopy.mf6.ModflowGwfrcha(gwf, recharge=0.00000001)
    
    
    ws_results = f"./model_files/{pars.name}"
    if not os.path.exists(ws_results):
        os.makedirs(ws_results)

    head_filerecord = f"{pars.name}.hds"
    budget_filerecord = f"{pars.name}.bud"
    if pars.outputs == "all":
        saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")], nper_trans: [("HEAD", "LAST"), 
                            ("BUDGET", "LAST")], nper_trans+1: [("HEAD", "FREQUENCY", 50), ("BUDGET", "FREQUENCY", 50)]}
    elif pars.outputs == "last":
    	saverecord = [("HEAD", "LAST"), ("BUDGET", "LAST")]
    	
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    flopy.mf6.ModflowGwfsto(gwf, ss=pars.Ss, sy=pars.sy, iconvert=1, transient=True)

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
    
    flopy.mf6.ModflowGwtcnc(
       gwt,
    stress_period_data=cnc_data,
       pname="CNC-1",
    ) 

    flopy.mf6.ModflowGwtmst(gwt, porosity=pars.porosity)

    strt = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    strt[:, :,-1] = 35.0
    flopy.mf6.ModflowGwtic(gwt, strt=strt)

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=pars.diff, alh=pars.alpha_L, ath1=pars.alpha_L*pars.alpha_anisT, alv=pars.alpha_L*pars.alpha_anisV)

    sourcerecarray = [
        # ("GHB-1", "AUX", "CONCENTRATION"),
        ("CHD-1", "AUX", "CONCENTRATION")
    ]

    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    if pars.outputs == "all":
    	saverecord = {0: [("CONCENTRATION", "LAST")], nper_trans: [("CONCENTRATION", "LAST")], nper_trans+1: [("CONCENTRATION", "FREQUENCY", 50)]}
    elif pars.outputs == "last":
    	saverecord = [("CONCENTRATION", "LAST")]
    	
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{pars.name}.cbc",
        concentration_filerecord=f"{pars.name}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=saverecord,
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    pars.sea_levels = sea_levels
    pars.aquifer_cells = aquifer_cells
    pars.aquitard_cells = aquitard_cells
    pars.top_cells = top_cells
    pars.bottom_cells = bottom_cells
    pars.save_parameters()

    if not just_for_pars:
        return sim
    else:
        return pars

def run_model(swt):
    swt.write_simulation()
    success, buff = swt.run_simulation(silent=False)
    if not success:
            print(buff)

    return success

def get_last_results(name):
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
    
    kstpkper = gwt.output.concentration().get_kstpkper()
    kstp = [k[0] for k in kstpkper]
    predev_k = kstpkper[len(kstp)-1-kstp[::-1].index(0)]
    postdev_k = kstpkper[-1]
    predev_k=(predev_k[0], predev_k[1]-1)
    
    conc = [gwt.output.concentration().get_data(kstpkper=predev_k), gwt.output.concentration().get_data(kstpkper=postdev_k)]
    head = [gwf.output.head().get_data(kstpkper=predev_k), gwt.output.concentration().get_data(kstpkper=postdev_k)]
    bud = gwf.output.budget()
    times = bud.get_times()
    qx = np.zeros_like(conc)
    qy = np.zeros_like(conc)
    qz = np.zeros_like(conc)

    for i, time in enumerate([predev_k, postdev_k]):

        qx[i], qy[i], qz[i] = flopy.utils.postprocessing.get_specific_discharge(
            bud.get_data(text="DATA-SPDIS", kstpkper=time)[0],
            gwf,
        )
 

    return conc, head, qx, qz, times-(times[-1]-365*200)*np.ones_like(times)

def get_results(name, scatter=False, backup=False):
    
    if backup == False:
        pars = load_parameters(name)
    else:
        pars = load_parameters(name, backup==True)
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
    # .set_trace()
    bud = gwf.output.budget()
    times = bud.get_times()
    times_to_save=[times[0]] + [times[i] for i in range(-74, 0)]
    # times_to_save=times
    conc = []
    head = []
    qx = []
    qz = []

    for time in times_to_save:
        conc.append(gwt.output.concentration().get_data(totim=time))
        head.append(gwt.output.head().get_data(totim=time))
        qxi, _, qzi = flopy.utils.postprocessing.get_specific_discharge(
                bud.get_data(text="DATA-SPDIS", totim=time)[0],
                gwf,
               )
        qx.append(qxi)
        qz.append(qzi)
  

    return conc, head, qx, qz, times_to_save-(times_to_save[1]-365*200)*np.ones_like(times_to_save)
    
def get_results_for_hc(name):
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
    bud = gwf.output.budget()
    times = bud.get_times()
    # pdb.set_trace()
    conc=[]
    conc.append(gwt.output.concentration().get_data(totim=times[-2]))
    conc.append(gwt.output.concentration().get_data(totim=times[-1]))
    _, _, qz = flopy.utils.postprocessing.get_specific_discharge(
                bud.get_data(text="DATA-SPDIS", totim=times[-2])[0],
                gwf,
               )

    return conc, qz

def remove_model_files(name):
    shutil.rmtree(f"/home/superuser/sloping_confined_aquifer/model_files/{name}")


def get_time_series_spline():
    arr = np.loadtxt("time_series.txt", delimiter=" ", skiprows=1)
    arr[1, 4] = -0.01
    arr[0, 4] = 0.0
    spl_time_to_sl = scipy.interpolate.interp1d(arr[:25, 0], arr[:25, 4], fill_value='extrapolate')
    #pdb.set_trace()
    spl_sl_to_time = scipy.interpolate.interp1d(arr[:25, 4], arr[:25, 0], fill_value='extrapolate')
    #pdb.set_trace()
    return spl_time_to_sl, spl_sl_to_time 

def build_model_start_from_modern(name, n_old, n_new, new_head, results, new_dt=None):

    pars = load_parameters(f"{n_old}", backup=True)
    spl_time_to_sl, spl_sl_to_time = get_time_series_spline()
    mf6exe = pars.exe_path 
    length_units = "m"
    time_units = "days"

    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97
    name = n_new
    ws = f"./model_files/{name}"
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name=mf6exe
    )
    time_spd = []

    nper = 1
    perlen = [pars.T_modern] 
    
    time_spd.append((pars.T_modern, pars.T_modern/10, 1))
    
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

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

    delc = 1
    delv= pars.dz
    delr=pars.dx
    top = pars.Lz-pars.z0
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)
    idomain = np.zeros((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)

    aquitard_cells = []
    aquifer_cells = []
    top_cells = []
    bottom_cells = []
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    offshore_boundary_cells_end = []

    for col in range(pars.ncol):
        first = True
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz)): 
            if first:
                first=False
                if col > int(pars.x_b/pars.dx):
                    offshore_boundary_cells.append([lay, 0, col])
            aquitard_cells.append([lay, 0, col])
            idomain[lay, 0, col]=1
            if col == 0:
                onshore_boundary_cells.append([lay, 0, col])
            elif col == pars.ncol-1 and not first:
                offshore_boundary_cells_end.append([lay, 0, col])
        first = True
        for lay in range(int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+int(pars.D/pars.dz)), int(np.floor(np.tan(pars.beta)*col*pars.dx/pars.dz)+pars.D/pars.dz+pars.H/pars.dz)):           
            if first:
                first=False
                top_cells.append([lay, 0, col])
            if lay<pars.nlay:
                aquifer_cells.append([lay, 0, col])
            idomain[lay, 0, col]=1
            if col == 0:
                onshore_boundary_cells.append([lay, 0, col])
            elif col == pars.ncol-1 and not first:
                offshore_boundary_cells_end.append([lay, 0, col])
        bottom_cells.append([lay, 0, col])
    if offshore_boundary_cells_end[-1] != [pars.nlay-1, 0, pars.ncol-1]:
        aquifer_cells.append([pars.nlay-1, 0, pars.ncol-1])
        idomain[pars.nlay-1, 0, pars.ncol-1]=1
        offshore_boundary_cells_end.append([pars.nlay-1, 0, pars.ncol-1])
            
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
    hk = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    k33 = np.zeros_like(hk)

    for cell in aquifer_cells: 
        hk[cell[0], cell[1], cell[2]] = pars.K_aquifer
        k33[cell[0], cell[1], cell[2]] = pars.K_aquifer/pars.anis_aquifer

    for cell in aquitard_cells: 
        hk[cell[0], cell[1], cell[2]] = pars.K_aquitard
        k33[cell[0], cell[1], cell[2]] = pars.K_aquitard/pars.anis_aquitard

    # pdb.set_trace()
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=hk,
        k33=k33,
        rewet_record="REWET WETFCT 1.0 IWETIT 1 IHDWET 0",
        wetdry=1
    )
    strt = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    # strt[:, :, -1] = 35.7
    flopy.mf6.ModflowGwfic(gwf, strt=results["head"][1][:, :, :]) # ???

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)
    
    sea_level = 0
    #pdb.set_trace()
    ghb_data = {}
    chd_data= {}
    cnc_data= {}
    ghbcond_v = 1e5# pars.K_aquitard/pars.anis_aquitard * delr * delc / (delv)
    ghbcond_ha = 1e5 # pars.K_aquifer * delv * delc / (delr)
    ghbcond_hb = 1e5# pars.K_aquitard * delv * delc / (delr)

    ghb_temp = []
    chd_temp = []
    cnc_temp = []
    for cell in offshore_boundary_cells:  
            # if pars.Lx-cell[2]*pars.dx<=(pars.Lx-((pars.x_b*np.tan(pars.beta))-sea_levels[0])/np.tan(pars.beta))+(sl0 + spl_time_to_sl(t0-(pars.T/(365000))*(i+1)))/np.tan(pars.beta):
        if cell[2] >=(pars.z_b-sea_level)/np.tan(pars.beta)/pars.dx:
            chd_temp.append([(cell[0], cell[1], cell[2]), sea_level, 35.0])
            cnc_temp.append([(cell[0], cell[1], cell[2]), 35.0])
    
    for cell in offshore_boundary_cells_end[1:]:
        if -(sea_level) <= cell[0]*pars.dz-(pars.Lz-pars.z0):
            if cell[0] <= pars.Lz-pars.H:
                ghbcond=ghbcond_hb
            else:
                ghbcond=ghbcond_ha
            chd_temp.append([(cell[0], cell[1], cell[2]), sea_level, 35.0])
            cnc_temp.append([(cell[0], cell[1], cell[2]), 35.0])
    
    for cell in onshore_boundary_cells:
        
        h_onshore = new_head

        if cell[0]*pars.dz-pars.Lz+pars.z0 >= -h_onshore and cell[0]*pars.dz>=pars.D:
            chd_temp.append([(cell[0], cell[1], cell[2]), h_onshore, 0.0])
            
    chd_data[0] = chd_temp
    ghb_data[0] = ghb_temp
    cnc_data[0] = cnc_temp

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )
    
    
    ws_results = f"./model_files/{name}"
    if not os.path.exists(ws_results):
        os.makedirs(ws_results)

    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    if pars.outputs == "all":
        saverecord = {0: [("HEAD", "FREQUENCY", 50), ("BUDGET", "FREQUENCY", 50)]}
    	
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    flopy.mf6.ModflowGwfsto(gwf, ss=pars.Ss, sy=pars.sy, iconvert=1, transient=True)

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
    
    flopy.mf6.ModflowGwtcnc(
       gwt,
    stress_period_data=cnc_data,
       pname="CNC-1",
    ) 

    flopy.mf6.ModflowGwtmst(gwt, porosity=pars.porosity)

    strt = np.zeros((pars.nlay, pars.nrow, pars.ncol))
    strt[:, :,-1] = 35.0
    flopy.mf6.ModflowGwtic(gwt, strt=results["conc"][1][:, :, :]) # ???

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=pars.diff, alh=pars.alpha_L, ath1=pars.alpha_L*pars.alpha_anisT, alv=pars.alpha_L*pars.alpha_anisV)

    sourcerecarray = [
        # ("GHB-1", "AUX", "CONCENTRATION"),
        ("CHD-1", "AUX", "CONCENTRATION")
    ]

    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    saverecord = {0: [("CONCENTRATION", "FREQUENCY", 50)]}
    	
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=saverecord,
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    pars.sea_levels = sea_level
    pars.aquifer_cells = aquifer_cells
    pars.aquitard_cells = aquitard_cells
    pars.top_cells = top_cells
    pars.bottom_cells = bottom_cells
    pars.save_parameters()

    return sim


def get_results_modern(name, old_name, backup=False):

    if backup == False:
        pars = load_parameters(old_name)
    else:
        pars = load_parameters(old_name, backup==True)
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
    # .set_trace()
    bud = gwf.output.budget()
    times = bud.get_times()
    times_to_save = times # [times[i] for i in range(-74, 0)]
    # times_to_save=times
    conc = []
    head = []
    qx = []
    qz = []

    for time in times_to_save:
        conc.append(gwt.output.concentration().get_data(totim=time))
        head.append(gwt.output.head().get_data(totim=time))
        qxi, _, qzi = flopy.utils.postprocessing.get_specific_discharge(
                bud.get_data(text="DATA-SPDIS", totim=time)[0],
                gwf,
               )
        qx.append(qxi)
        qz.append(qzi)
  

    return conc, head, qx, qz