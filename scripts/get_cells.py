import numpy as np

def get_cells(pars):

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
            if col == 0:
                onshore_boundary_cells.append([lay, 0, col])
            elif col == pars.ncol-1 and not first:
                offshore_boundary_cells_end.append([lay, 0, col])
        bottom_cells.append([lay, 0, col])
    if offshore_boundary_cells_end[-1] != [pars.nlay-1, 0, pars.ncol-1]:
        aquifer_cells.append([pars.nlay-1, 0, pars.ncol-1])
        offshore_boundary_cells_end.append([pars.nlay-1, 0, pars.ncol-1])

    return aquitard_cells, aquifer_cells, top_cells, bottom_cells, onshore_boundary_cells, offshore_boundary_cells, offshore_boundary_cells_end