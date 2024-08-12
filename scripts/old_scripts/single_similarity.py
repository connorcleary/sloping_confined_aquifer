from results_mf6 import get_pre_concs, get_pars, get_outputs
from sklearn.cluster import MiniBatchKMeans
import matplotlib.pyplot as plt
import numpy as np 
from sklearn.linear_model import LinearRegression
import heapq
from pars import load_parameters

def clean_images(images):
    # remove nans, placeholders and replace with 0
    for n, image in enumerate(images):
        for i in range(image.shape[0]):
            for j in range(image .shape[1]):
                if image[i, j] == 1e30 or image[i,j] < 0:
                   image[i, j] = -35*1e-6
                image[i,j] = image[i,j]/35
        images[n] = image.flatten()

    return images 

def find_neighbors(images, index, n):
    # find n closest neighbors to image at index
    distances = []
    for image in images:
        distances.append(np.linalg.norm(image - images[index]))
    indices_best = np.array(distances).argsort()[1:n+1]
    return indices_best

def mask_inactive_cells(image):
    for i in range(len(image)):
        if image[i]== -1e-6:
            image[i] = np.nan
    return image

def plot_spread_and_furthest_neighbor(name, images, index, indices, shape):
    pars = load_parameters(f"{name}{0}", backup=True)
    x = np.linspace(-pars.x_b, pars.L, pars.ncol)/1000
    y = np.linspace(-pars.z0, pars.Lz-pars.z0, pars.nlay)
    f, axs = plt.subplots(3, 1, figsize=(4, 7), gridspec_kw={'height_ratios': [1, 1, 3]}, layout="tight")
    plt.rcParams.update({'font.size': 9})
    axs[0].pcolormesh(x, y, np.flipud(mask_inactive_cells(images[indices[-1]]).reshape(shape[0], shape[1])), cmap="viridis")
    axs[1].pcolormesh(x, y, np.flipud(mask_inactive_cells(images[index]).reshape(shape[0], shape[1])), cmap="viridis")
    axs[0].get_shared_x_axes().join(axs[0], axs[1])
    axs[0].set_box_aspect(0.3)
    axs[1].set_box_aspect(0.3)
    axs[1].set_xlabel("Distance offshore [km]")
    axs[0].set_ylabel("Depth [masl]")
    axs[1].set_ylabel("Depth [masl]")
    axs[0].set_title("Selected model")
    axs[1].set_title("Furthest neighbor")
    axs[2].set_ylabel("Fraction OFG Salinized")
    n = len(images)
    inputs = get_pars(name, n)
    outputs, _, _, _, _ = get_outputs(name, n, inputs)
    outputs_ = []
    outputs_.append(outputs)
    outputs_.append(np.array(outputs)[indices])
    axs[2].boxplot(outputs_, labels=["Total", "Neighbors"])
    axs[2].set_ylim(0, 1)
    axs[2].set_box_aspect(1)
    plt.show()

def find_greatest_effect(name, n, best, scales, pars):

    labels = np.array([int(i in best) for i in range(n)])

    inputs = get_pars(name, n)
    outputs, _, _, _, _ = get_outputs(name, n, inputs)
    scores = np.zeros((2, inputs.shape[1]))
    best_pars = [] 
    # find greatest effect on each cluster
    pars.pop(4)
    for i in range(2):
        inputs_i = inputs[labels == i]
        outputs_i = inputs[labels == i]
        for j in range(inputs.shape[1]):
            if scales[j] == "log":
                inputs_ = np.log10(inputs[labels == i, j])
            else:
                inputs_ = inputs[labels == i, j]

            reg = LinearRegression().fit(inputs_.reshape(-1, 1), np.log10(np.array(outputs)[labels == i]))
            scores[i, j] = reg.score(inputs_.reshape(-1, 1), np.log10(np.array(outputs)[labels == i]))
        best_pars.append(np.array(pars)[scores[i, [0, 1, 2, 3, 5, 6]].argsort()[-2:]])
    return best_pars

def all_single_similarity(name, n, scales, pars, index, n_neighbors):
    images = get_pre_concs(name, n)
    shape = images[0].shape
    images = clean_images(images)
    best = find_neighbors(images, index, n_neighbors)
    plot_spread_and_furthest_neighbor(name, images, index, best, shape)
    best_pars = find_greatest_effect(name, n, best, scales, pars)
   

