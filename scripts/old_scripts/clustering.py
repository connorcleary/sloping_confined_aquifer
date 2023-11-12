from results_mf6 import get_pre_concs, get_pars, get_outputs
from sklearn.cluster import MiniBatchKMeans
import matplotlib.pyplot as plt
import numpy as np 
from sklearn.linear_model import LinearRegression
import heapq

def clean_images(images):
    # remove nans, placeholders and replace with 0
    for n, image in enumerate(images):
        for i in range(image.shape[0]):
            for j in range(image .shape[1]):
                if image[i, j] == 1e30 or image[i,j] < 0:
                   image[i, j] = 0
                image[i,j] = image[i,j]/35
        images[n] = image.flatten()

    return images 

def find_clusters(images, n_clusters):
    # perform k_means clustering
    kmeans = MiniBatchKMeans(n_clusters)
    kmeans.fit(images)
    return kmeans

def plot_centroids(kmeans, shape):
    # plot centre of each cluster
    centroids = kmeans.cluster_centers_
    f, axs = plt.subplots(len(centroids), figsize=(6, 4), layout="tight")
    for i, centroid in enumerate(centroids):
        axs[i].pcolormesh(np.flipud(centroid.reshape(shape[0], shape[1])), cmap="viridis")
        axs[i].set_box_aspect(0.3)

    plt.show()
    pass

def find_greatest_effect(name, n, kmeans, n_clusters, scales, pars):
    labels = kmeans.labels_
    inputs = get_pars(name, n)
    outputs, _, _, _, _ = get_outputs(name, n, inputs)
    scores = np.zeros((n_clusters, inputs.shape[1]))
    best_pars = [] 
    # find greatest effect on each cluster
    pars.pop(4)
    for i in range(n_clusters):
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

def plot_results(name, n, kmeans, n_clusters, best_pars):
    # plot distribution for each cluster etc.
    labels = kmeans.labels_
    inputs = get_pars(name, n)
    outputs, _, _, _, _ = get_outputs(name, n, inputs)
    outputs_ = []
    for i in range(n_clusters):
        outputs_.append(np.array(outputs)[labels == i])

    f, ax = plt.subplots(figsize=(4, 4), layout="tight")
    ax.boxplot(outputs_, labels=[str(i) for i in range(n_clusters)])
    ax.set_ylim(0, 1)
    plt.show()

def find_plot_clusters(name, n, scales, pars, n_clusters=4):
    # put it all together
    images = get_pre_concs(name, n)
    shape = images[0].shape
    images = clean_images(images)
    kmeans = find_clusters(images, n_clusters)
    plot_centroids(kmeans, shape)
    best_pars = find_greatest_effect(name, n, kmeans, n_clusters, scales, pars)
    print(best_pars)
    plot_results(name, n, kmeans, n_clusters, best_pars)
    pass