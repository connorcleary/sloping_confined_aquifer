import numpy as np

def save_smaller(name, results):
    np.save(f"/home/superuser/sloping_confined_aquifer/results/base10/{name}_results.npy", results)

def load_smaller(name):
    results = {}
    results_raw = np.load(f"/home/superuser/sloping_confined_aquifer/results/base10/{name}_results.npy")  
    results["conc"] = results_raw[0]
    results["head"] = results_raw[1]
    results["qx"] = results_raw[2]
    results["qz"] = results_raw[3]

    return results

