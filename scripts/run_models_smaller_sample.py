from pars import ModelParameters
from mf6 import build_steady_model, build_model_start_from_modern, run_model, get_last_results, get_results_modern, remove_model_files
from results_smaller import load_smaller, save_smaller
from multiprocessing import Pool

def model_run(inputs):
    Ka, KbV = inputs
    pars = ModelParameters(f"10km_{Ka}_{KbV}", K_aquifer=Ka, K_aquitard=KbV*100)
    sim = build_steady_model(pars)
    run_model(sim)
    conc, head, qx, qz = get_last_results(f"10km_{Ka}_{KbV}")
    save_smaller(f"10km_{Ka}_{KbV}", [conc, head, qx, qz])
    pass

def run_base_models():
    Kas = [10, 10]
    KbVs = [1e-3, 1e-3]
    inputs = []

    for Ka, KbV, in zip(Kas, KbVs):
        inputs.append([Ka, KbV])

    p=Pool(processes=8)
    p.map(model_run, inputs)


def modern_model_run(inputs):
    n_old, n_new, new_head, results = inputs
    sim = build_model_start_from_modern("base10", n_old, n_new, new_head, results, new_dt=None)
    run_model(sim)
    conc, head, qx, qz = get_results_modern(n_new, n_old)
    save_smaller(n_new, [conc, head, qx, qz])
    remove_model_files(n_new)
    pass

def run_modern_models():

    Kas = [10, 10]
    KbVs = [1e-3, 1e-3]
    heads = [1, 0, -1, -2, -5]
    inputs = []

    for Ka, KbV, head in zip(Kas, KbVs, heads):
        n_old = f"10km_{Ka}_{KbV}"
        results =  load_smaller(n_old)
        n_new = f"{Ka}_{KbV}_{head}"
        inputs.append([n_old, n_new, head, results])

    p=Pool(processes=8)
    p.map(modern_model_run, inputs)

def main():
    run_modern_models()

if __name__=="__main__":
    main()