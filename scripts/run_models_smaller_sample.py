from pars import ModelParameters
from mf6 import build_model_start_from_modern, run_model, get_results_modern, remove_model_files
from results_smaller import load_smaller, save_smaller
from multiprocessing import Pool


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
    heads = [0, -5]
    inputs = []

    for Ka, KbV, head in zip(Kas, KbVs, heads):
        n_old = f"10km_{Ka}_{KbV}"
        results =  load_smaller(n_old)
        n_new = f"{Ka}_{KbV}_{head}"
        inputs.append([n_old, n_new, head, results])

    p=Pool(processes=)
    p.map(modern_model_run, inputs)

def main():
    run_modern_models()

if __name__=="__main__":
    main()