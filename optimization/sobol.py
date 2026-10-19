import numpy as np
import time
from datetime import datetime
from pathlib import Path
from SALib.sample.sobol import sample
from SALib.analyze.sobol import analyze
from optimization.shared import evaluate_sample, parallel_evaluation


def run_sobol_analysis(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, N=2, threshold=0.01, confidence_factor=2, num_workers=4, seed=1234):
    
    start = time.perf_counter()
    print("Beginning Sobol' Sensitivity Analysis...")

    # Parameter names and their initial values
    parameter_names = list(independent_parameters.keys())

    # Construct Sobol' problem definition
    problem = {
        "num_vars": len(parameter_names),
        "names": parameter_names,
        "bounds": [[independent_bounds[name]["lower"], independent_bounds[name]["upper"]] for name in parameter_names]
    }

    # Verify explorable range
    for name,(lower,upper) in zip(problem["names"],problem["bounds"]):
        if not np.isfinite(lower) or not np.isfinite(upper):    raise ValueError(f"Non-finite bounds for '{name}'.")
        if lower >= upper:  raise ValueError(f"Parameter '{name}' must have lower < upper; got [{lower},{upper}].")

    # Generate Sobol' samples
    X = sample(problem, N=N, calc_second_order=False, seed=seed)

    print(f"Independent parameters: {len(parameter_names)}")
    print(f"Base sample size: {N}")
    print(f"Model evaluations: {len(X)}")

    # Evaluate samples in parallel
    Y = parallel_evaluation(X, base_element, model, reactions, physical, independent_parameters, constraints, objective_targets, num_workers=4)

    # Calculate sensitivity indices
    results = analyze(problem, Y, calc_second_order=False, num_resamples=100, seed=seed)

    # Screen parameters
    retained, removed, relative_importance = screen_parameters(results, parameter_names, threshold, confidence_factor)

    # Write output file
    log_sobol_results(results, parameter_names, retained, removed, relative_importance)

    print("Sobol' Sensitivity Analysis complete.")
    elapsed = time.perf_counter() - start
    hours, remainder = divmod(elapsed, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}\n")

    return X, Y, results, retained, removed, relative_importance


# def run_sobol_analysis(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, N=2, threshold=0.01, confidence_factor=2, seed=1234):

#     start = time.perf_counter()
#     print("Beginning Sobol' Sensitivity Analysis...")

#     # Parameter names and their initial values
#     parameter_names = list(independent_parameters.keys())

#     # Construct Sobol' problem definition
#     problem = {
#         "num_vars": len(parameter_names),
#         "names": parameter_names,
#         "bounds": [[independent_bounds[name]["lower"], independent_bounds[name]["upper"]] for name in parameter_names]
#     }

#     # Verify explorable range
#     for name,(lower,upper) in zip(problem["names"],problem["bounds"]):
#         if not np.isfinite(lower) or not np.isfinite(upper):    raise ValueError(f"Non-finite bounds for '{name}'.")
#         if lower >= upper:  raise ValueError(f"Parameter '{name}' must have lower < upper; got [{lower},{upper}].")

#     # Generate Sobol' samples
#     X = sample(problem, N=N, calc_second_order=False, seed=seed)

#     print(f"Independent parameters: {len(parameter_names)}")
#     print(f"Base sample size: {N}")
#     print(f"Model evaluations: {len(X)}")

#     # Test to skip GLOBE simulations
#     rng = np.random.default_rng(seed)
#     Y = rng.random(X.shape[0])
    
#     # # Run GLOBE for every sample
#     # Y = np.empty(X.shape[0], dtype=float)

#     # for i, sample_i in enumerate(X):
#     #     print(f"...evaluating sample {i+1} of {len(X)}")
#     #     Y[i] = evaluate_sample(base_element, model, reactions, physical, independent_parameters, constraints, sample_i, objective_targets)

#     #     if not np.isfinite(Y[i]):   raise ValueError(f"Non-finite objective for sample {i}: {Y[i]}")

#     # Calculate sensitivity indices
#     results = analyze(problem, Y, calc_second_order=False, num_resamples=100, seed=seed)

#     # Screen parameters
#     retained, removed, relative_importance = screen_parameters(results, parameter_names, threshold, confidence_factor)

#     # Write output file
#     log_sobol_results(results, parameter_names, retained, removed, relative_importance)

#     print("Sobol' Sensitivity Analysis complete.")
#     elapsed = time.perf_counter() - start
#     hours, remainder = divmod(elapsed, 3600)
#     minutes, seconds = divmod(remainder, 60)
#     print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}\n")

#     return X, Y, results, retained, removed, relative_importance


def screen_parameters(results, parameter_names, threshold=0.01, confidence_factor=1):
    """
    Definition:: Determines which parameters to remove from parameter space (post-Sobol) based on relative importance.
                 Relative importance is calculated as {ST + confidence_factor * ST_conf}
                 Parameters are retained if relative importance > 0.01, i.e. greater than 1% of output variance.

    ----------------------------------------------------------------------
    Parameters
    ----------
    results : dict
        Dictionary containing parameter names, mean elementary effect (mu), mean absolute elementary effect (mu_star), and standard deviation of elementary effect (sigma)
    parameter_names : list
        List containing names of parameters that were evaluated in Morris sampling and analysis
    threshold : float
        Threshold value used for detemining which parameters are necessary to retain for optimization. Parameter retained if combined sensitivity score >= threshold.
    confidence_factor: float
        Multiplier applied to ST confidence interval width
    ----------------------------------------------------------------------
    Returns
    ----------
    retained : list
        List of retained parameters
    removed : list
        List of removed parameters
    relative_importance : dict
        Dictionary containing relative importance for each sampled parameter
    """
    st = np.asarray(results["ST"], dtype=float)
    st_conf = np.asarray(results["ST_conf"], dtype=float)

    # rng = np.random.default_rng(seed=1234)
    # st = rng.uniform(0.005, 0.02, size=len(parameter_names))
    # st_conf = 1E-05 * np.ones(len(parameter_names))
    
    # Initialize lists
    removed = []
    retained = []
    importance = np.zeros(len(st), dtype=float)

    for i, name in enumerate(parameter_names):
        # Calculate upper limit of uncertainty
        importance[i] = st[i] + confidence_factor * st_conf[i]

        # Determine if parameter is candidate for removal
        remove = (np.isfinite(st[i]) and np.isfinite(st_conf[i]) and importance[i] < threshold)

        if remove:  removed.append(name)
        else:       retained.append(name)

    # Create dictionary for relative impotance
    relative_importance = dict(zip(parameter_names, importance))

    return retained, removed, relative_importance


def log_sobol_results(results, parameter_names, retained, removed, relative_importance):
    """
    Definition:: Writes final Morris results and retained/removed parameters to a log

    Parameters
    ----------
    result : SALib.analyze.sobol.analyze
        Result returned by SALib.analyze.sobol.analyze().
    parameter_names : list[str]
        Parameter names in the same order as the optimization vector
    retained : list
        Retained parameters after screening
    removed : list
        Removed parameters after screening
    relative_importance : dict
        Relative importance values that were used to retain/remove parameters
    """

    # log_path = os.getcwd() + "/sobol_output.log"
    log_path = Path.cwd() / 'sobol_output.log'
    # output_dir = Path(output_dir)
    # output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # log_path = output_dir / f"morris_results_{timestamp}.log"

    s1 = results["S1"]
    st = results["ST"]
    st_conf = results["ST_conf"]

    # order = sorted(relative_importance.items(), key=lambda item: item[1], reverse=True)

    # Sort parameters by relative importance from highest to lowest
    relative_importance = dict(sorted(relative_importance.items(), key=lambda item: item[1], reverse=True))

    # Get parameter names and indices in descending order
    sorted_parameters = list(relative_importance.keys())
    order = [parameter_names.index(parameter) for parameter in sorted_parameters]

    # # Sort results in order of relative importance
    # s1 = results["S1"][sorted_indices]
    # st = results["ST"][sorted_indices]
    # st_conf = results["ST_conf"][sorted_indices]

    with log_path.open("w", encoding="utf-8") as log:
        log.write("SOBOL' SENSITIVITY ANALYSIS RESULTS\n")
        log.write(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        log.write(f"Parameters screened: {len(parameter_names)}\n")
        log.write(f"Parameters retained: {len(retained)}\n")
        log.write(f"Parameters removed: {len(removed)}\n\n")

        log.write("SENSITIVITY RESULTS\n")
        log.write("-" * 180 + "\n")
        log.write(
            f"{'Parameter':<35}"
            f"{'S1':>24}"
            f"{'ST':>24}"
            f"{'ST_conf':>24}"
            f"{'Rel. Importance':>24}\n"
        )
        log.write("-" * 180 + "\n")

        for i in order:
            log.write(
                f"{parameter_names[i]:<35}"
                f"{s1[i]:>24.5e}"
                f"{st[i]:>24.5e}"
                f"{st_conf[i]:>24.5e}"
                f"{relative_importance[parameter_names[i]]:>24.5e}\n"
            )

        log.write("\nRETAINED PARAMETERS\n")
        log.write("-" * 80 + "\n")
        for name in sorted(retained):
            log.write(f"{name}\n")

        log.write("\nREMOVED PARAMETERS\n")
        log.write("-" * 80 + "\n")
        for name in sorted(removed):
            log.write(f"{name}\n")

    print(f"Sobol' Sensitivity Analysis results saved to: {log_path}")
