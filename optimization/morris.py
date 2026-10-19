import numpy as np
import time
from datetime import datetime
from pathlib import Path
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze
from optimization.shared import evaluate_sample, parallel_evaluation


def run_morris_screening(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, num_trajectories=10, num_levels=4, num_workers=4, threshold=0.5, seed=1234):
    
    start = time.perf_counter()
    print("Beginning Morris Screening...")
    
    # Parameter names and their initial values
    parameter_names = list(independent_parameters.keys())

    # Construct Morris problem definition
    problem = {
        "num_vars": len(parameter_names),
        "names": parameter_names,
        "bounds": [[independent_bounds[name]["lower"], independent_bounds[name]["upper"]] for name in parameter_names]
    }

    # Verify explorable range
    for name,(lower,upper) in zip(problem["names"],problem["bounds"]):
        if not np.isfinite(lower) or not np.isfinite(upper):    raise ValueError(f"Non-finite bounds for '{name}'.")
        if lower >= upper:  raise ValueError(f"Parameter '{name}' must have lower < upper; got [{lower},{upper}].")

    # Generate Morris samples
    X = sample(problem, N=num_trajectories, num_levels=num_levels, seed=seed)

    print(f"Independent parameters: {len(parameter_names)}")
    print(f"Morris trajectories: {num_trajectories}")
    print(f"Model evaluations: {len(X)}")

    # Evaluate samples in parallel
    Y = parallel_evaluation(X, base_element, model, reactions, physical, independent_parameters, constraints, objective_targets, num_workers)

    # Calculate Morris sensitivity indices
    results = analyze(problem, X, Y, num_levels=num_levels, num_resamples=100, seed=seed)

    # Screen parameters
    retained, removed, sensitivity_scores = screen_parameters(results, parameter_names, threshold)

    # Write output file
    log_morris_results(results, parameter_names, retained, removed)

    print("Morris Screening complete.")
    elapsed = time.perf_counter() - start
    hours, remainder = divmod(elapsed, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}\n")

    return X, Y, results, retained, removed, sensitivity_scores


# def run_morris_screening(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, num_trajectories=10, num_levels=4, threshold=0.5, seed=1234):

#     start = time.perf_counter()
#     print("Beginning Morris Screening...")
    
#     # Parameter names and their initial values
#     parameter_names = list(independent_parameters.keys())

#     # Construct Morris problem definition
#     problem = {
#         "num_vars": len(parameter_names),
#         "names": parameter_names,
#         "bounds": [[independent_bounds[name]["lower"], independent_bounds[name]["upper"]] for name in parameter_names]
#     }

#     # Verify explorable range
#     for name,(lower,upper) in zip(problem["names"],problem["bounds"]):
#         if not np.isfinite(lower) or not np.isfinite(upper):    raise ValueError(f"Non-finite bounds for '{name}'.")
#         if lower >= upper:  raise ValueError(f"Parameter '{name}' must have lower < upper; got [{lower},{upper}].")

#     # Generate Morris samples
#     X = sample(problem, N=num_trajectories, num_levels=num_levels, seed=seed)

#     print(f"Independent parameters: {len(parameter_names)}")
#     print(f"Morris trajectories: {num_trajectories}")
#     print(f"Model evaluations: {len(X)}")

#     # Test to skip GLOBE simulations
#     rng = np.random.default_rng(seed)
#     Y = rng.random(X.shape[0])

#     # # Run GLOBE for every sample
#     # Y = np.empty(X.shape[0], dtype=float)

#     # for i, sample_i in enumerate(X):
#     #     print(f"... evaluating sample {i+1}/ of {len(X)}")
#     #     Y[i] = evaluate_sample(base_element, model, reactions, physical, independent_parameters, constraints, sample_i, objective_targets)

#     #     if not np.isfinite(Y[i]):   raise ValueError(f"Non-finite objective for sample {i}: {Y[i]}")

#     #     # if (i + 1) % 10 == 0 or i == len(X) - 1:    print(f"Completed {i + 1}/{len(X)} evaluations")

#     # Calculate Morris sensitivity indices
#     results = analyze(problem, X, Y, num_levels=num_levels, num_resamples=100, seed=seed)

#     # Screen parameters
#     retained, removed, sensitivity_scores = screen_parameters(results, parameter_names, threshold)

#     # Write output file
#     log_morris_results(results, parameter_names, retained, removed)

#     print("Morris Screening complete.")
#     elapsed = time.perf_counter() - start
#     hours, remainder = divmod(elapsed, 3600)
#     minutes, seconds = divmod(remainder, 60)
#     print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}\n")

#     return X, Y, results, retained, removed, sensitivity_scores


def screen_parameters(results, parameter_names, threshold):
    """
    Definition:: Determines which parameters to remove from parameter space (pre-Sobol) based on combined sensitivity score.
                 Sensitivity score calculated as the Euclidian distance defined by Ciric et al. (2012)

    ----------------------------------------------------------------------
    Parameters
    ----------
    results : dict
        Dictionary containing parameter names, mean elementary effect (mu), mean absolute elementary effect (mu_star), and standard deviation of elementary effect (sigma)
    parameter_names : list
        List containing names of parameters that were evaluated in Morris sampling and analysis
    threshold : float
        Threshold value used for detemining which parameters are necessary for further analysis in Sobol' SA. Parameter retained if combined sensitivity score >= threshold.

    ----------------------------------------------------------------------
    Returns
    ----------
    retained : list
        List of retained parameters
    removed : list
        List of removed parameters
    sensitivity_scored : dict
        Dictionary containing combined sensitivity score for each sampled parameter
    """
    mu_star = np.asarray(results["mu_star"], dtype=float)
    sigma = np.asarray(results["sigma"], dtype=float)

    # Combined Morris sensitivity score (RMS of mean absolute elementary effect and standard deviation of elementary effect)
    scores = np.sqrt(mu_star**2 + sigma**2)

    # Normalize scores by the maximum sensitivity score to [0, 1]
    max_score = np.max(scores)

    if max_score > 0:   normalized_scores = scores / max_score
    else:   normalized_scores = np.zeros_like(scores)

    # Create dictionary for sensitivity scores
    sensitivity_scores = dict(zip(parameter_names, normalized_scores))

    # Create list of retained parameters and their sensitivity score
    retained = [name for name, score in sensitivity_scores.items() if score >= threshold]

    # Create list of removed parameters and their sensitivity score
    removed = [name for name, score in sensitivity_scores.items() if score < threshold]

    return retained, removed, sensitivity_scores


def log_morris_results(results, parameter_names, retained, removed):
    """
    Definition:: Writes final Morris results and retained/removed parameters to a log

    Parameters
    ----------
    result : SALib.analyze.morris.analyze
        Result returned by SALib.analyze.morris.analyze().
    parameter_names : list[str]
        Parameter names in the same order as the optimization vector
    retained : list
        Retained parameters after screening
    removed : list
        Removed parameters after screening
    """

    # log_path = os.getcwd() + "/morris_output.log"
    log_path = Path.cwd() / 'morris_output.log'
    # output_dir = Path(output_dir)
    # output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # log_path = output_dir / f"morris_results_{timestamp}.log"

    mu = results["mu"]
    mu_star = results["mu_star"]
    sigma = results["sigma"]
    mu_star_conf = results["mu_star_conf"]

    score = np.sqrt(mu_star**2 + sigma**2)
    normalized_score = (score / np.max(score)if np.max(score) > 0 else np.zeros_like(score))

    order = np.argsort(normalized_score)[::-1]

    with log_path.open("w", encoding="utf-8") as log:
        log.write("MORRIS SCREENING RESULTS\n")
        log.write(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        log.write(f"Parameters screened: {len(parameter_names)}\n")
        log.write(f"Parameters retained: {len(retained)}\n")
        log.write(f"Parameters removed: {len(removed)}\n\n")

        log.write("SENSITIVITY RESULTS\n")
        log.write("-" * 180 + "\n")
        log.write(
            f"{'Parameter':<35}"
            f"{'mu':>24}"
            f"{'mu_star':>24}"
            f"{'sigma':>24}"
            f"{'mu_star_conf':>24}"
            f"{'Score':>24}"
            f"{'Norm. Score':>24}\n"
        )
        log.write("-" * 180 + "\n")

        for i in order:
            log.write(
                f"{parameter_names[i]:<35}"
                f"{mu[i]:>24.5e}"
                f"{mu_star[i]:>24.5e}"
                f"{sigma[i]:>24.5e}"
                f"{mu_star_conf[i]:>24.5e}"
                f"{score[i]:>24.5e}"
                f"{normalized_score[i]:>24.5f}\n"
            )

        log.write("\nRETAINED PARAMETERS\n")
        log.write("-" * 80 + "\n")
        for name in sorted(retained):
            log.write(f"{name}\n")

        log.write("\nREMOVED PARAMETERS\n")
        log.write("-" * 80 + "\n")
        for name in sorted(removed):
            log.write(f"{name}\n")

    print(f"Morris results saved to: {log_path}")
