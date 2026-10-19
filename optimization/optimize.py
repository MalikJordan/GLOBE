import numpy as np
import time
from datetime import datetime
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from scipy.optimize import minimize
from optimization.shared import evaluate_sample


def run_optimization(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, maxiter=100, tolerance=1E-08, num_workers=4):
    
    start = time.perf_counter()
    print("Beginning Parameter optimization...")

    # Parameter names and their initial values
    parameter_names = list(independent_parameters.keys())

    # Initialize values for parameter set
    x0 = np.array([independent_parameters[name] for name in parameter_names], dtype=float)

    # Extract parameter bounds
    bounds = [[independent_bounds[name]["lower"], independent_bounds[name]["upper"]] for name in parameter_names]

    # Call objective function
    objective = partial(objective_function, base_element=base_element, model=model, reactions=reactions, physical=physical, independent_parameters=independent_parameters, constraints=constraints, objective_targets=objective_targets)

    # Run bounded gradient-based optimization
    with Pool(processes=num_workers) as pool:   result = minimize(fun=objective_function, x0=x0, method="L-BFGS-B", jac="2-point", bounds=bounds, options={"maxiter": maxiter, "ftol": float(tolerance), "workers": pool.map})

    # Convert optimized vector back into a named dictionary
    optimized_parameters = dict(zip(parameter_names, result.x))

    # Write output file
    log_optimization_results(result, parameter_names, independent_parameters, optimized_parameters, bounds, x0)

    print("Parameter optimization complete.")
    elapsed = time.perf_counter() - start
    hours, remainder = divmod(elapsed, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}\n")

    return optimized_parameters, result


def objective_function(base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets):
    return evaluate_sample(base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets)


# call_count = 0

# def run_optimization(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, maxiter=100, tolerance=1E-08):
    
#     start = time.perf_counter()
#     print("Beginning Parameter optimization...")

#     # Parameter names and their initial values
#     parameter_names = list(independent_parameters.keys())

#     # Initialize values for parameter set
#     x0 = np.array([independent_parameters[name] for name in parameter_names], dtype=float)

#     # Extract parameter bounds
#     bounds = [[independent_bounds[name]["lower"], independent_bounds[name]["upper"]] for name in parameter_names]

#     # Pass objective function into SciPy
#     def objective_function(sample):
#         global call_count
#         call_count += 1

#         # Print current evaluation
#         print(f"...evaluating sample {call_count}")
#         return evaluate_sample(base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets)

#     # Run bounded gradient-based optimization
#     result = minimize(fun=objective_function, x0=x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": maxiter, "ftol": float(tolerance)})

#     # Convert optimized vector back into a named dictionary
#     optimized_parameters = dict(zip(parameter_names, result.x))

#     # Write output file
#     log_optimization_results(result, parameter_names, independent_parameters, optimized_parameters, bounds, x0)

#     print("Parameter optimization complete.")
#     elapsed = time.perf_counter() - start
#     hours, remainder = divmod(elapsed, 3600)
#     minutes, seconds = divmod(remainder, 60)
#     print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}\n")

#     return optimized_parameters, result


def log_optimization_results(result, parameter_names, initial_parameters, optimized_parameters, bounds, initial_objective):
    """
    Definition:: Write parameter optimization results to a text log.

    Parameters
    ----------
    result : scipy.optimize.minimize
        Result returned by scipy.optimize.minimize().
    parameter_names : list[str]
        Parameter names in the same order as the optimization vector.
    initial_parameters : dict
        Initial parameter values.
    optimized_parameters : dict
        Optimized parameter values.
    bounds : list[tuple]
        Lower and upper bounds in parameter_names order.
    initial_objective : float
        Objective function evaluated at the initial parameter values.
    """

    log_path = Path.cwd() / 'parameter_optimization.log'

    with log_path.open("w", encoding="utf-8") as log:
        log.write("PARAMETER OPTIMIZATION RESULTS\n")
        log.write(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        log.write(f"Method: {result.get('method', 'L-BFGS-B')}\n")
        log.write(f"Success: {result.success}\n")
        log.write(f"Message: {result.message}\n")
        log.write(f"Iterations: {result.nit}\n")
        log.write(f"Function evaluations: {result.nfev}\n")

        log.write("\nOBJECTIVE FUNCTION\n")
        log.write("-" * 180 + "\n")
        log.write(f"Initial objective: {initial_objective:.8e}\n")
        log.write(f"Final objective:   {result.fun:.8e}\n")

        improvement = initial_objective - result.fun
        log.write(f"Absolute improvement: {improvement:.8e}\n")

        if initial_objective != 0:
            percent_improvement = (improvement / abs(initial_objective) * 100)
            log.write(f"Relative improvement: {percent_improvement:.4f}%\n")

        log.write("\nOPTIMIZED PARAMETERS\n")
        log.write("-" * 180 + "\n")
        log.write(
            f"{'Parameter':<35}"
            f"{'Initial':>24}"
            f"{'Optimized':>24}"
            f"{'Lower Bound':>24}"
            f"{'Upper Bound':>24}"
            f"{'Change (%)':>24}\n"
        )
        log.write("-" * 180 + "\n")

        for i, name in enumerate(parameter_names):
            initial = initial_parameters[name]
            optimized = optimized_parameters[name]
            lower, upper = bounds[i]

            if initial != 0:
                change_percent = ((optimized - initial) / abs(initial) * 100)
                change_text = f"{change_percent:24.5f}"
            else:
                change_text = f"{'N/A':>24}"

            log.write(
                f"{name:<35}"
                f"{initial:24.5e}"
                f"{optimized:24.5e}"
                f"{lower:24.5e}"
                f"{upper:24.5e}"
                f"{change_text}\n"
            )

    print(f"Optimization log written to: {log_path}")