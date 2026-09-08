# import pyfabm
# import os

# yaml_file = "/Users/malikjordan/fabm/testcases/fabm-jrc-med_ergom.yaml"

# print("Loading medERGOM...")
# model = pyfabm.Model(yaml_file)

# print("Model loaded.")
# print("Number of state variables:", len(model.state_variables))

# print("\nState variables:")
# for var in model.state_variables:
#     print(f"  {var.name}: {var.value}")

# # print("\nDependencies:")
# # for var in model.dependencies:
# #     print(
# #         f"  {var.name}: "
# #         f"value={var.value}, "
# #         f"units={var.units}"
# #     )

# print("\nDependency objects:")
# for dep in model.dependencies:
#     print(
#         dep.name,
#         type(dep),
#         "value =", dep.value
#     )

# print("\nInterior dependencies:")
# for dep in model.interior_dependencies:
#     print(dep.name, type(dep), dep.value)

# print("\nHorizontal dependencies:")
# for dep in model.horizontal_dependencies:
#     print(dep.name, type(dep), dep.value)

# print("\nScalar dependencies:")
# for dep in model.scalar_dependencies:
#     print(dep.name, type(dep), dep.value)

# print("\nDiagnostics:")
# for var in model.diagnostic_variables:
#     print(f"  {var.name}: {var.units}")

# print("\nSetting dependencies...")

# for dep in model.interior_dependencies:
#     print(dep.name, type(dep.value), dep.value)

# print("\nSetting dependencies...")

# for dep in model.interior_dependencies:
#     print(dep.name, type(dep.value), dep.value)

# print("\nStarting FABM...")
# model.start()

# print("FABM started.")

import pyfabm
import numpy as np

yaml_file = "/Users/malikjordan/fabm/testcases/fabm-jrc-med_ergom.yaml"

print("Loading medERGOM...")
model = pyfabm.Model(yaml_file)

print("Model loaded.")
print("Number of state variables:", len(model.state_variables))


# ------------------------------------------------------------
# Initial state
# ------------------------------------------------------------

initial_state = {
    "med_ergom/pp": 0.001,
    "med_ergom/ff": 0.001,
    "med_ergom/bb": 0.001,
    "med_ergom/zz": 0.001,
    "med_ergom/dd": 0.01,
    "med_ergom/aa": 0.001,
    "med_ergom/nn": 0.0001,
    "med_ergom/po": 0.0001,
    "med_ergom/o2": 300.0,
    "med_ergom/pw": 0.001,
    "med_ergom/fl": 100.0,
    "med_ergom/pb": 0.001,
}

for name, value in initial_state.items():
    model.findStateVariable(name).value = value


# ------------------------------------------------------------
# Set 0D geometry
# ------------------------------------------------------------

model.setCellThickness(10.0)


# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------

dependencies = {
    "med_ergom/dic": 2100.0,
    "downwelling_photosynthetic_radiative_flux": 50.0,
    "temperature": 20.0,
    "practical_salinity": 38.0,
    "surface_downwelling_photosynthetic_radiative_flux": 100.0,
    "wind_speed": 5.0,
    "bottom_stress": 0.0,
}

for name, value in dependencies.items():
    dep = model.findDependency(name)
    print(f"Setting {name} = {value}")
    dep.value = value


# ------------------------------------------------------------
# Check dependencies BEFORE start
# ------------------------------------------------------------

print("\nDependencies before start:")

for dep in model.dependencies:
    print(f"{dep.name}: {dep.value}")


# ------------------------------------------------------------
# Start FABM
# ------------------------------------------------------------

print("\nStarting FABM...")
model.start()

print("FABM started.")


# ------------------------------------------------------------
# Calculate rates
# ------------------------------------------------------------

print("\nCalculating rates...")

rates = model.getRates()

print("\nDiagnostics:")

for diagnostic in model.diagnostic_variables:
    print(
        f"{diagnostic.name}: "
        f"value={diagnostic.value}, "
        f"units={diagnostic.units}"
    )

print("Rates calculated.\n")

for state, rate in zip(model.state_variables, rates):
    print(f"{state.name}: {rate}")