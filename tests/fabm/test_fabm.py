import pyfabm

# --------------------------------------------------
# Load FABM NPZD
# --------------------------------------------------

yaml_file = "/Users/malikjordan/fabm/testcases/fabm-gotm-npzd.yaml"

model = pyfabm.Model(yaml_file)

print("Model loaded.")
print("Number of state variables:", len(model.state_variables))


# --------------------------------------------------
# State variables
# --------------------------------------------------

print("\nState variables:")

for state in model.state_variables:
    print(f"{state.name}: {state.value}")


# --------------------------------------------------
# Set cell thickness
# --------------------------------------------------

model.cell_thickness = 10.0

# print("\nCell thickness:", model.cell_thickness)


# --------------------------------------------------
# Dependencies
# --------------------------------------------------

print("\nDependencies:")

for dependency in model.dependencies:
    print(
        f"{dependency.name}: "
        f"value={dependency.value}, "
        f"units={dependency.units}"
    )


# --------------------------------------------------
# Set dependencies
# --------------------------------------------------

model.dependencies[
    "surface_downwelling_photosynthetic_radiative_flux"
].value = 100.0

model.dependencies[
    "downwelling_photosynthetic_radiative_flux"
].value = 100.0

model.dependencies[
    "npzd/dic"
].value = 2100.0


# --------------------------------------------------
# Start FABM
# --------------------------------------------------

print("\nStarting FABM...")

model.start()

print("FABM started.")


# --------------------------------------------------
# Calculate rates
# --------------------------------------------------

print("\nCalculating rates...")

rates = model.getRates()

print("Rates calculated.")


# --------------------------------------------------
# Print rates
# --------------------------------------------------

print("\nNPZD rates:")

for state, rate in zip(model.state_variables, rates):
    print(f"{state.name}: {rate}")

print("\nFABM NPZD test completed successfully.")