# for the initialization of the test Signac project
import signac

# Initialize the project in the current working directory
project = signac.init_project("ideal_gas_project")

# Iterate over multiple pressures to define the state space
for p in range(1,10):
    sp = {"p": p, "kT": 1.0, "N": 1000}
    job = project.open_job(sp).init()