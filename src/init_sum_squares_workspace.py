# Populate the workspace for the sum of squares example Row Tutorials

import signac

N = 10

project = signac.get_project("project_sum_squares")

for x in range(N):
    job = project.open_job({'x': x}).init()

    