# Actions to be done on the sum of squares project

import argparse
import os

import signac

def square(*jobs):
    """ Implement the squaring operation.

    Retrieves the value stored in the 'x' parameter from a job's statepoint, squares it, and then
    writes the result to a file called 'square.out' when done.
    """

    for job in jobs:
        # If the output file already exists, don't do anything
        if job.isfile('square.out'):
            continue

        # Open a temporary file so that the action is not completed early or if an error is thrown.
        # If the .out file exists but is incorrectly generated, then we will never know!
        # this saves us!

        with open(job.fn('square.out.in_progress'), 'w') as file:
            x = job.cached_statepoint['x']
            # grab the value stored in x
            file.write(f'{x**2}')

        # Assuming no errors are thrown, we finish by renaming temp file to out file
        os.rename(job.fn('square.out.in_progress'), job.fn('square.out'))

def compute_sum(*jobs):
    """ Implementes the compute_sum action
    Prints the sum of 'square.out' from each job directory.
    """
    total = 0
    for job in jobs:
        # Read in the value stored in square.out
        with open(job.fn('square.out')) as file:
            total += int(file.read())

    print(total)

if __name__ =='__main__':
    # Parse the command line arguments: python action.py --action <ACTION> [DIRECTORIES]
    parser = argparse.ArgumentParser()
    parser.add_argument('--action', required=True)
    parser.add_argument('directories', nargs='+')
    args = parser.parse_args()

    # Open the signac jobs
    # Don't need to pass a path here I guess?
    # Probably bc this code really gets called while you are INSIDE of the project directory
    # so the current directory houses the project
    project = signac.get_project()
    jobs = [project.open_job(id=directory) for directory in args.directories]

    # Call the action
    globals()[args.action](*jobs)

