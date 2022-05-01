"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
from pathlib import Path
from typing import Callable, Dict

from instance import Instance
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper

#my imports
from itertools import chain, combinations
import pulp as p 


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )


def solve(instance: Instance) -> Solution:
    #parsing possible locations
    possible_loc = set()
    for city in instance.cities:
        x = city.x
        y = city.y
        possible_loc.union(all_points_in_radius(x,y))




    possible_loc_pts = {}
    for loc in possible_loc:
        possible_loc_pts.add(Point(loc[0], loc[1]))






    #all possible tower placements
    tower_settings = list(powerset(possible_loc_pts))

    #run LP on all the possible solution sets.
    for tower_setting in tower_settings:
        #test whether we cover all cities with our towers
        if not covers_all(tower_setting, instance.cities):
            continue

        #since we can cover all the cities now we will use these settings to minimize penalty_radius
        #run LP








def covers_all(tower_setting, city_setting):
    '''checks to see if the towers in tower_setting covers all the cities in city_settings'''
    all_reachable_from_tower = set()

    for tower in tower_setting:
        all_reachable_from_tower.union(all_points_in_radius(tower.x, tower.y))

    for city in city_setting:
        if not (city in all_reachable_from_tower):
            return False

    return True


def all_points_in_radius(x,y):
    ''' finds all the points in a 3 unit radius around x,y'''
    possible_loc = set()
    for i in range(4):
        for j in range(4):
            if (sqrt(i^2 + j^2) <= 3):
                lst = [(x + i,y + j), (x - i, y - j), (x + i, y - j), (x - i, y + j)]
                for coord in lst:
                    if (coord[0] >= 0) and (coord[0]< D) and (coord[1]< D) and (coord[1] >= 0):
                        possible_loc.add(coord)

    return possible_loc







SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive,
    "actual": solve,
}

# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                        "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
