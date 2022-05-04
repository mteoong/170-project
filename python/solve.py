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
from pulp import *
from math import *
from point import Point


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )


def solve_greedy(instance: Instance) -> Solution:
    #parsing possible locations
    all_towers = []
    possible_loc = {}
    list_of_cities = set()

    for city in instance.cities:
        x = city.x
        y = city.y
        construct_dict(x,y, instance.coverage_radius, instance, possible_loc, city)
        list_of_cities.add(city)

    while True:
        tower_potentials = max_points(possible_loc)
        if not tower_potentials:
            break
        while tower_potentials:
            max_penalty = -1
            max_tower = None
            for tower in tower_potentials:
                pt_tower = Point(tower[0],tower[1])
                all_towers.append(pt_tower)
                sol = Solution(all_towers, instance)
                penalty = sol.penalty()
                if penalty > max_penalty:
                    max_penalty  = penalty
                    max_tower = pt_tower 

                all_towers.remove(pt_tower)
            all_towers.append(max_tower)
            print(max_tower.x, max_tower.y)
            print("line 1")
            decrement(possible_loc, tower_potentials, (max_tower.x, max_tower.y), instance, 3)
            print("line 2")

    return Solution(
        instance=instance,
        towers=all_towers,
    )
                
                
        #we run through all the first tower potentials and calculate penalty
        #choose the tower with least penalty if two have same random
        #iterate through all the points around the chosen tower and decrement by one
        #remove decremented towers from tower potentials
        #if tower potentials is empty
        #rerun
        #else
        #choose one of the tower potentials left
           

def decrement(possible_loc, tower_potentials, winning_tower,instance, rad):
    '''for each city covered by WINNING_TOWER remove 1 from all towers in WINNING_TOWER Radius'''
    for city in possible_loc[winning_tower]:
        x = winning_tower[0]
        y = winning_tower[1]
        
        for i in range(rad + 1):
            for j in range(rad + 1):
                if (sqrt(i**2 + j**2) <= rad):
                    lst = [(x + i,y + j), (x - i, y - j), (x + i, y - j), (x - i, y + j), (x,y)]
                    for coord in lst:
                        if (coord[0] >= 0) and (coord[0]< instance.D) and (coord[1]< instance.D) and (coord[1] >= 0):
                            coord_to_pt = (coord[0], coord[1])
                            possible_loc[coord_to_pt].remove(city)
                            if coord_to_pt in tower_potentials:
                                tower_potentials.remove(coord_to_pt)

   

def max_points(possible_loc):
    maxs = 1
    set_of_keys = set()
    for point in possible_loc:
        if len(possible_loc[point]) == maxs:
            set_of_keys.add(point)

        elif len(possible_loc[point]) > maxs:
            set_of_keys = set()
            set_of_keys.add(point)
            maxs = len(possible_loc[point])
            
    return set_of_keys
        
def construct_dict(x, y, rad, instance, possible_loc, city):
    ''' constructs the dictionary'''

    for i in range(rad + 1):
        for j in range(rad + 1):
            if (sqrt(i**2 + j**2) <= rad):
                lst = [(x + i,y + j), (x - i, y - j), (x + i, y - j), (x - i, y + j), (x,y)]
                for coord in lst:
                    if (coord[0] >= 0) and (coord[0]< instance.D) and (coord[1]< instance.D) and (coord[1] >= 0):
                        coord_to_pt = (coord[0], coord[1])
                        if coord_to_pt in possible_loc:
                            possible_loc[coord_to_pt].add(city)
                        else:
                            possible_loc[coord_to_pt] = set()
                            possible_loc[coord_to_pt].add(city)

    return possible_loc



SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive,
    "greedy": solve_greedy,
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




'''
    def solve_not_naive(instance: Instance) -> Solution:
        #parsing possible locations
        possible_loc = set()
        for city in instance.cities:
            x = city.x
            y = city.y
            possible_loc = possible_loc.union(all_points_in_radius(x,y, instance.coverage_radius, instance))


        possible_loc_pts = set()
        for loc in possible_loc:
            #print(loc)
            possible_loc_pts.add(Point(loc[0], loc[1]))



        #all possible tower placements -> list of tuples
        #print(possible_loc_pts)


        LP_Prob = LpProblem("problem")
        LP_tower_var = {}
        coverage_var = coverage_functions(possible_loc_pts, instance.cities, instance.coverage_radius, instance)
        penalty_var = coverage_functions(possible_loc_pts, possible_loc_pts, instance.penalty_radius, instance)
        print(LP_tower_var, coverage_var, penalty_var)

        #Tower variables
        for tower in possible_loc_pts:
            tower_name = 'x' + str(tower.x) + ',y' + str(tower.y)
            LP_tower_var[tower_name] = LpVariable(tower_name, cat = 'Binary')

        #objective
        LP_Prob += lpSum(170 * exp(.17* lpSum(penalty_var[tower][tower_two] *  LP_tower_var['x' + str(tower_two.x) + ',y' + str(tower_two.y)]
        for tower_two in possible_loc_pts) - 1) * LP_tower_var['x' + str(tower.x) + ',y' + str(tower.y)] for tower in possible_loc_pts)

        #constraints
        LP_Prob += lpSum(lpSum(coverage_var[city][tower] * LP_tower_var['x' + str(tower.x) + ',y' + str(tower.y)]
        for tower in possible_loc_pts) >= 1 for city in instance.cities)

        LP_Prob.solve()
        print("Status:", LpStatus[LP_Prob.status])

        for s in LP_Prob.variables():
            print(s.name, '=', s.varValue)

        print("obj value = ", value(LP_Prob.objective))
'''

'''

    lst = set()

    for tower in tower_setting:
        all_reachable_from_tower.union(all_points_in_radius(tower.x, tower.y, radius)

    for city in instance.cities:
        locs = all_points_in_radius(city.x, city.y, radius)
        for tower in tower_setting:
            if tower in locs:
                lst.add((tower, 1))

            else:
                lst.add((tower,0))


    return lst'''

'''''''''     
def coverage_functions(tower_setting, city_setting, radius, instance):
    checks to see if the towers in tower_setting covers all the cities in city_settings
    #set of tuples that holds aij = 1 if tower j covers city i
    tower_dict = {}

    for city in city_setting:
        locs = all_points_in_radius(city.x, city.y, radius, instance)
        tower_dict[city] = {}
        for tower in tower_setting:
            if tower in locs:
                tower_dict[city][tower] = 1

            else:
                tower_dict[city][tower] = 0


    return tower_dict
'''
