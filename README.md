# schedEnergyCom

This repository contains the data and the formulations related to the paper: 

Sangaré, M, Bourreau, E., Fortz, B., Pachurka, A., and Poss, M. Loads scheduling for demand response in energy communities. Computers & OR, In press. [pdf](https://hal.science/hal-03880548)]

Notice:
The file named <code>pasEchange.jl</code> calculates the solution when the members do not exchange their surplus within the community. These solutions are used as initial solutions for the column generation-based heuristics and to ensure that members' situations do not deteriorate after joining a community.

Before solving the <code>MILP with SOS1 variables</code>, please set <code>CPX_PARAM_WORKDIR</code> in the file <code>MILPSOS.jl</code>.

While solving the problem with the column generation-based heuristic, if you set <code>act=1</code>, the subproblem returns the first feasible solution found by the solver.
