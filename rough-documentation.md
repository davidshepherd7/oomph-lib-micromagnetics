
Some notes on my code:


Potential gotchas
-------------------------

Some gotchas in my micromagnetics code aren't/can't be caught by paranoid
checks (some of these should be fixed eventually, some of these cannot
really be fixed):

* Bem by default assumes that it can automatically locate sharp edges/corners
  by looking for nodes with multiple boundaries. It also assumes that these
  edges/corners are right angles!

* Midpoint method by bdf takes a half time step, this is done by halving the
  step size and reducing time by dt/2 before implicit steps. Then an update
  is done afterwards to more to the final time. You should be careful with
  the order of operations to ensure that e.g. predictions happen before the
  time step is halved, boundary condition updates are correct, post step
  updates are consistent.... etc.

* Midpoint method by residual/Jacobian modification requires complex
  modifications to the residual and Jacobian calculations. Be sure that you

  1. Include the d_u_evaltime_by_du_np1 term in the Jacobian
  2. Properly time-interpolate time/values/x/derivatives in jacobian/residual
   (all of these MUST be evaluated at the midpoint)
  3. Make sure that boundary conditions are consistent


* Possibly others, as indicated by compile warnings...


* Make doesn't actually include all dependencies properly yet... Need to
  set it up to automatically rebuild meshes and build git version header!


Use of factory methods
--------------------

In constrast to most of oomph-lib I've tried to minimise duplication of
problem specification code. I've also tried to allow everything to be
controlled via the command line. To allow this I've used a lot of factory
methods. The idea is to have a function which creates objects for you, so
given some parameter the function will create and return a pointer to an
object.

Advantages:
* Far less time spent compliling
* Debug options are always there ready to use: saves a lot of time
* Generic options like solvers, mesh shape, time "
* Trivial to create scripts which vary parameters

Disadvantages:
* Object creation code is much more involved (it has to be able to handle
  all cases rather than just one case).
* ??ds anything else?

Computation speed is not an issue: time to run factories is trivial
compared to, e.g. Jacobian assembly, even for the smallest of problems.



Use of generic output (via MyProblem class)
-----------------------------

All my problem classes derive from one generic problem class (which derives
from the Problem base class). This class handles output in a consistent
way, which allows easy parsing by a general purpose script (parse.py) and
ensures that any information needed for analysis, debugging or reproduction
of results is always available.

As with the use of factory methods this comes with a small complexity
increase.

Additionally the time taken (per step) to compute some of the outputs
(especially energies) may be slowing down explicit timesteppers because
they take so many very small steps. The solve and Jacobian assembly almost
certainly dominates time taken by implicit timesteppers.
