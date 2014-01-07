#!/bin/bash

micromagroot="../.."
CONTROL_SCRIPTS="$micromagroot/control_scripts"
MM_LIB_DIR="$micromagroot/"


# source common commands for energy tests
. "$micromagroot/self_test/micromag_energy/init_validate.sh"


# In this test we check that an exact block preconditioner is exactly the
# same as the Jacobian by checking that gmres converges in one iteration.


# No magnetostatics
# ============================================================

valdir=$TPWD/Validation_no_ms
new_clean_dir $valdir

cd $CONTROL_SCRIPTS/driver/
./driver llg -dt 0.001 -tmax 0.0009 -ref 0 -mesh sq_square -solver gmres -prec blockexact \
    -ts midpoint-bdf -happ minus_z -initm z -output-jac always -disable-ms \
    -outdir $valdir 2>&1 > $valdir/stdout

# Check that gmres converges in one step (i.e. block-exact == exact
# solution of J)
cut -d\; -f6 < $valdir/trace | tail -n1 > $valdir/iters
wrapped_fpdiff $valdir/iters <(echo "[1]")

# Check that we got the expected number of blocks (49)
wrapped_fpdiff <(ls -l $valdir/ | grep "J_1_1_block_" | wc -l) <(echo "49")

# Check that phi blocks are empty
wrapped_fpdiff  $valdir/J_1_1_block_0_0 <(echo -n "")
wrapped_fpdiff  $valdir/J_1_1_block_1_1 <(echo -n "")
wrapped_fpdiff  $valdir/J_1_1_block_5_5 <(echo -n "") # boundary phis
wrapped_fpdiff  $valdir/J_1_1_block_6_6 <(echo -n "")


# explicit magnetostatics but implicit llg
# ============================================================
# this should be exactly the same as without ms

valdir=$TPWD/Validation_expl_ms
new_clean_dir $valdir

cd $CONTROL_SCRIPTS/driver/
./driver llg -decoupled-ms -dt 0.001 -tmax 0.0009 -ref 0 -mesh sq_square \
    -solver gmres -prec blockexact \
    -ts midpoint-bdf -happ minus_z -initm z -output-jac always \
    -outdir $valdir 2>&1 > $valdir/stdout

# Check that gmres converges in one step (i.e. block-exact == exact
# solution of J)
cut -d\; -f6 < $valdir/trace | tail -n1 > $valdir/iters
wrapped_fpdiff $valdir/iters <(echo "[1]")

# Check that we got the expected number of blocks (49)
wrapped_fpdiff <(ls -l $valdir/ | grep "J_1_1_block_" | wc -l) <(echo "49")

# Check that phi blocks are empty
wrapped_fpdiff  $valdir/J_1_1_block_0_0 <(echo -n "")
wrapped_fpdiff  $valdir/J_1_1_block_1_1 <(echo -n "")
wrapped_fpdiff  $valdir/J_1_1_block_5_5 <(echo -n "") # boundary phis
wrapped_fpdiff  $valdir/J_1_1_block_6_6 <(echo -n "")


# With magnetostatics (implicit, CR Jacobian)
# ============================================================

valdir=$TPWD/Validation_ms
new_clean_dir $valdir

cd $CONTROL_SCRIPTS/driver/
./driver llg -dt 0.001 -tmax 0.0009 -ref 0 -mesh sq_square -solver gmres -prec blockexact \
    -ts midpoint-bdf -happ minus_z -initm z -output-jac always \
    -outdir $valdir 2>&1 > $valdir/stdout

# Check that gmres converges in one step (i.e. block-exact == exact
# solution of J), two Newton steps expected here.
cut -d\; -f6 < $valdir/trace | tail -n1 > $valdir/iters
wrapped_fpdiff $valdir/iters <(echo "[1, 1]")

# Check that we got the expected number of blocks
wrapped_fpdiff <(ls -l $valdir/ | grep "J_1_1_block_" | wc -l) <(echo "49")

# Check that phi blocks are non-empty
[[ -s $valdir/J_1_1_block_0_0 ]]
[[ -s $valdir/J_1_1_block_1_1 ]]
[[ -s $valdir/J_1_1_block_5_5 ]]
[[ -s $valdir/J_1_1_block_6_6 ]]


# With m blocking (implicit, CR Jacobian)
# ============================================================

valdir=$TPWD/Validation_ms
new_clean_dir $valdir

cd $CONTROL_SCRIPTS/driver/
./driver llg -dt 0.001 -tmax 0.0009 -ref 0 -mesh sq_square -solver gmres -prec blockexact \
    -ts midpoint-bdf -happ minus_z -initm z -output-jac always \
    -blocking group-m \
    -outdir $valdir 2>&1 > $valdir/stdout

# Check that gmres converges in one step (i.e. block-exact == exact
# solution of J), two Newton steps expected here.
cut -d\; -f6 < $valdir/trace | tail -n1 > $valdir/iters
wrapped_fpdiff $valdir/iters <(echo "[1, 1]")

# Check that we got the expected number of blocks
wrapped_fpdiff <(ls -l $valdir/ | grep "J_1_1_block_" | wc -l) <(echo "25")

# Check that phi blocks are non-empty
[[ -s $valdir/J_1_1_block_0_0 ]]
[[ -s $valdir/J_1_1_block_1_1 ]]
[[ -s $valdir/J_1_1_block_3_3 ]]
[[ -s $valdir/J_1_1_block_4_4 ]]


# ??ds not working yet
# # With magnetostatics (implicit, SOM Jacobian)
# # ============================================================

# valdir=$TPWD/Validation_ms
# new_clean_dir $valdir

# cd $CONTROL_SCRIPTS/driver/
# ./driver llg -dt 0.001 -tmax 0.0009 -ref 0 -mesh sq_square -solver som-gmres \
#     -prec som-main-blockexact \
#     -ts midpoint-bdf -happ minus_z -initm z -implicit-ms -output-jac always \
#     -outdir $valdir 2>&1 > $valdir/stdout

# # Check that gmres converges in one step (i.e. block-exact == exact
# # solution of J), two Newton steps expected here.
# cut -d\; -f6 < $valdir/trace | tail -n1 > $valdir/iters
# wrapped_fpdiff $valdir/iters <(echo "[1, 1]")


exit 0
