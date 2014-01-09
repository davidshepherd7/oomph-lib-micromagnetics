# Run the real script and exit with it's error status
./validate.py &> Validation/validation.log
exit $?
