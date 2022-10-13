### Compiling the code:

For the compilation to succeed, one must specify the path to the MilliDAQ shared object file (libMilliDAQ.so) and `ConfigurationReader.cc` in the compile script 
`compile.sh`.
This involves editing `setup.sh` to have the correct environmental variables. The OFFLINESITE is used to define the site for the database

To compile the code, which produces an executable file (from the macro) that will be used for the tree-making, execute the following command

```bash
. setup.sh
./compile.sh (target macro name).exe

```

Keep track of the "target macro name" as this macro will be needed when making the offline trees.

##To use common scripts, log onto milliqan username on SL6 machine (e.g. cms1, cms3, cms6, cms29)

From the milliqanOffline directory you can run on a single file using the runMakeTree.py script.
### Making the Offline Trees
The code is run by executing the `runMakeTree.py` script. You may execute the `runMakeTree.py` script using the following command:
```bash
   python3 scripts/runOfflineFactory.py --inFile "input_filename".root --outputFile "output_filename".root --exe ./(target macro name).exe (--publish) ...
```
Where `(target macro name).exe` is the one specified when running the `compile.sh` script. The publish option will include the file details in the mongoDB database.