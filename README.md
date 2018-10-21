# MetabolicNetworks
Files for the work on metabolic networks with Leonid, Tamon, Cedric, Reza, Hooman, and Nafiseh

## Prerequistes

You need to have licensed MATLAB before running the scripts. 

Have CellNetAnalyzer in a folder with the same address of this repository folder in your computer.

Test and run FluxModeCalculator and make sure you have all its requirements installed.

If you want to convert SBML format to any other format you need to set MongooseGUI3 and install its requirements.

## Minimcal cut-set scripts
There are implemenations of four approaches for enumerating MCSs here. You can run either of them by running its script followed by model name followed by target reaction followed by 1 or 0 which indicates prior network reduction. -1 stands for all reactions.
```
$./runDual.sh BIOMD0000000004 1 0
```
This command enumerate all the MCSs for target reaction one on the uncompressed stoichiometry matrix of BIOMD0000000004 with Dual approach.

```
$./runRowSpace.sh BIOMD0000000005 -1 1
```
This command find all the MCSs in the compressed model BIOMD0000000005 for every target reaction with RowSpace approach.

The output would be CSV table format in the `log` folder with each row belongs to a target reaction containig time of each process and number of MCSs enumerated. You will see the total and average time at the last row if you commanded to enumerate MCSs for every target reaction.

## Converting formats
With a stoichiometry matrix and a boolean reversiblity vector you can create an SBML file. All you have to do run the function `network_to_sbml` and give it these inputs and your chosen model name.

By running `SBML_to_All_other_format.py` in MongooseGUI3 folder you will have many other formats of the network you may need. Some of these formats are needed for running the scripts.

## Superset removal

`SuperSetRemoval.jar` is written in java to remove superset and repeated sets and return the remaining sets.

## Other implementations
Flux mode calculator, CellNetAnalayzer, and MoongooseGUI3 are implementations done by other people and was added here.
