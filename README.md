# NCOurd

### Expectation-maximization (EM) algorithm to infer the length distribution of Non-crossover (NCO) events from gene conversion tracts.

NCOurd is implemented in python.

## Table of Contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [Citation](#citation)
* [Contact](#contact)

## Requirements

* python3:
    * numpy
    * pandas
    * scipy

To install python packages:
```
pip install numpy pandas scipy
```
## Installation

Clone the Git repository

```
git clone https://github.com/DecodeGenetics/NCOurd.git 
cd ncourd
```

## Usage

### Running the EM model

```
python bin/EMmodel.py --help
```

displays the command line interface:
```
usage: EMmodel.py [-h] [--detfile DETFILE] [--ngroups NGROUPS1]
                  [--maxiter MAXITER] [--verbose] [--ntracts NTRACTS1]
                  [--nskip NSKIP] [--multiprocess] [--paramfile PARAMFILE]
                  [--method METHOD] [--dclass DCLASS] [--minmeth MINMETH]
                  [--maxspan MAXSPAN] [--nbinom_loc NBINOM_LOC]
                  [--randomsample RANDOMSAMPLE] [--accwts ACCWTS]
                  [--cubicgrid] [--twidx TWIDX] [--mink MINK]
                  tractfiles [tractfiles ...]

Run EM model on tractsets

positional arguments:
  tractfiles            tract weight function files to use

optional arguments:
  -h, --help            show this help message and exit
  --detfile DETFILE     file containing detection function
  --ngroups NGROUPS     number groups to fit - Default is #tractfiles
  --maxiter MAXITER     maximum number of iterations - Default is 1000
  --verbose             more printing
  --ntracts NTRACTS     number of tracts to take from each tractfile
  --nskip NSKIP         number of tracts to skip from each tractfile
  --multiprocess        use multiprocessing
  --paramfile PARAMFILE
                        File containing initial parameters for the EM model
  --method METHOD       Simultaneously "simul" or separately "sep" evaluate
                        shape and scale
  --dclass DCLASS       Class of distributions excluding the default negative
                        binomial distributions "nbinom" gamma "gamma" and
                        geometric "geom" distributions are supported
  --minmeth MINMETH     Method to pass to scipy.optimize.minimize for
                        simultaneous evaluation must support bounds
  --maxspan MAXSPAN     Discard tracts which are at least this long
  --nbinom_loc NBINOM_LOC
                        If nbiom distributions are used this sets location
                        parameter of nbiom distributions - Default is 1
  --randomsample RANDOMSAMPLE
                        Sample random tracts from input tracts with given seed
  --accwts ACCWTS       Accelerate weights for the first N iterations
  --cubicgrid           use cubic wgrid instead of linear for integration
  --twidx TWIDX         indexes to use in the tractfile
  --mink MINK           minimum shape parameter for the equivalent gamma
                        distributions (variance is at least mean**2*k)
```

For usage examples see: 
```Makefile```

### Input files for EMmodel.py
####```tractfiles```:
Consist of a line for each gene conversion tracts.
Each line contains 2800 numerical values of the tract function at the "k-grid".
The tractfiles arguments are read with ```numpy.loadtxt```

####```--detfile DETFILE```:
Tab delimited values with header having fields "k" and "p".
"k" is an integer values and
"p" is the value of the detection function at "k".
The detfile argument is read with ```pandas.read_csv```
and interpolated with ```scipy.interpolate.interp1d``` with ```kind="cubic"``` 

### Generating detection function file from informative marker set

```
python bin/Markers2DF.py --help
```

displays the command line interface:

```
usage: Markers2DF.py [-h] [--regions REGIONS] [--dense] [--penetrance PENETRANCE] 
                     [--output OUTPUT] markerfile

Calculate detection function from informative markers writing the output to
[markerfile].p[PENETRANCE].detdf - The output file is suitable as input file for
EMmodel.py

positional arguments:
  markerfile            informative marker file

optional arguments:
  -h, --help            show this help message and exit
  --regions REGIONS     file with regions where the detecion functions 
                        should be calculated
  --dense               use denser grid
  --penetrance PENETRANCE
                        probability of NCO marker being gene converted
  --output OUTPUT       name of output file
```

#### Input file for Markers2DF.py

Tab separated values with header having columns "Chr" and "pos"

- "Chr" is a label for the chromosome containing the marker
- "pos" is the genomic position of the marker

### Generating tract function file from tract file

```
python bin/tracts2twdet.py --help
```

displays the command line interface:

```
usage: tracts2twdet.py [-h] [--penetrance PENETRANCE] tracts

calculate tract functions from tracts and writes them to the
file [tracts].p[PENETRANCE].twdet.gz - The output file is
suitable as input for EMmodel.py

positional arguments:
  tracts                file containing the tracts

optional arguments:
  -h, --help            show this help message and exit
  --penetrance PENETRANCE
                        probability of NCO marker becoming
                        gene converted
```

#### Input file for tracts2twdet.py

Tab separated values with header having columns "Chr", "lbd", "ubd", "lgap" and "rgap"

- "Chr" is a label for the chromosome containing the marker
- "lbd" is the genomic position of the first marker
- "ubd" is the genomic position of the last marker
- "lgap" is a comma separated list of distances from "lbd" to 1st,2nd,3rd,...,20th informative marker to the left
- "rgap" is a comma separated list of distances from "ubd" to 1st,2nd,3rd,...,20th informative marker to the right


## Citation

Coming soon

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.


