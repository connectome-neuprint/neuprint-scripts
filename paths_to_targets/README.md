# Paths to targets

**paths_to_targets** finds all paths from a set of source bodies to a set of downstream target bodies.
This how to install the necessary packages to run the script and walks through a few examples.


## Installation

This installation requires conda package manager. If you do not have conda installed you can find the install guide 
[here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda).


First, either create or activate your desired conda environment. Then install the python requirements in
`requirements.txt` and install **neuprint-python** from the **flyem-forge** channel.

```
source activate <env>
conda install --file requirements.txt --yes
conda install -c flyem-forge neuprint-python --yes
```


## Running the script

The easiest way to run **paths_to_targets.py** is using a TOML file. [job.toml](job.toml) is a mock
example connecting neurons 1, 2, and 3 to the FB (Fan-Shaped Body). The options are split into two tables: `[client]` 
and `[job]`. The `[client]` table contains the parameters relating to the neuprint-python client.  The `[job]` table 
contains parameters relating to generation of the paths.  Below is a description of each parameter in the TOML file.

- client
    - **server** (required): Neuprint server URL (https://emdata1.int.janelia.org:11000/).
    - **token** (required): Neuprint authentication token.
    
- job
    - **sources** (required): List of body IDs from which the paths start.
    - **depth** (required): Maximum length of the paths.  Loops and autapses excluded.
    - **connection_threshold** (optional): 
    - **targets** (optional): List of body IDs where paths should end. Either this value or **roi_target** should be defined. 
    - **roi_target** (optional): ROI from which to retrieve target body IDs. Either this value or **targets** should be defined. 
    - **roi_threshold** (optional): Minimum combined number of pre and postsynaptic sites for bodies in **roi_target**.
    - **roi_pre_threshold** (optional): Minimum number of presynaptic sites for bodies in **roi_target**.
    - **roi_post_threshold** (optional): Minimum number of postsynaptic sites for bodies in **roi_target**.
    - **graph_file** (optional): File where the connectivity graph information should be stored.
    - **path_file** (optional): File where all path permutations are stored.
    - **verbose** (optional): If true, or excluded, print runtime information.  If false, run silently 


## Examples
Once the TOML file has been configure, the script can be ran using the option `-f`.

```bash
python paths_to_targets.py -f job.toml
```

All the parameters that can be set in the TOML file can also be set in the command line.  Any parameters set in the
command line will override the corresponding value from the TOML file. Below is an example using the same TOML file, but
with varying connection thresholds and output file names.


```bash
python paths_to_targets.py -f job.toml -ct 10 -p path10.csv
python paths_to_targets.py -f job.toml -ct 20 -p path20.csv
```

For a full list of command line options and their description run the script with the options `-h or --help`.

```bash
$ python paths_to_targets.py -h
$ python paths_to_targets.py --help
```

