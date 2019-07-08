# ebids
Python code to translate electrophysiological data to BIDS format.

## Installation

Currently macOS/Linux only.

Add the ebids binary directory to your path. For example, if ebids is in your home directory under analysis:
```bash
export PATH=$PATH:$HOME/analysis/ebids/bin
```

Also add ebids to your python path:
```bash
export PYTHONPATH=$PYTHONPATH:$HOME/analysis/ebids
```

To run those commands automatically, put them in `$HOME/.bashrc`.

TODO: support for standard package installation.

If using python for analysis, you'll also want pybids: `pip install pybids`. If using Matlab, try [bids-matlab](https://github.com/bids-standard/bids-matlab).

## Converting behavioral data

For xmaze data, in the Terminal run:
```bash
xmaze-convert raw_dir bids_dir sub ses
```

where `raw_dir` is the path to the raw behavioral data for one session, `bids_dir` is the root directory of the bids data, sub is a subject code like `DS1902`, and ses is a session code like `01`.

## Reading BIDS data

For BIDS compatibility, in addition to running the relevant data conversion scripts, you also need `dataset_description.json` and `participants.tsv` files in the `bids_dir`.

### Python

In python, you should then be able to validate the BIDS formatting and load information about the dataset. At this point, you can easily do things like load behavioral data for a task run.

```python
from bids import BIDSLayout
layout = BIDSLayout(bids_dir) # validate dataset and load information
data = layout.get(task='xmaze', subject='DS1902', session='01', run=1).get_df() # load run events
data.loc[data.trial_type == 'reward',:] # get just the reward events
```

### MATLAB

Matlab support seems to be less well-developed, but also works:

```matlab
layout = bids.layout(bids_dir); # load dataset information
data = struct2table(layout.subjects(1).func(1).meta); # load run events
```
