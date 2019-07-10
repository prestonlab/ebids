# ebids
Code to translate electrophysiological data to BIDS format.

## Installation

### eBIDS

Currently macOS/Linux only.

Add the ebids binary directory to your path. For example, if ebids is in your home directory under analysis:
```bash
export PATH=$PATH:$HOME/analysis/ebids/bin
```

Also add ebids to your python path:
```bash
export PYTHONPATH=$PYTHONPATH:$HOME/analysis/ebids
```

To install the Matlab scripts, add `ebids/matlab` to your Matlab path.

To run those commands automatically, put them in `$HOME/.bashrc`.

### Dependencies

To work with Neuralynx data in python (used for some preprocessing scripts), download our [fork of the Neo project](https://github.com/prestonlab/python-neo) and install by adding the project directory to your PYTHONPATH (similar to above for ebids).

If using python for analysis, you'll also want pybids: `pip install pybids`. If using Matlab, try [bids-matlab](https://github.com/bids-standard/bids-matlab).

## Converting behavioral data

### XMaze

In the Terminal run:
```bash
xmaze-convert raw_dir bids_dir sub ses
```

where `raw_dir` is the path to the raw behavioral data for one session, `bids_dir` is the root directory of the bids data, sub is a subject code like `DS1902`, and ses is a session code like `01`.

### eCoupling

In Matlab, run:
```matlab
ecoupling_convert(raw_dir, bids_dir, sub, ses);
```
The inputs are similar to XMaze. The `raw_dir` should contain a file called `header.mat`.

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

```octave
layout = bids.layout(bids_dir); % load dataset information
data = struct2table(layout.subjects(1).func(1).meta); % load run events
```
