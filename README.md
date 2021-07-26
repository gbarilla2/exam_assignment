# Exam assignment
[![Documentation Status](https://readthedocs.org/projects/higgs-to4l-analysis/badge/?version=latest)](https://higgs-to4l-analysis.readthedocs.io/en/latest/?badge=latest)


This is the project for the exam of computing methods for experimental physics. The scope of this work is to find the Higgs boson in the decay channel H->ZZ->4l with CMS Open data.

This analysis is based on this [article](https://arxiv.org/abs/1202.1997).

## Requirements

Python3 and [PyRoot](https://root.cern/manual/python/)

## Usage

```bash
python3 -i -m higgs_to4l_analysis [-h] [-nofast] [-local] [-time]
```

Program that find the Higgs boson in the decay channel H->ZZ->4l with CMS Open data.

optional arguments:
  -h, --help  show this help message and exit
  -nofast     No-fast mode take data from raw file and it does the data selection. This saves also data selected to
              data path.
  -local      Local mode take data from raw file in local and it does the data selection.
  -time       Print the execution time.

## Object and Event selection

The events are selected by applying the following cuts:
  -4 particle for event (4&mu;,4e,2&mu;2e)  
  -
  -
  -
