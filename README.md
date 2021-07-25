# Exam assignment
[![Documentation Status](https://readthedocs.org/projects/cmepda20-giovanni/badge/?version=latest)](https://cmepda20-giovanni.readthedocs.io/en/latest/?badge=latest)

This is the project for the exam of computing methods for experimental physics. The scope of this work is to find the Higgs boson in the decay channel H->ZZ->4l with CMS Open data.

The python file is in python dir.

usage: exam_assignment.py [-h] [-nofast] [-local] [-time]

Program that find the Higgs boson in the decay channel H->ZZ->4l with CMS Open data.

optional arguments:
  -h, --help  show this help message and exit
  -nofast     No-fast mode take data from raw file and it does the data selection. This saves also data selected to
              data path.
  -local      Local mode take data from raw file in local and it does the data selection.
  -time       Print the execution time.
