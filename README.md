# type-II-A-leader-repeat

ZLH : see ZLH

This repository has code relevant to the paper "Spacer prioritization in CRISPR-Cas9 immunity is enabled by the leader RNA".  The scrips were used to generate Figures 

## Installation

The scripts require Python version 3, a recent version of Perl and the bash interpreter.  The only non-standard Python library is the Python bindings for the [ViennaRNA library](https://www.tbi.univie.ac.at/RNA/)

## About the data in this repository

Data under the directory `data/original` came from Omer Alkhnbashi's predictions of the relevant type.  Each line corresponds to one system and consists of three tab-delimited fields.  The first field is an identifier for the system.  The second field is the 180 nucleotides upstream of the predicted first repeat (i.e., the potential leader).  The third field is the sequence of the predicted first repeat.

As described in the paper, this data was split into training and test datasets, which are available in the directory `data/training-test`.  Each file name describes the sub-type and whether the data is training or test data.  For example, the file `II-A--training.tab` is the training data for type II-A systems.  These files are constructed using the script `py/ChaseTrainingValidationTestData.py` based on the data in `data/original`, such that the leader of each system within a (training or test) dataset is no more than 70% identical to another system, and each leader sequence in the training dataset is at most 50% identical to any leader in the test dataset.  The format of each dataset file is the same as the files in `data/original`. There are just fewer lines.

Our initial test dataset for type II-C systems is in the directory `data/training-test/old/old-II-C--test.tab`. Subsequent analysis suggested that some of these systems might not be correctly classified, so we conservatively removed them the test dataset (see file `data/training-test/II-C--test.tab`).  Neither dataset have significant aggregate p-values.

## About the scripts

The scripts are almost exclusively written in Python.  The scripts are prefixed with "Chase", which was my code for this project. In summary, there are the following:

- ChaseExtractKeyFieldsWithActualValue.py : This script is not necessary to reproduce the paper's results.  It was for visualizing specific statistics, and for drawing MFE structures for visualization.  The drawing 
- ChaseFigureFromFolding.py : Makes blue-dot figures, like in ZLH
- ChaseFolding.py : Main script to calculate folding statistics for a set of systems.  Aggregate p-values are calculated in `ChaseSummarize.py`.
- ChaseGC.py : This script is not necessary to reproduce the paper's results.  It was for looking at the G+C content of leaders and repeats.
- ChaseLibrary.py : This file provides functions used by other scripts.  It cannot be run directly.
- ChaseRemoveBadII-C.py : Script to remove type II-C systems which might have been incorrectly classified.
- ChaseSummarize.py : Summarize the data created by ChaseFolding.py to compute aggregate p-values and draw histograms for all statistics.
- ChaseTrainingValidationTestData.py : Cluster a set of systems (using one of the files in `data/original`
- ChaseVarious.sh
- ChaseVariousVarious.pl


## Example commands

