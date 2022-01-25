# type-II-A-leader-repeat

This repository has code relevant to the paper "The leader RNA prioritizes newly-acquired spacers for CRISPR-Cas9 immunity".  The scrips were used to generate Figures 

## Installation

## About the data

Data under the directory `data/original` came from Omer Alkhnbashi's predictions of the relevant type.  Each line corresponds to one system and consists of three tab-delimited fields.  The first field is an identifier for the system.  The second field is the 180 nucleotides upstream of the predicted first repeat (i.e., the potential leader).  The third field is the sequence of the predicted first repeat.

As described in the paper, this data was split into training and test datasets, which are available in the directory `data/training-test`.  Each file name describes the sub-type and whether the data is training or test data.  For example, the file `II-A--training.tab` is the training data for type II-A systems.  These files are constructed using the script `py/ChaseTrainingValidationTestData.py` based on the data in `data/original`, such that the leader of each system within a (training or test) dataset is no more than 70% identical to another system, and each leader sequence in the training dataset is at most 50% identical to any leader in the test dataset.  The format of each dataset file is the same as the files in `data/original`. There are just fewer lines.

## About the scripts

## Example commands

