# Markov-chain Monte Carlo Ground-Motion Selection Algorithms for Conditional Intensity Measure Targets

The code at this repository provides an implementation of the algorithms described in:

Shi Y and Stafford PJ (2018). "Markov-chain Monte Carlo Ground-Motion Selection Algorithms for Conditional Intensity Measure Targets", _Earthquake Engineering & Structural Dynamics_ (currently in press)

The code here presents an illustrative example of _Case 1_ from the above article.

The script ``` Main_records_selection.m ``` should be used as the entry point for the code while the main work is performed within ``` MCMCGMS_method1_for_case_study_SA_PLUS_DS595_upload``` and ```MCMCGMS_method2_for_case_study_SA_PLUS_DS595_upload``` for algorthims MCMCGMS method 1 and MCMCGMS method 2, respectively.
The scripts aim to be self-documenting, but the original manuscript should be consulted for the mathematical details of each step.

This code is intended to assist readers in coding record selection algorithms using MCMCGMS methods for their own purposes. It is simple to adjust this example code to consider other intensity measures cases not already considered in these example files.
