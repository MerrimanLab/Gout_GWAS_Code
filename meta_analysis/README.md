# Meta-analysis scripts

This directory contains the main scripts used to run single- and trans-ancestry
meta-analysis.

The overall process of the meta-analysis pipeline is as follows:
1. Clean up the raw input summary stats, before running proper QC
1. Run QC on per-study basis first, then across-studies (this is done within ancestry)
1. Run single-ancestry meta-analysis with the QCed input summary stats
1. Clean up the single-ancestry meta-analysis results and process it for trans-ancestry meta-analysis
1. Run trans-ancestry meta-analysis and clean up the results

Refer to the `README` documents in the relevant directory for specifics.

