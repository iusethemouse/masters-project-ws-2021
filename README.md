# A Comparison of Neural Network-based CRN Abstraction Methods: StochNetV2 and DeepCME
Master's project by Ivan Prigarin, University of Konstanz, Winter Semester 21/22

## Project aims
- Examine the performance and applicability of StochNetV2 [1] to additional example CRNs presented in the DeepCME paper [2].
- Implement and evaluate a method for estimating species concentration moments using StochNetV2, and compare the results with those of DeepCME.
- Implement and evaluate a method for estimating parameter sensitivities of the moment estimations using StochNetV2, and compare the results with those of DeepCME.

## Pipelines
In this repository, you can find the following Jupyter Notebooks and additional files:
- `1-stochnetv2-training.ipynb`: the pipeline to train a StochNetV2 neural network on a given CRN.
- `2-moment-estimation.ipynb`: the pipeline to estimate the first and second moments for a given CRN using a trained StochNetV2 neural network.
- `3-parameter-sensitivity-estimation.ipynb`: the pipeline to estimate the parameter sensitivities of the first and second moment estimations for a given CRN using a trained StochNetV2 neural network.
- `crn-definitions`: a folder containing the definitions of example CRNs using the `Gillespy` framework.
- `summary.pdf`: a set of slides outlining the findings and conclusions of the project.

## References
[1] StochNetV2: A Tool for Automated Deep Abstractions for Stochastic Reaction Networks, Denis Repin, Nhat-Huy Phung, Tatjana Petrov, QEST 2020
[2] DeepCME: A deep learning framework for computing solution statistics of the chemical master equation by Ankit Gupta, Christoph Schwab and Mustafa Khammash, PLoS Computational Biology 2021
