# Divergence_in_tensors
Alexander Litvinenko, RWTH Aachen, Germany
#
THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY
#
This project contains numerical examples for this article "Computing f-Divergences and Distances of High-Dimensional Probability Density Functions -- Low-Rank Tensor Approximations", available at https://arxiv.org/abs/2111.07164
#
The user needs, first, to install TT-Toolbox from here https://github.com/oseledets/TT-Toolbox .

Then to run TT-Toolbox/setup.m

Only after that the user can run single Matlabs files from this repository
#
The following functionality is available:
1. Tensor train (TT) approximation of a probability characteristic function (pcf)
2. TT approximation of a probability density function (pdf)
3. Computation of the log() function of a pdf and then the KL distance
4. Computation of the square root of a pdf and then the Hellinger distance
#

For example:

my_pdf_tensor_ex1.m was used for Example 6.1 and to generate data from Table 5, i.e.
for the computation of $D_{H}(\alpha_1,\alpha_2)$ between two \alpha-stable distributions (\alpha=1.5 and \alpha=0.9) for different AMEn tolerances. 

my_H2dist_tensor.m was used to compute $D_{H}(\alpha_1,\alpha_2)$ between two \alpha - stable distributions for different dimensions d and resolutions n.


my_KLD_tensor.m was used for computing $D_KL(\alpha_1,\alpha_2)$ between two \alpha-stable distributions for various \alpha with fixed d=8 and n=64.
 
Plese write me an email if you have any questions litvinenko@uq.rwth-aachen.de
