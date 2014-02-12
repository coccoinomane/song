DESCRIPTION OF SONG
===================

SONG (Second Order Non-Gaussianity) is a second-order Boltzmann code that computes the effect of non-linear dynamics on the observables of the Cosmic Microwave Background. The reason for writing SONG was not to provide a more accurate version of the already existing first-order Boltzmann codes, such as CAMB or CLASS. Rather, SONG is a tool that, given a cosmological model, provides predictions for "new" observables or probes that do not exist at first order, such as

* the intrinsic bispectrum of the CMB (Pettinari, Fidler et al. 2013),
* the angular power spectrum of the spectral distortions (Renaux-Petel, Pitrou, Fidler & Pettinari 2013),
* the power spectrum of the magnetic fields generated at recombination, and
* the angular power spectrum of the B-mode polarisation (Fidler, Pettinari et al. 2014).

So far, SONG only computes the intrinsic bispectrum. It is our intention to include the other effects in the near future. This task is achievable with a comparatively smaller effort, because all these observables can be built starting from the second-order transfer functions; as we shall describe in the rest of the chapter, SONG already implements the complex framework needed to compute the second-order transfer functions up to today.

SONG is able to compute the intrinsic bispectrum of the CMB in about 10 CPU-hours, which is roughly equivalent to 10 minutes on a 60-core machine or a few hours on a quad- core one. Once they are implemented, the other observables will take considerably less time, because they do not involve the computation of the non-separable bispectrum integral. These numbers have to be compared with the two weeks taken by CMBquick (Pitrou et al. 2010) and the few days needed by CosmoLib2nd (Huang & Vernizzi 2013) for a full bispectrum run. (Note that these are rough estimates based on private communications with the authors of the aforementioned codes.)

The structure of SONG is based on that of CLASS, a recently released first-order Boltzmann code (Blas, Lesgourgues & Tram 2011). In particular, SONG inherits the philosophy of CLASS, that is to provide an easy-to-use interface that builds on a modular and flexible internal structure. Special care is taken to avoid the use of hard-coded numerical values, or "magic numbers"; the physical and numerical parameters are controlled through two separate input files by the user, who needs to set only those parameters of their interest, the others taking default values. In writing SONG we have followed the principle of encapsulation, so that a programmer who wants to modify or add a feature to SONG has to “hack” the code only in a few localised portions of the source files. When in doubt, said programmer can resort to the internal documentation, that comprises more than 10,000 lines of comments.

We conclude with a summary of the most relevant properties of SONG:

* SONG is written in C using only freely distributed libraries.
* It inherits from CLASS a modular and flexible structure (work is in progress to implement a Python interface, also adapted from the one used by CLASS).
* It employs an ad hoc differential evolver designed for stiff systems to solve the BES. 123
* It is OpenMP parallelised.
* Its source code is extensively documented with more than 10.000 lines of comments.
* It uses novel algorithms for Bessel convolution, bispectrum integration and 3D inter- polation.
* It implements the concept of beta-moments, whereby the non-realitivistic and relativistic species are treated in a unified way in terms of the moments of the distribution function.
