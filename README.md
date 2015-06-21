## DESCRIPTION

SONG computes the non-linear evolution of the Universe in order to predict cosmological observables such as the bispectrum of the Cosmic Microwave Background (CMB).

More precisely, SONG (Second Order Non-Gaussianity) is a second-order Boltzmann code, as it solves the Einstein and Boltzmann equations up to second order in the cosmological perturbations; the physics, mathematics and numerics of SONG are described extensively in my publicly available [PhD thesis][10], especially in Chapters 5 and 6.

The reason for writing SONG was not to provide a more accurate version of the already existing first-order Boltzmann codes, such as [CAMB][1] or [CLASS][2]. Rather, SONG is a tool that, given a cosmological model, provides predictions for "new" observables or probes that do not exist at first order, such as

* the intrinsic bispectrum of the CMB ([[3]] and [[4]]),
* the angular power spectrum of the spectral distortions [[5]],
* the power spectrum of the magnetic fields generated at recombination, and
* the angular power spectrum of the B-mode polarisation [[6]].

The current version of SONG computes the intrinsic bispectrum and its correlations with any primordial bispectrum. It is our intention to include the other effects in the near future. This task is achievable with a comparatively smaller effort, because all these observables can be built starting from the second-order transfer functions, and SONG already implements the framework needed to compute them up to today.

SONG is able to compute the intrinsic bispectrum of the CMB in about 10 CPU-hours, which is roughly equivalent to 10 minutes on a 60-core machine or a few hours on a quad- core one. Once they are implemented, the other observables will take considerably less time, because they do not involve the computation of the non-separable bispectrum integral. These numbers have to be compared with the two weeks taken by [CMBquick][7] and the few days needed by [CosmoLib2nd][8] for a full bispectrum run. (Note that these are rough estimates based on private communications with the authors of the aforementioned codes in 2013.)

The structure of SONG is based on that of [CLASS][1], a linear Boltzmann code introduced in 2013 [[9]]. In particular, SONG inherits the philosophy of CLASS, that is to provide an easy-to-use interface that builds on a modular and flexible internal structure. Special care is taken to avoid the use of hard-coded numerical values, or "magic numbers"; the physical and numerical parameters are controlled through two separate input files by the user, who needs to set only those parameters of their interest, the others taking default values. In writing SONG we have followed the principle of encapsulation, so that a programmer who wants to modify or add a feature to SONG has to "hack" the code only in a few localised portions of the source files. When in doubt, said programmer can resort to the internal documentation, that comprises more than 10,000 lines of comments.


## MAIN FEATURES OF SONG

Here is a summary of the most relevant properties of SONG:

* SONG is written in C using only freely distributed libraries.
* It inherits from CLASS a modular and flexible structure (work is in progress to implement a Python interface, also adapted from the one used by CLASS).
* It employs an ad hoc differential evolver designed for stiff systems [[9]] to solve the Einstein and Boltzmann equations.
* It is OpenMP parallelised.
* Its source code is extensively documented with more than 10,000 lines of comments.
* It uses novel algorithms for Bessel convolution, bispectrum integration and 3D interpolation.
* It implements a Fisher matrix module to quantify the observability of the any intrinsic or primordial bispectrum.
* It fully includes polarisation of the cosmic microwave background.
* It implements the concept of beta-moments, whereby the non-realitivistic and relativistic species are treated in a unified way in terms of the moments of the distribution function.


## DOCUMENTATION
SONG's source code is extensively documented with more than 10,000 lines of comments. The physical and numerical steps are explained in the source (.c) files, while the employed variables are documented in the header (.h) files.  We are also working on a manual which we plan to release by the end of 2015.

The physics, mathematics and numerics of SONG are described extensively in my PhD thesis [[10]], especially in Chapters 5 and 6. For any doubt or enquiry, please email me at <guido.pettinari@gmail.com>.


## GETTING STARTED
To compile and make a test run of SONG, execute the following commands from your terminal:

    make song
    ./song ini/intrinsic.ini pre/quick_song_run.pre

The first command compiles the code using the instructions contained in the file `makefile`. This can be customised to match the configuration of your system. The main variables are the location of the C compiler (gcc by default, the only compiler SONG has been tested on) and of the OpenMP library (if you want parallelisation); SONG does not rely on any other external library.

The second command executes a test run of SONG using the input files `ini/intrinsic.ini` and `pre/quick_song_run.pre`. These are text-only files with a list of "key = value" settings; SONG will read their content and set its internal parameters accordingly.

To make a SONG run that will produce an actual physical result, try with

    ./song ini/intrinsic.ini pre/sn_pol_10percent_lmax2000.pre
    
SONG will output to file and to screen the Fisher matrix and signal-to-noise of the intrinsic bispectrum, with a moderate angular resolution, and a precision of about 10% for both temperature and polarisation; the calculation will take about 40 minutes on an 8-core machine.

Feel free to experiment with the parameter files! For example, you can remove from `ini/intrinsic.ini` the E-polarisation to make SONG run much faster. For a guide on what each parameter does, please refer to the (not yet) documented file `explanatory.ini` or to the (yet) documented `ini/intrinsic.ini`. Use these files as templates for creating your custom input files!


## PARALLEL COMPUTING

SONG is a parallel code, via the OpenMP standard. Therefore, it can use all the cores in a single node, but it cannot span multiple nodes.

Set the number of cores you want to use on the command line via the environment variable `OMP_NUM_THREADS`. For example, if you run SONG on a laptop with 8 cores on a bash shell, you might want to execute `export OMP_NUM_THREADS=8` before running SONG.


## DIRECTORY STRUCTURE
The directory structure of SONG is important to learn how the code works:

* The 'source' directory contains the main source files in C. Each file corresponds to a module in SONG.

* The 'tools' directory contains accessory source files in C with purely numerical functions or utility functions.

* The 'main' directory contains the main source files, i.e. the executable files, including song.c.

* The 'python' directory contains Python scripts to launch SONG (not implemented yet).

* The 'include' directory contains the declaration files (.h) for all the C files in the 'source', 'main' and 'tools' directories.

* The 'test' directory contains executable programs to test the outputs of SONG.

* The 'scripts' directory contains bash, python and gnuplot scripts to run SONG efficiently or to process and plot SONG's results.

* The 'ini' and 'pre' directories contain, respectively, parameter and precision files that can be fed to SONG.

* The 'output' folder, initially empty, will containt the products of SONG computation usually in the form of text files.


## CONTRIBUTE
Please feel free to contribute to the development of SONG! You can do so via a pull request using the Github project page: <https://github.com/coccoinomane/song.git>.

If you are not familiar with Git and Github, please consider sending me the proposed modifications via email at <guido.pettinari@gmail.com>. Otherwise, there are very good [Github tutorials] available.

If you do not have a feature in mind, but nonetheless want to contribute to SONG, you are very welcome to do so! Just have fun addressing the many TODO comments in the code :-)


## CONTACT
Feel free to send an email to Guido Pettinari <guido.pettinari@gmail.com> or Christian Fidler <christian.fidler@port.ac.uk> for any enquiry or suggestion. We realise that bispectra and second-order perturbation theory may be a bit confusing, so we are happy to help.


[1]: http://camb.info/ "The CAMB Boltzmann code"
[2]: http://class-code.net/ "The CLASS Boltzmann code"
[3]: http://arxiv.org/abs/1302.0832 "The intrinsic bispectrum of the Cosmic Microwave Background"
[4]: http://arxiv.org/abs/1406.2981 "Impact of polarisation on the intrinsic CMB bispectrum"
[5]: http://arxiv.org/abs/1312.4448 "Spectral distortions in the cosmic microwave background polarization"
[6]: http://arxiv.org/abs/1401.3296 "The intrinsic B-mode polarisation of the Cosmic Microwave Background"
[7]: http://www2.iap.fr/users/pitrou/cmbquick.htm "The CMBquick 2nd-order Boltzmann code"
[8]: http://arxiv.org/abs/1212.3573 "The CMB bispectrum from recombination"
[9]: http://arxiv.org/abs/1104.2933 "The Cosmic Linear Anisotropy Solving System (CLASS)"
[10]: http://arxiv.org/abs/1405.2280 "The intrinsic bispectrum of the Cosmic Microwave Background (Ph.D. thesis)"
[Github tutorials]: https://help.github.com/articles/good-resources-for-learning-git-and-github/ "Good Resources for Learning Git and GitHub"
