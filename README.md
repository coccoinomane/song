## UPDATE: SONG is looking for a maintainer

Find out more at

* https://github.com/coccoinomane/song/issues/10

## DESCRIPTION

SONG computes the non-linear evolution of the Universe in order to predict cosmological observables such as the bispectrum of the Cosmic Microwave Background (CMB).

More precisely, SONG (Second Order Non-Gaussianity) is a second-order Boltzmann code, as it solves the Einstein and Boltzmann equations up to second order in the cosmological perturbations; the physics, mathematics and numerics of SONG are described extensively in my publicly available [PhD thesis][10], especially in Chapters 5 and 6.

The reason for writing SONG was not to provide a more accurate version of the already existing first-order Boltzmann codes, such as [CAMB][1] or [CLASS][2]. Rather, SONG is a tool that, given a cosmological model, provides predictions for "new" observables or probes that do not exist at first order, such as

* the intrinsic bispectrum of the CMB ([[3]] and [[4]]),
* the angular power spectrum of the spectral distortions [[5]],
* the power spectrum of the magnetic fields generated around and after recombination [[11]],
* the angular power spectrum of the B-mode polarisation [[6]],
* the Fourier bispectrum of the cold dark matter density fluctuations.

SONG's structure is based on that of [CLASS][1], a linear Boltzmann code introduced in 2013 [[9]]. In particular, SONG inherits the philosophy of CLASS, that is to provide an easy-to-use interface built on a modular and flexible internal structure. Special care is taken to avoid the use of hard-coded numerical values, or "magic numbers"; the physical and numerical parameters are controlled through two separate input files by the user, who needs to set only those parameters of their interest, the others taking default values.

In writing SONG we have followed the principle of encapsulation, so that if you want to modify or add a feature to SONG, you have to edit the code in few localised portions. When in doubt, you can resort to the internal documentation, which comprises more than 10,000 lines of comments.


## MAIN FEATURES OF SONG

Here is a summary of the most relevant properties of SONG:

* SONG is an open-source project written in C using only freely distributed libraries.
* It is fast: it takes about an hour on a modern laptop to compute the intrinsic bispectrum at 10% precision.
* It employs an ad-hoc differential evolver designed for stiff systems [[9]].
* It fully includes CMB polarisation at second order.
* It implements the tight coupling approximation for faster computation.
* It implements a Fisher matrix module to quantify the observability of any intrinsic or primordial bispectrum.
* It uses novel algorithms for Bessel convolution, bispectrum integration and 3D interpolation.
* It implements the concept of beta-moments, whereby the non-realitivistic and relativistic species are treated in a unified way in terms of the moments of the distribution function.
* It inherits from CLASS a modular and flexible structure.
* It is OpenMP parallelised.
* Its source code is extensively documented with more than 10,000 lines of comments.

## DOCUMENTATION
SONG's source code is extensively documented with more than 10,000 lines of comments. The physical and numerical steps are explained in the source (.c) files, while the employed variables are documented in the header (.h) files.  We are also working on a series of tutorials which I plan to release early in 2016.

The physics, mathematics and numerics of SONG are described extensively in my PhD thesis [[10]], especially in Chapters 5 and 6. For any doubt or enquiry, please email me at <guido.pettinari@gmail.com>.


## GETTING STARTED
First, download SONG. You can do so either from the project's page <https://github.com/coccoinomane/song/releases> or from the command line with:

    git clone --recursive https://github.com/coccoinomane/song.git song.git

The `--recursive` flag is important because it allows to download both SONG and CLASS in one command. To compile and make a test run of SONG, enter SONG's directory and execute the following commands from your terminal:

    make song
    ./song ini/intrinsic.ini pre/quick_song_run.pre

The first command compiles the code using the instructions contained in the file `makefile`. This can be customised to match the configuration of your system. The main variables are the location of the C compiler (gcc by default) and of the OpenMP library (if you want parallelisation); SONG does not rely on any other external library.

The second command executes a test run of SONG using the input files `ini/intrinsic.ini` and `pre/quick_song_run.pre`. These are text-only files with a list of "key = value" settings; SONG will read their content and set its internal parameters accordingly.

To make a SONG run that will produce an actual physical result, try with

    ./song ini/intrinsic.ini pre/sn_pol_10percent_lmax2000.pre
    
SONG will output to file and to screen the Fisher matrix and signal-to-noise of the intrinsic bispectrum, with a moderate angular resolution, and a precision of about 10% for both temperature and polarisation; the calculation will take about 40 minutes on an 8-core machine.

Feel free to experiment with the parameter files! For example, if you delete `eBisp` from `ini/intrinsic.ini`, SONG will neglect E-polarisation and run much faster. For a guide on what each parameter does, please refer to the (*not yet*) documented file `explanatory.ini` or to the documented `ini/intrinsic.ini`. Use these files as templates for creating your custom input files!


## PARALLEL COMPUTING

SONG is a parallel code, via the OpenMP standard. Therefore, it can use all the cores in a single node, but it cannot span multiple nodes.

Set the number of cores you want to use on the command line via the environment variable `OMP_NUM_THREADS`. For example, if you run SONG on a laptop with 8 cores on a bash shell, you might want to execute `export OMP_NUM_THREADS=8` before running SONG.

Apple Os X default compiler, `clang`, does not natively support OpenMP. I suggest using the GNU `gcc` compiler, which can be downloaded from [Macports], [Homebrew] or [HPC].


## DIRECTORY STRUCTURE
The directory structure of SONG is important to learn how the code works:

* The 'source' directory contains the main source files in C. Each file corresponds to a module in SONG.

* The 'tools' directory contains accessory source files in C with purely numerical functions or utility functions.

* The 'main' directory contains the main source files, i.e. the executable files, including song.c.

* The 'python' directory contains Python scripts to inspect some of SONG outputs (credits to Thomas Tram).

* The 'include' directory contains the declaration files (.h) for all the C files in the 'source', 'main' and 'tools' directories.

* The 'test' directory contains executable programs to test the outputs of SONG.

* The 'scripts' directory contains bash and gnuplot scripts to run SONG iteratively and to plot its results.

* The 'ini' and 'pre' directories contain, respectively, parameter and precision files that can be fed to SONG.

* The 'output' folder, initially empty, will contain the products of SONG computations, usually in the form of text files.


## CONTRIBUTE!
Please feel free to contribute to the development of SONG! You can do so via the Github project page: <https://github.com/coccoinomane/song.git>. Fork the repository, make your modifications and send a pull request.

The [master branch] has the latest stable SONG release, while the [develop branch] contains the more advanced (and potentially buggy) development version.

If you are not familiar with Git and Github, please consider sending me the proposed modifications via email at <guido.pettinari@gmail.com>. Otherwise, there are very good [Github tutorials] available.

If you do not have a feature in mind, but nonetheless want to contribute to SONG, you are very welcome to do so! Just have fun addressing the many TODO comments in the code :-)


## CONTACT
Feel free to send an email to Guido Pettinari <guido.pettinari@gmail.com> or Christian Fidler <christian.fidler@port.ac.uk> for any enquiry or suggestion. We realise that bispectra and second-order perturbation theory may be a bit confusing, so we are happy to help.


## CITATIONS
Writing SONG took more than four years. If you found SONG useful, please acknowledge our effort and that of our collaborators by citing one or more of the following references. If in doubt, please use the first reference.

* G. W. Pettinari, C. Fidler, R. Crittenden, K. Koyama, and D. Wands. "The intrinsic bispectrum of the cosmic microwave background". J. Cosmology Astropart. Phys., 04(2013)003, doi: 10.1088/1475-7516/2013/04/003 [[3]]
```
@article{pettinari:2013a,
       author = {{Pettinari}, G.~W. and {Fidler}, C. and {Crittenden}, R. and {Koyama}, K. and {Wands}, D.},
        title = "{The intrinsic bispectrum of the cosmic microwave background}",
      journal = {J. Cosmology Astropart. Phys.},
         year = 2013,
        month = apr,
       volume = 4,
          eid = {003},
        pages = {3},
          doi = {10.1088/1475-7516/2013/04/003},
       adsurl = {http://adsabs.harvard.edu/abs/2013JCAP...04..003P},
archivePrefix = "arXiv",
       eprint = {1302.0832}
}
```

* G. W. Pettinari, "The intrinsic bispectrum of the cosmic microwave background [PhD Thesis]". Springer Theses Series, Volume XXIII, 2015 [[10]]
```
@book{pettinari:2015a,
	   Author = {Pettinari, G.~W.},
	Publisher = {Springer International Publishing},
	   Series = {Springer Theses},
	    Title = {The Intrinsic Bispectrum of the Cosmic Microwave Background},
	      Url = {http://www.springer.com/gp/book/9783319218816},
	     Year = {2015},
          doi = {10.1088/1475-7516/2013/04/003},
         Isbn = {9783319218823},
archivePrefix = "arXiv",
       eprint = {1405.2280}
}
```

* G. W. Pettinari, C. Fidler, R. Crittenden, K. Koyama, A. Lewis & D. Wands. "Impact of polarization on the intrinsic cosmic microwave background bispectrum". Phys. Rev. D., 90, 103010, doi: 10.1103/PhysRevD.90.103010 [[4]]
```
@article{pettinari:2014a,
	   author = {Pettinari, G.~W. and Fidler, C. and Crittenden, R. and Koyama, K. and Lewis, A and Wands, D.},
	    issue = {10},
	  journal = {Phys. Rev. D},
	    month = {Nov},
	 numpages = {6},
	    pages = {103010},
	publisher = {American Physical Society},
	    title = {Impact of polarization on the intrinsic cosmic microwave background bispectrum},
	      url = {http://link.aps.org/doi/10.1103/PhysRevD.90.103010},
	   volume = {90},
	     year = {2014},
archivePrefix = "arXiv",
       eprint = {1406.2981},
          doi = {10.1103/PhysRevD.90.103010}
}
```


[1]: http://camb.info/ "The CAMB Boltzmann code"
[2]: http://class-code.net/ "The CLASS Boltzmann code"
[3]: http://arxiv.org/abs/1302.0832 "The intrinsic bispectrum of the Cosmic Microwave Background"
[4]: http://arxiv.org/abs/1406.2981 "Impact of polarisation on the intrinsic CMB bispectrum"
[5]: http://arxiv.org/abs/1312.4448 "Spectral distortions in the cosmic microwave background polarization"
[6]: http://arxiv.org/abs/1401.3296 "The intrinsic B-mode polarisation of the Cosmic Microwave Background"
[7]: http://www2.iap.fr/users/pitrou/cmbquick.htm "The CMBquick 2nd-order Boltzmann code"
[8]: http://arxiv.org/abs/1212.3573 "The CMB bispectrum from recombination"
[9]: http://arxiv.org/abs/1104.2933 "The Cosmic Linear Anisotropy Solving System (CLASS)"
[10]: http://arxiv.org/abs/1405.2280 "The intrinsic bispectrum of the Cosmic Microwave Background [PhD thesis]"
[11]: http://arxiv.org/abs/1511.07801 "A precise numerical estimation of the magnetic field generated around recombination"
[Macports]: https://www.macports.org
[Homebrew]: http://brew.sh
[HPC]: http://hpc.sourceforge.net
[Github tutorials]: https://help.github.com/articles/good-resources-for-learning-git-and-github/ "Good Resources for Learning Git and GitHub"
[Github page]: https://github.com/coccoinomane/song.git
[master branch]: https://github.com/coccoinomane/song/tree/master
[develop branch]: https://github.com/coccoinomane/song/tree/develop
