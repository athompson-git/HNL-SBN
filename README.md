# HNLs at DUNE, DarkQuest, and the SBN Program

[![DOI](https://zenodo.org/badge/265690633.svg)](https://doi.org/10.5281/zenodo.14894889)

Code written for the production of $Z^\prime$ and Heavy Neutral Leptons (HNL) for [https://arxiv.org/abs/2410.08981](https://arxiv.org/abs/2410.08981), "Enhancing the Sensitivity to Seesaw Predictions in Gauged $B−L$ Scenarios", Francesco Capozzi, Bhaskar Dutta, Gajendra Gurung, Wooyoung Jang, Ian M. Shoemaker, Adrian Thompson, and Jaehoon Yu

Contact: Adrian Thompson ```a.thompson@northwestern.edu```

### Requirements:
* alplib: [https://github.com/athompson-git/alplib](https://github.com/athompson-git/alplib)
* git LFS (for large flux files)

### Files include
* HNL branching ratios: ```hnl_decays.py```
* Zprime branching ratios: ```zprime_decay_widths.py```
* Flux simulation: ```hnl_flux.py```


### Constants
* ```darkquest_constants.py```
* ```dune_constants.py```


### Data
* R ratio: ```data/r_ratio_by_sqrts-GeV.txt```

We also acknowledge and use constraints on the HNL parameter space from [https://www.hep.ucl.ac.uk/~pbolton/index.html](https://www.hep.ucl.ac.uk/~pbolton/index.html)