---
title: "CRELES: Migration models"
date: "October 2021"
author: "Paola Arguello Pascualli"
output: 
  html_document: 
    keep_md: true
---

# Migration models
A new github repository to host the different analysis of the demographic history of Costa Rica was created and its located here: **https://github.com/santiago1234/DemogrHistoryCostaRica**. Santiago shared his previous codes on the different models he tried for his analysis, his codes are located here: https://github.com/santiago1234/mxb-genomes/tree/main/analysis-doc/210621-TRACTS-3pops

We will test the different models for the (1) whole CRELES cohort, for (2) only the Nicoyans and for (3) all the non-Nicoyans/

## Tracts: ppx_xxp

The first model we tried is the simplest one which has a first pulse of NAT and EUR ancestry and then a second impulse of AFR ancestry. The code and the raw results can be found here: *https://github.com/santiago1234/mxb-genomes/tree/main/analysis-doc/210531-TRACTS-ppx_xxp*



### Nicoyans


```r
library(tidyverse)
```

```
## Warning in as.POSIXlt.POSIXct(Sys.time()): unknown timezone 'zone/tz/2021a.3.0/
## zoneinfo/America/Mexico_City'
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
## ✓ tibble  3.1.4     ✓ dplyr   1.0.7
## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
## ✓ readr   2.0.2     ✓ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(scales)
```

```
## 
## Attaching package: 'scales'
```

```
## The following object is masked from 'package:purrr':
## 
##     discard
```

```
## The following object is masked from 'package:readr':
## 
##     col_factor
```

```r
library(ggthemes)
theme_set(theme_tufte(base_family = "Helvetica"))
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/1_Model_ppx_xxp/results_Nicoyans/fits.csv")
```

```
## Rows: 153 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (2): Ancestry, model
## dbl (7): bins, dat, pred, boot, ci_l, ci_u, loglkl
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
fit$cM <- fit$bins * 100
colors <- c("#e41a1c", "#4daf4a", "#377eb8")
fit %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, alpha = .9, color = "black") +
  geom_rangeframe(color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.03, 0)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(
    x = "Tract Length (cM)",
    y = "Relative frequency"
  ) +
  theme(
    legend.position = c(.85, .7),
    legend.box.background = element_rect(colour = "black", size = 0.1)
  )
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
##ggsave("fit.pdf", height = 3, width = 4)
##Collapse
```

### No-Nicoyans

```r
library(tidyverse)
library(scales)
library(ggthemes)
theme_set(theme_tufte(base_family = "Helvetica"))
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/1_Model_ppx_xxp/results_Non-Nicoyans/fits.csv")
```

```
## Rows: 153 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (2): Ancestry, model
## dbl (7): bins, dat, pred, boot, ci_l, ci_u, loglkl
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
fit$cM <- fit$bins * 100
colors <- c("#e41a1c", "#4daf4a", "#377eb8")
fit %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, alpha = .9, color = "black") +
  geom_rangeframe(color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.03, 0)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(
    x = "Tract Length (cM)",
    y = "Relative frequency"
  ) +
  theme(
    legend.position = c(.85, .7),
    legend.box.background = element_rect(colour = "black", size = 0.1)
  )
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
##ggsave("fit.pdf", height = 3, width = 4)
##Collapse
```



## Tracts: ppx_xxp_pxx

This model describes two pulses of NAT and EUR ancestry, followed by a AFR pulse and a final discrete pulse of EUR ancestry. 
Two codes  needed to run in order to complete the analysis are outlined below:
 
1. Fit the model


```bash

### ------------------------ 1_fit_ppx_xxp_pxx.py ------------------------------ ###

#!/usr/bin/env python

import os
import sys
import scipy
import numpy
import pylab
# the path to tracts if not in your default pythonpath
tractspath = "/data/users/parguello/LocalAncestryInference/modules/tracts-python3"  
sys.path.append(tractspath)
import tracts
import ppx_xxp_pxx


from warnings import warn

# optimization method. Here we use brute force; other options are implemented
# in tracts but are not implemented in this driver script.
method = "brute"

# demographic models to use

# ppx_xxp_fix has an initial pulse of migration from populations 0 and 1,
# followed by a pulse from population 2. "fix" referes to the fact that the
# migration rates in the model are fixed to the observed global ancestry
# proportions--then, we only have to optimize the timing of the migrations.
func = ppx_xxp_pxx.ppx_xxp_pxx_fix

# this function keeps track of whether parameters are in a "forbidden" region:
# whether mproportions are between 0 and 1, times positive, etc.
bound = ppx_xxp_pxx.outofbounds_ppx_xxp_pxx_fix

# defines the values of parameters to loop over in the brute force
# step. Times are in units of 100 generations.
slices = (
        slice(.08, .14, .01),  # range to explore for init_Eu
        slice(.08, .14, .01),  # " tstart
        slice(0.04, .1, .02),  # " afam_time
        slice(0.04, .1, .02)   # " nuEu_time
    )


# absolute bounds that parameters are not allowed to cross: times
# must be between 0 and 100 generations.
bounds = [(0, 1), (0, 1)]

# choose order of populations: the labels in our ancestry files are
# strings, and we need to tell tracts which string correspond to which
# population in the model. Here the population labels are (somewhat
# confusingly) numbers that do not match the order in the population.
# "labels" will tell tracts that model population 0 has label '0', model
# population 1 has label '2', and model population 2 has label '1' in the
# local ancestry files.
labels = ['EUR', 'AMR', 'AFR']

# directories in which to read and write
directory = str(sys.argv[1])
##"/data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/CRELES_all/"
outdir = str(sys.argv[2])
## "results/"

if not os.path.exists(outdir):
    os.makedirs(outdir)


# string between individual label and haploid chromosome id in input file
inter = "_CRELES"
# string at the end of input file. Note that the file should end in ".bed"
end = "_cM.bed"


# Get a list of all individuals in directory.
_files = os.listdir(directory)
files = [file
        for file in _files
        if file.split('.')[-1] == "bed"]  # only consider bed files

# Get unique individual labels
names = list(set(file.split('_CRELES_')[0] for file in files))

if len(_files) != len(files):
    warn("some files in the bed directory were ignored, since they do not "
            "end with `.bed`.")

# Load the population using the population class's constructor. It
# automatically iterates over individuals and haploid copies (labeled _A"
# and "_B" by default
pop = tracts.population(names=names, fname=(directory, inter, end))

# Rather than creating a new population for each bootstrap instance, we
# just replace the list of individuals to iterate over. We need to save a
# copy of the initial list of individuals to do this!
indivs = pop.indivs

# generate the histogram of tract lengths
(bins, data) = pop.get_global_tractlengths(npts=50)

data = [data[poplab] for poplab in labels]

# Calculate ancestry proportions
props = list(pop.get_mean_ancestry_proportions(labels))


Ls = pop.Ls
nind = pop.nind

cutoff = 2  # ignore particular bins

bootorder = range(len(indivs))


xopt = tracts.optimize_brute_fracs2(
    bins, Ls, data, nind, func, props, slices, outofbounds_fun=bound, cutoff=cutoff)
print(xopt)
optmod = tracts.demographic_model(func(xopt[0], props))
optpars = xopt[0]
liks = xopt[1]
maxlik = optmod.loglik(bins, Ls, data, pop.nind, cutoff=cutoff)

expects = []
for popnum in range(len(data)):
    expects.append(optmod.expectperbin(Ls, popnum, bins))

expects = nind * numpy.array(expects)

bootnum = 0
outf = outdir + "boot%d_%2.2f" % (bootnum, maxlik,)

fbins = open(outf + "_bins", 'w')
fbins.write("\t".join(map(str, bins)))
fbins.close()

fliks = open(outf + "_liks", 'w')
fliks.write("\t".join(map(str, liks)))
fliks.close()

fdat = open(outf + "_dat", 'w')
for popnum in range(len(data)):
    fdat.write("\t".join(map(str, data[popnum])) + "\n")

fdat.close()
fmig = open(outf + "_mig", 'w')
for line in optmod.mig:
    fmig.write("\t".join(map(str, line)) + "\n")

fmig.close()
fpred = open(outf + "_pred", 'w')

for popnum in range(len(data)):
    fpred.write("\t".join(map(str, pop.nind * numpy.array(optmod.expectperbin(Ls, popnum, bins)))) + "\n")


fpred.close()
fpars = open(outf + "_pars", 'w')

fpars.write("\t".join(map(str, optpars)) + "\n")
fpars.close()

# bootstrap order in case we need to rerun/check something.
ford = open(outf + "_ord", 'w')
ford.write("\t".join(map(lambda i: "%d" % (i,), bootorder)))
ford.close()

```
    
2. Collect the data
    

```bash

import pandas as pd
import numpy as np
import sys
sys.path.append("/data/users/parguello/LocalAncestryInference/modules/")
import tractsresults


subpop = str(sys.argv[1])
## "CRELES"
dir_to_res = str(sys.argv[2])
## "/data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Nicoyans/"

# same order as it was given to fit the data
labels = ['EUR', 'AMR', 'AFR']
# load the PUR data
res = tractsresults.load_boot_rest(dir_to_res, boot_n=0, labels=labels)
res['population'] = subpop

# load the estimated parameters
pars = tractsresults._load_boot_file(dir_to_res, filetype="pars")
pars = pars.round(4) * 100
pars = "init_Eu = {}, tstart = {}\nafam_time = {}, nuEu_time = {}".format(*pars)
res['parameters'] = pars

res.to_csv("results.csv", index=False)

```

3. Migration model plot


```r
library(tidyverse)
library(scales)
library(ggthemes)
plot.migration.plot <- function(fit, plot_title){
  
  theme_set(theme_tufte(base_family = "Helvetica"))
  fit$cM <- fit$bins * 100
  colors <- c("#e41a1c", "#4daf4a", "#377eb8")
  
  fit %>% 
    ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
    geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
    geom_line(show.legend=FALSE) +
    geom_line(show.legend=FALSE, color = "black", size = 0.2) +
    geom_point(aes(y=dat),shape = 21, alpha = .9, color = "black") +
    geom_rangeframe(color = "black") +
    ggtitle(plot_title) +
    scale_y_log10(
      oob = scales::squish_infinite,
      expand = c(0, 0),
      breaks = 10^(-1:5),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(expand = c(.03, 0)) +
    coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(
      x = "Tract Length (cM)",
      y = "Relative frequency"
    ) +
    theme(
      legend.position = c(.85, .7),
      legend.box.background = element_rect(colour = "black", size = 0.1)
    )

}
```


### Nicoyans

#### Migration model for only Nicoyans


```bash
## Nicoyans
cd /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/
python3 1_fit_ppx_xxp_pxx.py  /data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/Nicoyans/ /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Nicoyans/

cd /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Nicoyans/
python3 ../2_collect_results.py Nicoyan /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Nicoyans/

```

#### Plot for Nicoyans


```r
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/2_Model_ppx_xxp_pxx/results_Nicoyans/results.csv")
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 153 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (3): Ancestry, population, parameters
## dbl (6): bins, dat, pred, boot, ci_l, ci_u
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
plot.migration.plot(fit, "ppx_xxp_pxx: Nicoyans")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


### Non-Nicoyans

#### Migration models for Non-Nicoyans

```bash

## Non-Nicoyans
cd /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/
python3 1_fit_ppx_xxp_pxx.py  /data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/Non-Nicoyans/ /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Non-Nicoyans/

cd /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Non-Nicoyans/
python3 ../2_collect_results.py No_Nicoyan /data/users/parguello/LocalAncestryInference/2_Model_ppx_xxp_pxx/results_Non-Nicoyans/

```


#### Plot for Non-Nicoyans


```r
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/2_Model_ppx_xxp_pxx/results_Non-Nicoyans/results.csv")
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 153 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (3): Ancestry, population, parameters
## dbl (6): bins, dat, pred, boot, ci_l, ci_u
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
plot.migration.plot(fit, "ppx_xxp_pxx: Non-Nicoyans")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-9-1.png)<!-- -->





## Tracts: ccx_xxp

This model describes a continuos contribution of NAT and EUR ancestry, followed by a AFR pulse.
Two codes  needed to run in order to complete the analysis are outlined below:
 
1. Fit the model
    

```bash

#!/usr/bin/env python

import os
import sys
import scipy
tractspath = "/data/users/parguello/LocalAncestryInference/modules/tracts-python3"  # the path to tracts if not in your default pythonpath
sys.path.append(tractspath)
import tracts
import ccx_xxp
import numpy
import pylab

from warnings import warn

# demographic models to use

# ppx_xxp_fix has an initial pulse of migration from populations 0 and 1,
# followed by a pulse from population 2. "fix" referes to the fact that the
# migration rates in the model are fixed to the observed global ancestry
# proportions--then, we only have to optimize the timing of the migrations.
func = ccx_xxp.ccx_xxp

# this function keeps track of whether parameters are in a "forbidden" region:
# whether mproportions are between 0 and 1, times positive, etc.
bound = ccx_xxp.outofbounds_ccx_xxp


# absolute bounds that parameters are not allowed to cross: times
# must be between 0 and 100 generations.
bounds = [(0, 1), (0, 1)]

# choose order of populations: the labels in our ancestry files are
# strings, and we need to tell tracts which string correspond to which
# population in the model. Here the population labels are (somewhat
# confusingly) numbers that do not match the order in the population.
# "labels" will tell tracts that model population 0 has label '0', model
# population 1 has label '2', and model population 2 has label '1' in the
# local ancestry files.
labels = ['EUR', 'AMR', 'AFR']

# directories in which to read and write
directory = str(sys.argv[1])
##"/data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/CRELES_all/"

outdir = str(sys.argv[2])
## "/data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_CRELES/"

if not os.path.exists(outdir):
    os.makedirs(outdir)

# string between individual label and haploid chromosome id in input file
inter = "_CRELES"
# string at the end of input file. Note that the file should end in ".bed"
end = "_cM.bed"


# Get a list of all individuals in directory.
_files = os.listdir(directory)
files = [file
        for file in _files
        if file.split('.')[-1] == "bed"]  # only consider bed files

# Get unique individual labels
names = list(set(file.split('_CRELES_')[0] for file in files))

if len(_files) != len(files):
    warn("some files in the bed directory were ignored, since they do not "
            "end with `.bed`.")

# Load the population using the population class's constructor. It
# automatically iterates over individuals and haploid copies (labeled _A"
# and "_B" by default
pop = tracts.population(names=names, fname=(directory, inter, end))

# Rather than creating a new population for each bootstrap instance, we
# just replace the list of individuals to iterate over. We need to save a
# copy of the initial list of individuals to do this!
indivs = pop.indivs

# generate the histogram of tract lengths
(bins, data) = pop.get_global_tractlengths(npts=50)

data = [data[poplab] for poplab in labels]

# Calculate ancestry proportions
props = list(pop.get_mean_ancestry_proportions(labels))

Ls = pop.Ls
nind = pop.nind

cutoff = 2  # ignore particular bins

bootorder = range(len(indivs))


# the starting parameters for this models
startparams = [
    0.04,  # frac1
    0.04,  # frac2
    0.16,  # t1
    0.06,  # frac3
    0.10,  # t2
]

xopt = tracts.optimize_cob_fracs2(startparams,bins,Ls,data,nind,func,props,outofbounds_fun=bound,cutoff=cutoff,epsilon=1e-2)

print(xopt)

optmod = tracts.demographic_model(func(xopt))
optpars = xopt
maxlik = optmod.loglik(bins, Ls, data, pop.nind, cutoff=cutoff)

expects = []
for popnum in range(len(data)):
    expects.append(optmod.expectperbin(Ls, popnum, bins))

expects = nind * numpy.array(expects)

bootnum = 0
outf = outdir + "boot%d_%2.2f" % (bootnum, maxlik,)

fbins = open(outf + "_bins", 'w')
fbins.write("\t".join(map(str, bins)))
fbins.close()


fdat = open(outf + "_dat", 'w')
for popnum in range(len(data)):
    fdat.write("\t".join(map(str, data[popnum])) + "\n")

fdat.close()
fmig = open(outf + "_mig", 'w')
for line in optmod.mig:
    fmig.write("\t".join(map(str, line)) + "\n")

fmig.close()
fpred = open(outf + "_pred", 'w')
for popnum in range(len(data)):
    fpred.write(
        "\t".join(map(str, pop.nind * numpy.array(optmod.expectperbin(Ls, popnum, bins)))) + "\n")
fpred.close()

fpars = open(outf + "_pars", 'w')

fpars.write("\t".join(map(str, optpars)) + "\n")
fpars.close()

# bootstrap order in case we need to rerun/check something.
ford = open(outf + "_ord", 'w')
ford.write("\t".join(map(lambda i: "%d" % (i,), bootorder)))
ford.close()

```
    
    
2. Collect results
    

```bash

import pandas as pd
import numpy as np
import sys
sys.path.append("/data/users/parguello/LocalAncestryInference/modules/")
import tractsresults

subpop = str(sys.argv[1])
## CRELES
dir_to_res = str(sys.argv[2])
## "/data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_CRELES/"

# same order as it was given to fit the data
labels = ['EUR', 'AMR', 'AFR']

# load the PUR data
res = tractsresults.load_boot_rest(dir_to_res, boot_n=0, labels=labels)
res['population'] = subpop

# load the estimated parameters
pars = tractsresults._load_boot_file(dir_to_res, filetype="pars")
pars = pars.round(4) * 100
pars = "init_Eu = {}, tstart = {}\nafam_time = {}, nuEu_time = {}".format(*pars)
res['parameters'] = pars

res.to_csv("results.csv", index=False)

# ancestry proportions

#migmat = np.loadtxt("results/boot0_-312.37_mig")
#ancp = tractsresults.ancestry_prop_over_time(migmat, labels)
#ancp.to_csv("results/ancprop.csv", index=False)
#pulses = tractsresults.migration_pulses(migmat, labels)
#pulses.to_csv("results/migpulses.csv", index=False)

```
    
    
3. Plot the migration model
    
    
    ```r
    library(tidyverse)
    library(scales)
    library(ggthemes)
    
    plot.migration.plot <- function(fit, plot_title){
      
      theme_set(theme_tufte(base_family = "Helvetica"))
      fit$cM <- fit$bins * 100
      colors <- c("#e41a1c", "#4daf4a", "#377eb8")
      
      fit %>% 
    ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
    geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
    geom_line(show.legend=FALSE) +
    geom_line(show.legend=FALSE, color = "black", size = 0.2) +
    geom_point(aes(y=dat),shape = 21, alpha = .9, color = "black") +
    geom_rangeframe(color = "black") +
    ggtitle(plot_title) +
    scale_y_log10(
      oob = scales::squish_infinite,
      expand = c(0, 0),
      breaks = 10^(-1:5),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(expand = c(.03, 0)) +
    coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(
      x = "Tract Length (cM)",
      y = "Relative frequency"
    ) +
    theme(
      legend.position = c(.85, .7),
      legend.box.background = element_rect(colour = "black", size = 0.1)
    )
    
    }
    ```

### Nicoyan

#### Migration model

```bash

cd /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp
python3 1_fit_ccx_xxp.py /data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/Nicoyans/ /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_Nicoyans/

cd /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_Nicoyans
python3 ../2_collect_results.py Nicoyan /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_Nicoyans/

```


#### Plot migration

```r
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/3_Model_ccx_xxp/results_Nicoyans/Nicoyans_results.csv")
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 153 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (3): Ancestry, population, parameters
## dbl (6): bins, dat, pred, boot, ci_l, ci_u
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
plot.migration.plot(fit, "ccx_xxp: Nicoyans")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


  
### Non-Nicoyan
#### Migration model

```bash

cd /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp
python3 1_fit_ccx_xxp.py /data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/Non-Nicoyans/ /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_Non-Nicoyans/

cd /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_Non-Nicoyans
python3 ../2_collect_results.py Non_Nicoyan /data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_Non-Nicoyans/

```


#### Plot migration

```r
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/3_Model_ccx_xxp/results_Non-Nicoyans/Non-Nicoyans_results.csv")
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 153 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (3): Ancestry, population, parameters
## dbl (6): bins, dat, pred, boot, ci_l, ci_u
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
plot.migration.plot(fit, "ccx_xxp: Non-Nicoyans")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


    
    

## Tracts: xxp_ccp 

This model was created by Santiago and describes an initial AFR pulse (xxp), followed by a continuous migration of NAT and EUR and a second AFR pulse (ccp). 

1. Fit the model


```r
#!/usr/bin/env python

import os
import sys
import scipy
tractspath = "/data/users/parguello/LocalAncestryInference/modules/tracts-python3/"  # the path to tracts if not in your default pythonpath
sys.path.append(tractspath)
import tracts
modelpath = "/data/users/parguello/LocalAncestryInference/4_Model_xxp_ccp/"
sys.path.append(modelpath)
import ccp
import numpy
import pylab

from warnings import warn

# demographic models to use

# ppx_xxp_fix has an initial pulse of migration from populations 0 and 1,
# followed by a pulse from population 2. "fix" referes to the fact that the
# migration rates in the model are fixed to the observed global ancestry
# proportions--then, we only have to optimize the timing of the migrations.
func = ccp.ccp

# this function keeps track of whether parameters are in a "forbidden" region:
# whether mproportions are between 0 and 1, times positive, etc.
bound = ccp.outofbounds_ccp

# local ancestry files.
labels = ['AFR', 'AMR', 'EUR']


# directories in which to read and write
directory = "/data/users/parguello/LocalAncestryInference/0_CRELES_Tracts_input/CRELES_all/"
outdir = "/data/users/parguello/LocalAncestryInference/4_Model_xxp_ccp/results_CRELES/"
os.makedirs(outdir)

inter = "_CRELES"
end = "_cM.bed"


# Get a list of all individuals in directory.
_files = os.listdir(directory)
files = [file
        for file in _files
        if file.split('.')[-1] == "bed"]  # only consider bed files

# Get unique individual labels
names = list(set(file.split('_CRELES')[0] for file in files))

if len(_files) != len(files):
    warn("some files in the bed directory were ignored, since they do not "
            "end with `.bed`.")

# Load the population using the population class's constructor. It
# automatically iterates over individuals and haploid copies (labeled _A"
# and "_B" by default
pop = tracts.population(names=names, fname=(directory, inter, end))

# Rather than creating a new population for each bootstrap instance, we
# just replace the list of individuals to iterate over. We need to save a
# copy of the initial list of individuals to do this!
indivs = pop.indivs

# generate the histogram of tract lengths
(bins, data) = pop.get_global_tractlengths(npts=51)

data = [data[poplab] for poplab in labels]

# Calculate ancestry proportions
props = list(pop.get_mean_ancestry_proportions(labels))

Ls = pop.Ls
nind = pop.nind

cutoff = 3  # ignore particular bins

bootorder = range(len(indivs))


# the starting parameters for this models

startparams = [
    0.001,  # frac1
    0.01,  # frac2
    2.05,  # frac3
    0.096,  # T
    0.005
]

xopt = tracts.optimize_cob_fracs2(startparams,bins,Ls,data,nind,func,props,outofbounds_fun=bound,cutoff=cutoff,epsilon=1e-2)

print(xopt)

optmod = tracts.demographic_model(func(xopt))
optpars = xopt
maxlik = optmod.loglik(bins, Ls, data, pop.nind, cutoff=cutoff)

expects = []
for popnum in range(len(data)):
    expects.append(optmod.expectperbin(Ls, popnum, bins))

expects = nind * numpy.array(expects)

bootnum = 0
outf = outdir + "boot%d_%2.2f" % (bootnum, maxlik,)

fbins = open(outf + "_bins", 'w')
fbins.write("\t".join(map(str, bins)))
fbins.close()


fdat = open(outf + "_dat", 'w')
for popnum in range(len(data)):
    fdat.write("\t".join(map(str, data[popnum])) + "\n")

fdat.close()
fmig = open(outf + "_mig", 'w')
for line in optmod.mig:
    fmig.write("\t".join(map(str, line)) + "\n")


fmig.close()
fpred = open(outf + "_pred", 'w')
for popnum in range(len(data)):
    fpred.write(
        "\t".join(map(str, pop.nind * numpy.array(optmod.expectperbin(Ls, popnum, bins)))) + "\n")
fpred.close()

fpars = open(outf + "_pars", 'w')

fpars.write("\t".join(map(str, optpars)) + "\n")
fpars.close()

# bootstrap order in case we need to rerun/check something.
ford = open(outf + "_ord", 'w')
ford.write("\t".join(map(lambda i: "%d" % (i,), bootorder)))
ford.close()
```

2. Collect results


```r
import pandas as pd
import numpy as np
import sys
sys.path.append("/data/users/parguello/LocalAncestryInference/modules/")
import tractsresults

subpop = str(sys.argv[1])
## CRELES
dir_to_res = str(sys.argv[2])
## "/data/users/parguello/LocalAncestryInference/3_Model_ccx_xxp/results_CRELES/"

# same order as it was given to fit the data
labels = ['EUR', 'AMR', 'AFR']

# load the PUR data
res = tractsresults.load_boot_rest(dir_to_res, boot_n=0, labels=labels)
res['population'] = subpop

# load the estimated parameters
pars = tractsresults._load_boot_file(dir_to_res, filetype="pars")
pars = pars.round(4) * 100
pars = "init_Eu = {}, tstart = {}\nafam_time = {}, nuEu_time = {}".format(*pars)
res['parameters'] = pars

res.to_csv("results.csv", index=False)

# ancestry proportions

migmat = np.loadtxt("boot0_-707.52_mig") ## This has to be changed with the mig file created on the fit step 
ancp = tractsresults.ancestry_prop_over_time(migmat, labels)
ancp.to_csv("ancprop.csv", index=False)
pulses = tractsresults.migration_pulses(migmat, labels)
pulses.to_csv("migpulses.csv", index=False)
```



3. Plot the migration model


### Nicoyans


```r
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/4_Model_xxp_ccp/results_Nicoyan/results.csv")
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 156 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (3): Ancestry, population, parameters
## dbl (6): bins, dat, pred, boot, ci_l, ci_u
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
plot.migration.plot(fit, "xxp_ccp: Non-Nicoyans")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-19-1.png)<!-- -->


### Non-Nicoyans



```r
fit <- read_csv("/Users/paola.arguello/Documents/GitHub/DemogrHistoryCostaRica/analysis/4_Model_xxp_ccp/results_Non-Nicoyan/results.csv")
```

```
## Warning: One or more parsing issues, see `problems()` for details
```

```
## Rows: 156 Columns: 9
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (3): Ancestry, population, parameters
## dbl (6): bins, dat, pred, boot, ci_l, ci_u
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
plot.migration.plot(fit, "xxp_ccp: Non-Nicoyans")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Transformation introduced infinite values in continuous y-axis

## Warning: Transformation introduced infinite values in continuous y-axis
```

![](Tracts_analysis_recap_files/figure-html/unnamed-chunk-20-1.png)<!-- -->
