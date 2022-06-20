# Worflow to install the R package from github using conda

1. Tag the repository
```
git tag -a v0.1
git push origin --tags
```

2. Create conda recipe

Suggestion: use a conda environment to build packages
```
conda create -n build-env conda-build
conda activate build-env
mamba install -c conda-forge r-base
mamba install -c conda-forge r-essentials
```

Then,
```
conda skeleton cran git@github.com:PatZeis/scATAC.Peak.annotation.git
```

This creates a directory named `r-scatac.peak.annotation` which contains three files: `meta.yaml`, `build.sh` and `bld.bat`.

**NOTES**

Two minor modifications to the skeletof file `meta.yaml` are required before building the package. (i) Modify license entry by adding quotes, as YAML cannot parse backticks:
```
license: "`use_mit_license()`, `use_gpl3_license()` or friends to pick a license"
```
(ii) change `r-genomicranges` by `bioconductor-genomicranges`

3. Build conda package
```
conda build -c conda-forge -c bioconda r-scatac.peak.annotation --R=<r_version>
```

4. Install package
```
mamba install -c /<CONDA_ENV_PATH>/conda-bld/  -c conda-forge -c bioconda r-scatac.peak.annotation
```
