<div align="center">
    <img src="artic-logo.png?raw=true?" alt="artic-logo" width="250">
    <h1>ARTIC - align_trim</h1>
    <h3>A standalone version of the artic fieldbioinformatics align_trim script for primer trimming and normalising viral amplicon sequencing</h3>
</div>

---

## Installation

### Via conda - !Not submitted yet!

```sh
conda install -c bioconda align_trim
```

### Via source

#### 1. downloading the source:

Download a [release](https://github.com/artic-network/align_trim/releases) or use the latest master (which tracks the current release):

```sh
git clone https://github.com/BioWilko/align_trim/
cd align_trim
```

#### 2. installing dependencies:

The `align_trim` has several [software dependencies](https://github.com/BioWilko/align_trim/blob/master/environment.yml). You can solve these dependencies using the minimal conda environment we have provided:

```sh
conda env create -f environment.yml
conda activate align_trim
```

#### 3. installing the pipeline:

```sh
python setup.py install
```

#### 4. testing the pipeline:

First check the pipeline can be called.

```
align_trim --version
```