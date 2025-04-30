## Permissions
For owner only: check reading permissions on input folder and make parent directory traversable:
```
chmod o+x ~
```

## Conda issues
```
conda config --set ssl_verify false
```

## Install Python3.9
```
conda install python=3.9
```

## Create an environment for multiqc which uses python 3.9 and activate it
```
conda create -n multiqc
```
```
conda activate multiqc
```

## Download and install MultiQC
```
git clone https://github.com/MultiQC/MultiQC.git
```
```
cd MultiQC
```
```
pip install .
```
