## Permissions
For owner only: check reading permissions on input folder and make parent directory traversable:
```
chmod o+x ~
```
Root directory needs to be made accessible by IT. If just for data sharing: share from storage node.

## Download from Zenodo (if soft linking does not work)
```
wget https://zenodo.org/records/14258052/files/BEN_CI16_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_CI16_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/LGS1_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/LGS1_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_CI18_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_CI18_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW10_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW10_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW12_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW12_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW13_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_NW13_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI18_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI18_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI9_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI9_sub_2.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI19_sub_1.fq.gz
wget https://zenodo.org/records/14258052/files/BEN_SI19_sub_2.fq.gz
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
