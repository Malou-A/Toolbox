# Toolbox
This is a project I did during a project week for the course Binp16 - Bioinformatics: Programming in Python

The project week allowed us to choose any project and I choose one where I could explore the possibilities to make programs with a graphical user interface, using the module **tkinter**.

## Description
I have created a graphical user interface toolbox for DNA sequences. The program provides 5 different buttons, with different funtions.

**1. Permute** -
This button will take a sequence and return all different permutations possible. One can then choose any of sequences and use as a primer in the primer function.
  
**2. Primer** -
This button will allow you to enter a primer, and search a sequence for that primer.
  
**3. GC-content** -
This button will calculate the GC-content of a sequence and plot it.
  
**4. Blast** -
This button will take a nucleotide sequence and run it against the BLAST database
  
**5. Translate** -
This button will take a nucleotide sequence and translate it to the corresponding aminoacid sequence in all frames and directions
  

## Configuration

In the file **Toolbox.py** The variable "HOME_DIR" needs to be changed to the absolute path of where the file is located.

Install miniconda to create a conda environment:
### Windows
Download and install Miniconda3:

https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe

### Linux
Download

https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Install:
```bash
bash ~/Downloads/Miniconda3-latest-Linux-x86_64.sh
```
### Create environment

Create a conda environment


Using python version 3.8.3, these additional packages needs to be installed:
```bash
$ conda create --name toolbox-env python=3.8.3
$ conda activate toolbox-env
$ conda install matplotlib
$ pip install Bio
```

The toolbox will open when the script is executed.

<p align="center">
<img src="images/Toolbox_img.png" align="center">
</p>

