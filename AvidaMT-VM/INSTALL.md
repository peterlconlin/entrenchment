# AvidaMT Singularity Vagrant Box Recipe

This folder contains all necessary files to build a [Singularity](https://sylabs.io/) [Vagrant Box](https://developer.hashicorp.com/vagrant/docs/boxes) capable of running [AvidaMT](https://github.com/kgskocelas/AvidaMT), a digital evolution platform for the study of the major transition to multicellularity and entrenchment in it. 

To do this, the container installs [Boost v1.71.0](https://boost.org), [Clang v9.0.0](https://clang.llvm.org/), and [ealib](https://github.com/dknoester/ealib) (the library AvidaMT is build on) in an Ubuntu 16.04 container. 


## Why Use a Singularity Vagrant Box?

AvidaMT requires a number of legacy dependencies and configuration with administrator privileges. Full-scale experiments with it are best run on a high-performance computing cluster, where most users are not administrators. As such, the easiest way to install and use it is via a virtual machine built on your local computer and then transferred to the computing cluster of your choice. 


## How to Make Your AvidaMT Singularity Vagrant Box

### Step 0: Create or sign in to your GitHub account.

### Step 1: Get and configure your copy of AvidaMT

Get a copy of AvidaMT by going to [the AvidaMT repository on GitHub](https://github.com/kgskocelas/AvidaMT.git) and clicking the "Fork" button in the upper righthand corner.

On your fork, you can find the configuration file at AvidaMT/etc/major_transitions.cfg. Edit it as desired with the parameters for your experiment.

Click "commit changes."

### Step 2: Get and configure your copy of the AvidaMT Singularity Recipe

1. Get a copy of the AvidaMT Singularity Recipe by clicking the "Fork" button in the upper righthand corner of this current repository (the one that you are reading this document in).

2. On your fork, open AvidaMT-VM/entrenchment_recipe. On line 63, replace "kgskocelas" in the hyperlink with your GitHub username so that it directs to your AvidaMT software fork from Step 1.

3. Click "commit changes."


### Step 3: Set up your Singularity Vagrant Box

Singularity runs natively on Linux and can be run on Mac and Windows through virtual machines (VMs). While the directions below are for Mac, they are largely copied from [Singularity's documentation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) which includes directions for other operating systems.

The following directions were validated on February 25, 2025 on macOS 13.7.4

1. Install virtualbox, vagrant, and vagrant-manager with the package manager of your choice via terminal. This example uses [Homebrew](https://brew.sh/).
```bash
$ brew install virtualbox --cask && \
    brew install vagrant --cask && \
    brew install vagrant-manager --cask
```

2. Create and enter a directory to be used with your Vagrant VM. Since we're making a Singularity Vagrant Box, I call mine vm-singularity.
```bash
$ mkdir vm-singularity && \
    cd vm-singularity
```

If you've already used this folder for another VM, you will need to destroy the VM and delete the Vagrantfile.
```bash
$ vagrant destroy && \
    rm Vagrantfile
```

3. Bring up a Virtual Machine.
```bash
$ export VM=sylabs/singularity-3.0-ubuntu-bionic64 && \
    vagrant init $VM && \
    vagrant up && \
    vagrant ssh
```

4. Check the installed version of Singularity with the following:
```bash
vagrant@vagrant:~$ singularity version
3.0.3-1
```

5. In your VM, initialize a Git repository.
```bash
vagrant@vagrant:~$ git init
```

6. Clone your AvidaMT Singularity Recipe from step 2 into your VM. Be sure to edit the GitHub link to match your fork.
```bash
vagrant@vagrant:~$ git clone https://github.com/YOUR_USERNAME/entrenchment.git
```

7. Build your SIMG file. This is a slow process and is best started and left to run instead of trying to use your computer while it's building.
```bash
vagrant@vagrant:~$ sudo singularity build avidaMT.simg entrenchment_recipe
```

8. Move your SIMG file to the high performance computing cluster you're using. Be sure to edit the server location with your information.
```bash
vagrant@vagrant:~$ scp avidaMT.simg YOUR_USERNAME@YOUR_SERVER_ADDRESS:
```


## Using AvidaMT

To run AvidaMT, navigate to the folder in which you want the data to appear and execute the SIMG file with:  
```bash
$ ./avidaMT.simg
```

To override a value on the command line:
```bash
$ ./avidaMT.simg --full.option.name=new_value
```
For example, if you wanted to override per-site mutation rates to be 1%, you might type:
```bash
$ ./avidaMT.simg --ea.mutation.site.p=0.01
```

To continue an existing run:
```bash
./avidaMT.simg -l <path_to_checkpoint_file>/checkpoint-1000000.xml.gz
```

To perform further analyses, load the checkpoint file and then run an analysis tool. This example also loads a line of descent file which can be analyzed.
```bash
./avidaMT.simg -l <path_to_checkpoint_file>/checkpoint-1000000.xml.gz --analyze lod_fitness --ea.analysis.input.filename <path_to_checkpoint_file>/lod-1000000.xml.gz
```

