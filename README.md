# DeepHit Model
R implementation of python [DeepHit competing risk model](https://github.com/cmclean5/PublicHealthModels/issues/1)

---

### Run DeepHit model on example dataset in R

Run DeepHit model over our Breast Cancer anthracyline exposure cohort:

```bash
git clone -b DeepHit https://github.com/cmclean5/PublicHealthModels.git
```

```bash
cd PublicHealthModels
```

Start R, and from within R run:

```R
   source('analysisAnthracylineExposure.R')
```

---

### Install R, python, tensorflow on Mac

- **R      version 4.2 (or higher), and**
- **Python version 3.9 (or higher).**

<details>

<summary>Seting up Tensorflow with R on a Mac ARM64-M1 chip </summary>

### .bashrc file

environment variables to define in .bashrc file to set-up R (4.2), python (3.9)   

```bash
   ##TERMINAL SETTINGS and Terminal Aliases
   ##change to Unix HOME area
   export HOME=/Users/cmclean

   ## set JAVA_HOME for ARM64-based M1
   export JAVA_HOME="/Library/Java/JavaVirtualMachines/zulu-17.jdk/Contents/Home"

   ## requied environment variables for,
   ## brew/minforge installation of tensorflow
   export GRPC_PYTHON_BUILD_SYSTEM_OPENSSL=1
   export GRPC_PYTHON_BUILD_SYSTEM_ZLIB=1

   ## GSL environment variables setup
   export GSL_HOME="/opt/homebrew/Cellar/gsl/2.7.1"
   export GSL_CFLAGS="${GSL_HOME}/include"
   export GSL_LIBS="${GSL_HOME}/lib"
   export GSL_CONFIG="${GSL_HOME}/bin/gsl-config"
   export PATH="${GSL_HOME}/bin:${PATH}"
   export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GSL_HOME}/lib"

   ## openblas environment variable setup
   ## For compilers to find openblas you may need to set:
   export BLAS_HOME="/opt/homebrew/opt/openblas"
   export BLAS_CFLAGS="${BLAS_HOME}/include"
   export BLAS_LIBS="${BLAS_HOME}/lib"
   export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BLAS_HOME}/lib"

   ## lapack environment variable setup
   ## For compilers to find lapack you may need to set:
   export LAPACK_HOME="/opt/homebrew/opt/lapack"
   export LAPACK_CFLAGS="${LAPACK_HOME}/include"
   export LAPACK_LIBS="${LAPACK_HOME}/lib"
   export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${LAPACK_HOME}/lib"

   ## what is my python set-up?
   ## 1) I have ARM64-based M1 chip python 3.9 setup, using brew
   ## 2) I have ARM64-based M1 chip python 3.9 setup, using miniforge

   #export my_py_setup="brew_9" ## for brew      python 3.9 installation
   export my_py_setup="mini_9" ## for miniforge python 3.9 installation

   if [ "$my_py_setup" == "brew_9" ] || [ "$my_py_setup" == "brew_10" ]; then

       ##-----------------------------------
       ## macs ARM64-based M1 chip python setup using brew
       ##-----------------------------------
    
       ## set the python version to 3.9
       if [ "$my_py_setup" == "brew_9" ]; then
           #export py_ver=python@3.9        
           export py_ver=3.9
       fi

       ## set the python version to 3.10
       if [ "$my_py_setup" == "brew_10" ]; then
           #export py_ver=python@3.10
           export py_ver=3.10        
       fi
        
       ## export paths
       PATH="/opt/homebrew/opt/python@${py_ver}/bin:$PATH"         ## path to python
       PATH="/Users/cmclean/Library/Python/${py_ver}/bin:$PATH"    ## path to pip3 and virtualenv
       export PATH

       ##-----------------------------------
       ## macs ARM64-based M1 chip python 3.9 setup, using brew
       ##-----------------------------------
       export PYTHONPATH="$PYTHONPATH:/opt/homebrew/opt/python@${py_ver}/bin/python3"
       export PYTHONSTARTUP=".pythonstartup.py"
       export WORKON_HOME=$HOME/.virtualenvs
       export PROJECT_HOME=$HOME/projects
       export VIRTUALENVWRAPPER_PYTHON=`which python3`
       export VIRTUALENVWRAPPER_VIRTUALENV=`which virtualenv`
       source `which virtualenvwrapper.sh`
       export LDFLAGS="-L/opt/homebrew/opt/${py_ver}/lib"
       export PKG_CONFIG_PATH="/opt/homebrew/opt/${py_ver}/lib/pkgconfig"

    ## install tensorflow for apple M1
    ## python3 -m pip install --upgrade tensorflow-macos
    ## work around for matlibplot lib
    function frameworkpython {
        if [[ ! -z "$VIRTUALENV" ]]; then
            PYTHONHOME=$VIRTUALENV $VIRTUALENVWRAPPER_PYTHON "$@"
        else
            $VIRTUALENVWRAPPER_PYTHON "$@"
        fi
    }

    echo ": python $py_ver setup using brew"
    
fi

if [ "$my_py_setup" == "mini_9" ]; then

    ##-----------------------------------
    ## macs ARM64-based M1 chip python setup using miniforge
    ##-----------------------------------
    
    ## requied environment variables for,
    ## minforge installation of tensorflow
    #export GRPC_PYTHON_BUILD_SYSTEM_OPENSSL=1
    #export GRPC_PYTHON_BUILD_SYSTEM_ZLIB=1

    ## init bash shell for miniforge
    source ${HOME}/SCRIPTS/UTILITIES/conda_init.sh

    ## start our shell without conda base environment activated
    conda deactivate
    
    ##-----------------------------------
    ## which python to use for reticulate in R
    ##-----------------------------------
    export RETICULATE_PYTHON="/opt/homebrew/Caskroom/miniforge/base/envs/r-reticulate/bin/python3"

    echo ": python 3.9 setup using miniforge"
    
fi
    
##-----------------------------------

##-----------------------------------
## R setup
##-----------------------------------
## Keep an eye for updates at: https://mac.r-project.org/
## Refs:
## [1] https://mac.r-project.org/tools/
## [2] https://www.r-bloggers.com/2021/02/fully-native-m1-apple-silicon-r-setup/
## [3] https://colinfay.me/r-installation-administration/installing-r-under-macos.html
## [4] https://cran.r-project.org/bin/macosx/

#export my_r_setup="x84_64" ## R 4.1 intel    x84_64 setup
export my_r_setup="arm_64" ## R 4.2 apple M1 arm_64 setup       

if [ "$my_r_setup" == "x84_64" ]; then

    ## gfortran for R (intel)
    PATH="/usr/local/gfortran/bin:$PATH"
    PATH="/usr/local/tcl-tk/8.6.12/bin:$PATH"
    export TCLTK_LIBS="/usr/local/tcl-tk/8.6.12/lib"
    export TCLTK_CPPFLAGS="/usr/local/tcl-tk/8.6.12/include"
    PATH="/Library/Frameworks/R.framework/Versions/4.1/Resources/bin:$PATH"
    export PATH

    ## location to where R installs packages
    export R_LIBS_USER="$HOME/.R/R-4.1/library"

    echo ": R 4.1 setup for mac intel x84_64"

    ##-----------------------------------------------------
    ## Make use of this new BLAS library:
    ## [5] https://pat-s.me/transitioning-from-x86-to-arm64-on-macos-experiences-of-an-r-user/
    ## 1) cd /Library/Frameworks/R.framework/Resources/lib/
    ## create a symbolic link pointing libRblas.dylib to the optimized BLAS implementation
    ## 2) ln -s -i -v libRblas.vecLib.dylib libRblas.dylib
    ## If you ever want to revert this, do
    ## 1) cd /Library/Frameworks/R.framework/Resources/lib/
    ## 2) ln -s -i -v libRblas.0.dylib libRblas.dylib
    ##-----------------------------------------------------
    
fi

if [ "$my_r_setup" == "arm_64" ]; then    

    ## gfortran for R (arm64)
    ## https://mac.r-project.org/
    ## wget https://mac.r-project.org/monterey/R-devel/arm64/R-devel.tar.gz
    ## tar fvxz R*.tar.gz -C /
    PATH="/opt/R/arm64/gfortran/bin:$PATH"
    PATH="/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/bin:$PATH"
    export TCLTK_LIBS="/opt/R/arm64/tcl-tk/8.6.12/lib"
    export TCLTK_CPPFLAGS="/opt/R/arm64/tcl-tk/8.6.12/include"
    export PATH

    ## location to where R installs packages
    export R_LIBS_USER="$HOME/.R/R-4.2.0/library"
```
</details>


<details>

<summary>Seting up R on a mac ARM64-M1 chip </summary>

### Installing R version 4.2 (or higher)

The easiest way is to install R on Mac is through [CRAN](https://cran.r-project.org) by going to the CRAN downloads page and following the links For Apple silicon (M1/M2) Macs. The next step is to click on the "R-4.3.1-arm64.pkg" (or newer version) file to begin the installation.

### Installing RStudio

To download RStudio, go to the [RStudio downloads page](https://posit.co/download/rstudio-desktop/#download) and get the .dmg for Mac OS, remember to keep default installation options.

### R packages to install for tensorflow

Start R from bash or RStudio

```R
install.packages("base64enc")
install.packages("reticulate")
install.packages("keras")
install.packages("tensorflow")
```

### Note install tensorflow and keras in R, might have to download & install manually
```bash
cd ~/Downloads
wget https://www.stats.bris.ac.uk/R/bin/macosx/big-sur-arm64/contrib/4.2/tensorflow_2.7.0.tgz
R CMD INSTALL tensorflow_2.7.0.tgz
wget https://www.stats.bris.ac.uk/R/bin/macosx/big-sur-arm64/contrib/4.2/keras_2.7.0.tgz
R CMD INSTALL keras_2.7.0.tgz
```

</details>



<details>

<summary>Installing Tensorflow on a Mac ARM64-M1 chip </summary>

### install miniforge via brew

```bash   
arch -arm64 brew install miniforge
```

### (1) create a new environment called r-reticulate
```bash
 conda create --name r-reticulate python=3.9
```

You'll find this environment create in 
```bash
/opt/homebrew/Caskroom/miniforge/base/envs/r-reticulate
```

### (2) to activate conda, first init it
```bash
source conda_init.sh
```
Where the `conda_init.sh` script is:
```bash
#!/bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/homebrew/Caskroom/miniforge/base/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh" ]; then
        . "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
    else
        export PATH="/opt/homebrew/Caskroom/miniforge/base/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

### (3) now can activate our environment
```bash
conda activate r-reticulate
```

### (4) install packages need for plotting keras model
```bash
 pip install pydot
 arch -arm64 brew install graphviz
```

### (5) install tensorflow
```bash
 conda install -c apple tensorflow-deps
 python -m pip install tensorflow-macos
 python -m pip install tensorflow-metal
 python -m pip install tensorflow-addons
```
### (6) You'll find tensorflow installed at
```bash
/opt/homebrew/Caskroom/miniforge/base/pkgs/
```

</details>

---

### Task-list

- [x] translate [DeepHit](https://github.com/cmclean5/PublicHealthModels/issues/1) from python to R
- [ ] translate [Dynamic-DeepHit](https://github.com/cmclean5/PublicHealthModels/issues/3) from python to R
- [ ] hyperparameter tuning of models
