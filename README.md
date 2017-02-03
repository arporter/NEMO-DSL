# NEMO-DSL

Investigation into the use of Domain Specific Languages with the NEMO ocean model

## Obtaining the Code ##

The benchmark codes contained in this project use the
[dl_timer](https://bitbucket.org/apeg/dl_timer) library for execution
timing. This code is included within the NEMO-DSL repository as a git submodule.
As such, when cloning the repository you need to use the `--recursive` flag in
order to get the source code of that submodule, e.g.:

    git clone --recursive https://github.com/arporter/NEMO-DSL

If you forget to do this then your cloned repository will contain an empty
`NEMO-DSL/src/shared/dl_timer` directory. To populate it you can do:

    cd NEMO-DSL
    git submodule init
    git submodule update
