# Agent-based model of the competition for IL-2 between conventional and regulatory T cells

## Project description

This project contains the Python code of an agent-based model for the competititon for IL-2 between two T cell populations. In particular, this code encompasses the one-population agent-based models described in the chapter 5 of my thesis: two-attributes deterministic dynmics, hybrid dynamics with death only, hybrid dynamics with death and activation, hybrid dynamics with death, division and activation. The extensions, such as the starvation process and the competition between two populations, are also implemented. This code does not distinguish the two cell cohorts used in the model analysis in chapter 5. 

## Folder content
The code is composed of 5 files:
- `globals.py`          which defines two lists as global variables, representing the two cell populations
- `params.py`           in which all the parameters of the model are defined
- `cells.py`            in which the two cell classes are defined and initialised
- `functions.py`        in which all the rules of the model are implemented
- `main.py`             which is the main routine: it calls the other files, initialises the model, runs the routine and saves the outputs

A detailed description on how the code works is provided in the appendix of the thesis.

## How to install and run the project

1. Download the folder and unzip it; or in command line, type 
> git clone https://github.com/leasta/ABM_thesis
2. Choose the parameter values by opening params.py, in particular choose whether the code should create an animation or run faster (boolean `animate`).
3. Run main.py with Python 3.
> Python 3 main.py

The outputs are created in folder *_Figure_DATE, where * must be replaced by the number of instances the code was run today (0 if it is the first time) and DATE is today's date in format DDMMYYYY. For instance, if one runs the code for the first time on February 1st 2023, the outputs will be found in the folder 0_Figure_01022023. 

## How to use this code
This code is, first, provided for illustration purposes, as supporting material of chapter 5 of my thesis. When choosing the same parameter values as in chapter 5, one can use this code to reproduce most figures of the chapter (without cohort distinction). A curious modeller can use the code to explore different parameter regimes or different models (with different rules) than those discussed in my thesis. This code can also be used for any further work on the competition for IL-2 

## Credits
This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 764698.


