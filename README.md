# RMT: Robust Simulations

This repository contains the data and simulation scripts for Chapter 6 in the thesis Robust Statistical Methods in the Credibility Movement of Psychological Science. For supplementary materials graphically summarising the results from all the combinations of simulated conditions, please visit [robustsim.netlify.app](https://robustsim.netlify.app/) . The associated OSF repository can be accessed here: <https://osf.io/jvdrp/>.

The repository is organised as follows:

-   `designs_between` - contains all the between-subjects designs with categorical predictors

-   `designs_numeric` - contains all cross-sectional designs with numeric only predictors

-   `designs_mixed` - contains cross-sectional designs with combined categorical and numeric predictors

-   `designs_within` - contains within-subjects and mixed designs

-   `simulation_summaries` - contains scripts for generating plots in the manuscript and on [robustsim.netlify.app](https://robustsim.netlify.app/)

Because of the wide range of simulated conditions, different designs were run in separate projects. This means that for example `designs_between/simulations_CAT-2B` is itself an R project. I recommend that the overarching R project (`osf_robust_simulations`) is not used for anything other than storing and working with github (otherwise you're working with a project within a project set-up which, as we know, can tear holes into the fabric of the universe).

If you want to re-produce the summaries but *not* the individual simulations, you need to retain the relative structure of all the folders in this repository. The `simulation_summaries` project grabs relevant files from all the relevant directories. This means that if you want to re-create summaries for between-subjects designs, it needs to be able to navigate into `designs_between`.

## Working with individual simulation designs 

Each top level directory starting with `designs_` contains R projects for specific designs. For example `designs_between/simulations_CAT-2B` is a project running the simulation for between-subjects categorical designs with 2 levels. The full list of designs is below.

| Design label  | Description                                                                                                                                                                                       |
|-------------------|-----------------------------------------------------|
| CAT-2B        | 1 categorical predictor, 2 groups, between-subjects.                                                                                                                                              |
| CAT-3B        | 1 categorical predictor, 3 groups, between-subjects.                                                                                                                                              |
| CAT-2B-2B     | 2 categorical predictors, 2 groups by 2 groups, between-subjects, no interaction.                                                                                                                 |
| CAT-2B-2B-INT | 2 categorical predictors, 2 groups by 2 groups, between-subjects, with interaction.                                                                                                               |
| CAT-2B-3B     | 2 categorical predictors, 2 groups by 3 groups, between-subjects, no interaction.                                                                                                                 |
| CAT-2B-3B-INT | 2 categorical predictors, 2 groups by 3 groups, between-subjects, with interaction.                                                                                                               |
| CAT-3B-3B     | 2 categorical predictors, 3 groups by 3 groups, between-subjects, no interaction.                                                                                                                 |
| CAT-3B-3B-INT | 2 categorical predictors, 3 groups by 3 groups, between-subjects, with interaction.                                                                                                               |
| N-NUM         | N numeric predictors (N being either 1, 2, 3 or 6), cross-sectional designs                                                                                                                       |
| 2-NUM-INT     | 2 numeric predictors, cross-sectional designs, with 2-way interaction.                                                                                                                            |
| 3-NUM-INT     | 3 numeric predictors, cross-sectional designs, with three 2-way interactions between the combinations of the predictors.                                                                          |
| N-NUM-CAT     | N-numeric predictors (N-being either 1, 3, or 8) plus additional 1 or 2 categorical predictors, cross-sectional designs, no interactions.                                                         |
| N-NUM-CAT-INT | N-numeric predictors (N-being either 1, 3, or 8) plus additional 1 or 2 categorical predictors, cross-sectional designs, with 1 or 2 interactions between a numeric and a categorical predictor.  |
| CAT-2W        | 1 categorical predictor, 2 groups, within-subjects.                                                                                                                                               |
| CAT-3W        | 1 categorical predictor, 3 groups, within-subjects.                                                                                                                                               |
| CAT-2W-2W     | 2 categorical predictors, 2 (within) x 2 (within), within-subjects, full factorial.                                                                                                               |
| CAT-2W-2W-2W  | 3 categorical predictors, 2 (within) x 2 (within) x 2 (within), within-subjects, full factorial.                                                                                                  |
| CAT-2W-2B     | 2 categorical predictors, 2 (within) x 2 (between), mixed design, full factorial.                                                                                                                 |
| CAT-2W-2W-2B  | 3 categorical predictors, 2 (within) x 2 (within) x 2 (between), mixed-design, full factorial.                                                                                                    |
| CAT-2W-2B-2B  | 3 categorical predictors, 2 (within) x 2 (between) x 2 (between), mixed-design, full factorial.                                                                                                   |
| CAT-2W-1NUM   | 1 categorical predictor with 2 groups, 1 numeric predictor, within-subjects. Not included in the final simulation.                                                                                |

: List of designs included in the simulation

For each designs, the project had the following structure:

-   `data`

    -   `processed_data` - contains the processed simulation results saved in an `.rds` file, including results for false positives, statistical power, accuracy, efficiency and convergence errors.

    -   `simulation_exports` - this contains partial exports generated during each loop during a simulation. `simulations_CAT-2B` contains an example of what this directory is supposed to look like. However, for most designs, these files were too large to upload to github. They are hosted on The University Of Sussex's Box instead and can be accessed here: <https://sussex.box.com/s/wgsemizl7ulv4zt53zxi0tyr38n6myeb> . So instead of `.rds` files, most repositories contain a text file sending you to this directory. You'll need to locate the relevant directory on Box (folde structure is identical to the one here), and download the export files into your local directory if you want to use these files.

    -   `simulation_grids` contains the settings for skewness, kurtosis, heteroscedasticity, etc., which each design reads in order to run the simulation and generate different combinations of conditions. Each design contains all of the grids, but typically only works with one. E.g. CAT-2B works with the grid `CAT-2B.csv` (the reason why there are all grids in the folder is that I needed to create a lot of copies of the projects for different designs and it was easier to leave them all in as the files are quite small anyway)

-   `r_docs`

    -   `data_processing_*.qmd` is the file that generates the `.rds` file in `data/processed_data` with all the results. The processing file grabs all the relevant files from the `data/simulation_exports` folder to do this.

    -   `full_simulation_code.qmd` is a physical manifestation of misplaced hope and ambition. It contains the simulation code that *can* , technically speaking, run all of the designs one by one. However, after running a few of the designs, it became clear that it would take several months to complete this way, so the code was re-written in script (located in a `scripts` folder) that runs each design separately, and additionally incorporates parallel processing (the bottom line is, you shouldn't need to use this file unless you have a computer that can quantum tunnel through time and space, but if you just want to see what the code looks like across designs, this is the easiest way to access it)

-   `scripts`

    -   `helpers.R` is a collection of helper functions used across all the designs, including helpers for generating datasets or fitting the quantile smoother for quantifying heteroscedasticity.

    -   `parallel_setup_script.R` should be the same across designs. It loads up some settings for the simulation and prepares some additional helper functions.

    -   `parallel_[design-name].R` is unique to each design and is the script used to actually run the simulation. It loads up the previous scripts, prepares design-specific functions and runs the simulation. Each script uses parallel processing - at the moment it's set up as `multisession` because the simulations were run on Posit Cloud, but the way the script is written should also work with the `multicore` setting.

This is the general structure used for most designs. For some designs, the computations were taking too long so the design itself was split into multiple child projects. For example, each child project would handle one effect size. In such scenario, there would be 4 child projects associated with a single design (for effect sizes 0, 0.1, 0.3 and 0.5). Whenever this happened, there is an additional `_merge` folder. This is used only for processing the data - it collects the exports from across the folders `simulation_exports` for a given design. Additionally, there is sometimes a `_parent` folder. This would have a handful of exports which were generated first time the design was set off (and we realised it would take too long to run, so we split it into child projects). The take-away for these types of designs is:

-   use individual child projects to run the simulations

-   use `_merge` project to generate summaries

-   use `_parent` projects only to access existing data (or you can use these projects to run specific designs if you have a fast computer in which case you don't need to split the design into child projects)

Designs that are affected by this:

-   3-way mixed designs in `designs_within`

-   Numeric-only designs (N-NUM) without interactions in `designs_numeric`

-   All combined predictor designs in `designs_mixed` . Generally, these designs take an extremely long time to run so they're chunked up into many smaller sub-projects.

## Working with simulation summaries 

The easiest way to access the summaries for all designs is through the website: [robustsim.netlify.app](https://robustsim.netlify.app/) . If you want to re-run the summaries yourself, the project `simulation_summaries` collates the results files from across the designs to create summary plots. The structure is the following:

-   `r_docs` - all files within this folder contain a `.qmd` version and a `.html` file, which is a rendered version of the Quarto file. Most summary files also contain plots stored in identically names folder ending with `_files` . This is generated automatically by the `.qmd` file. Resources are not embedded - for some designs, embedding the resources causes quarto to run out of memory because the generated file is so massive.

    -   `results_overview.qmd` - this is just a working note file with observed trends in the data that I compiled when going over all the summary files. This was a baseline used for writing up the results.

    -   `simulation_summaries_[design-name].qmd` renders a summary file for each set of designs.

    -   `simulation_summaries.qmd` - was an attempt to bind all the results in a single document. As mentioned above, the document was too large, so I ended up splitting it by design. Use `simulation_summaries_[design-name].qmd` instead.

-   `scripts`

    -   `helpers.R` - helper functions for creating plots

    -   `generate_[design-name].R` - these files collate all the relevant summary files and loop over the designs to generate code for summary plots for each. For each group of designs, this is written into `summary_qmds/[design-name].qmd` which is then read as a child document by `r_docs/simulation_summaries_[design-name].qmd` .
