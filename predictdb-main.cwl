{
    "class": "Workflow",
    "cwlVersion": "v1.2",
    "id": "dave/predixscan/predictdb-main/1",
    "label": "Predictdb_Split_by_Chromsome",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "inputs": [
        {
            "id": "Omics_Data_File",
            "type": "File",
            "label": "Omics Data File",
            "doc": "A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized.",
            "sbg:x": -711.9014892578125,
            "sbg:y": 191.67303466796875
        },
        {
            "id": "genotype_file",
            "type": "File",
            "doc": "Input file is expected to be a tab-delimited text file including a\nheader row, where the header field is ID (for snp varID) and then a\nvariable number of fields with the sample ID numbers.  The first column\nhas the snp ID in the format of (chr_pos_refAll_effAll_build) and the\ndosages are encoded on a 0-2 scale representing the number or imputed\nnumber of the effect alleles the sample possesses.",
            "sbg:x": -367.69915771484375,
            "sbg:y": -222.9619140625
        },
        {
            "id": "snp_annotation",
            "type": "File",
            "doc": "Input snp annotation file is expected to be a tab-delimited text file,\nwith a header row, with fields for chromosome, position, variantID,\nreference allele, alternative allele, rsid_label1, rsid_label2, and \nnumber of alternative alleles per site.",
            "sbg:x": -230.09579467773438,
            "sbg:y": -412.2391052246094
        },
        {
            "id": "gene_annotation",
            "type": "File",
            "sbg:x": 33.34193420410156,
            "sbg:y": -73.43011474609375
        },
        {
            "id": "Hidden_Factors",
            "type": "int",
            "label": "NUM_FACTORS",
            "doc": "Number of hidden factors to estimate. PEER uses automatic relevance determination to choose a suitable effective number of factors, so this parameter needs only to be set to a sufficiently large value. Without prior information available, a general recommendation is to use 25% of the number of samples but no more than 100 factors.",
            "sbg:exposed": true
        },
        {
            "id": "Output_File_Prefix",
            "type": "string",
            "label": "Output File Prefrix",
            "doc": "File name prefix for output files.",
            "sbg:x": -727.713134765625,
            "sbg:y": -129.49014282226562
        },
        {
            "id": "Covariates",
            "type": "File?",
            "label": "Covariate_Data",
            "doc": "--Covariates_File\n\n[REQUIRED]: A tab delimited file with a .tab file extension containing a matrix of size M + 1 × K + 1, where M >= N and is the number of samples for which covariate data is provided. .",
            "sbg:x": -352.8101806640625,
            "sbg:y": 297.5672912597656
        },
        {
            "id": "tissue_name",
            "type": "string",
            "sbg:x": 346.3424987792969,
            "sbg:y": -463.78082275390625
        },
        {
            "id": "population_name",
            "type": "string",
            "sbg:x": 297.26922607421875,
            "sbg:y": 177
        },
        {
            "id": "chr_list",
            "type": "string[]",
            "doc": "Enter the chromosome number for which you want to perform the nested elastic models for every gene in that desired chromosome. Please enter at least 1 chromosome and each chromosome will be in its own line.",
            "sbg:exposed": true
        },
        {
            "id": "model_prefix",
            "type": "string?",
            "doc": "This prefix will be added to the beginning of all the summary, covariance, and weights files produced.",
            "sbg:exposed": true
        },
        {
            "id": "seed",
            "type": "int?",
            "sbg:exposed": true
        },
        {
            "id": "RAM",
            "type": "int?",
            "doc": "In  MB",
            "sbg:exposed": true
        }
    ],
    "outputs": [
        {
            "id": "weights",
            "outputSource": [
                "model_bulding/weights"
            ],
            "type": "File[]?",
            "sbg:x": 720.01806640625,
            "sbg:y": -386.60723876953125
        },
        {
            "id": "summary",
            "outputSource": [
                "model_bulding/summary"
            ],
            "type": "File[]?",
            "sbg:x": 771.1414184570312,
            "sbg:y": -1.123328447341919
        },
        {
            "id": "covariances",
            "outputSource": [
                "model_bulding/covariances"
            ],
            "type": "File[]?",
            "sbg:x": 762.0408325195312,
            "sbg:y": 193.58444213867188
        },
        {
            "id": "db_output",
            "outputSource": [
                "database_summary/db_output"
            ],
            "type": "File[]?",
            "sbg:x": 1228.90966796875,
            "sbg:y": -252.17791748046875
        }
    ],
    "steps": [
        {
            "id": "peer",
            "in": [
                {
                    "id": "Omics_Data_File",
                    "source": "Omics_Data_File"
                },
                {
                    "id": "Output_File_Prefix",
                    "source": "Output_File_Prefix"
                },
                {
                    "id": "Hidden_Factors",
                    "source": "Hidden_Factors"
                }
            ],
            "out": [
                {
                    "id": "std_out"
                },
                {
                    "id": "peer_covariates"
                },
                {
                    "id": "peer_weights"
                },
                {
                    "id": "peer_precisions"
                },
                {
                    "id": "peer_residuals"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.1",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "dave/peer-development/peer/27",
                "baseCommand": [
                    "Rscript run_new_peer.R"
                ],
                "inputs": [
                    {
                        "id": "Omics_Data_File",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--omic_data_file",
                            "shellQuote": false,
                            "position": 0
                        },
                        "label": "Omics Data File",
                        "doc": "A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized."
                    },
                    {
                        "id": "Output_File_Prefix",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--output_prefix",
                            "shellQuote": false,
                            "position": 0
                        },
                        "label": "Output File Prefrix",
                        "doc": "File name prefix for output files."
                    },
                    {
                        "id": "Output_Directory",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--output_dir",
                            "shellQuote": false,
                            "position": 0
                        },
                        "label": "Output Directory",
                        "doc": "To specify an output directory"
                    },
                    {
                        "id": "Hidden_Factors",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--num_factors",
                            "shellQuote": false,
                            "position": 0
                        },
                        "label": "NUM_FACTORS",
                        "doc": "Number of hidden factors to estimate. PEER uses automatic relevance determination to choose a suitable effective number of factors, so this parameter needs only to be set to a sufficiently large value. Without prior information available, a general recommendation is to use 25% of the number of samples but no more than 100 factors."
                    },
                    {
                        "id": "Covariates_File",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--cov_file",
                            "shellQuote": false,
                            "position": 0
                        },
                        "label": "Covariates File",
                        "doc": "A tab delimited file with a .tab file extension containing a matrix of size M + 1 × C + 1, where M >= N and is the number of samples for which covariate data is provided. If this file is input, the set of samples used in the hidden factor estimation procedure will be the intersection of samples in the covariate matrix and omic data matrix. C is the number of known covariates to be included in association test regression models of downstream analyses. Examples of common covariates include sex, age, batch variables, and quality metrics. Categorical variables (e.g., batch number) have to be encoded as D - 1 indicator/binary variables, where D is the number of categories for a given categorical variable. For the indicator variables, a value of 1 signifies membership in the category and a value of 0 indicates otherwise. The first row and column must contain sample IDs and covariate IDs respectively. \n[Default = NULL]"
                    },
                    {
                        "id": "Alpha_Prior_A",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--alphaprior_a",
                            "shellQuote": false,
                            "position": 0
                        },
                        "doc": "Shape parameter of the gamma distribution prior of the model noise distribution [Default=0.001].",
                        "default": 0.001
                    },
                    {
                        "id": "Alpha_Prior_B",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--alphaprior_b",
                            "shellQuote": false,
                            "position": 0
                        },
                        "doc": "Scale parameter of the gamma distribution prior of the model noise distribution. [default=0.01]",
                        "default": 0.01
                    },
                    {
                        "id": "Eps_Prior_A",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--epsprior_a",
                            "shellQuote": false,
                            "position": 0
                        },
                        "doc": "Shape parameter of the gamma distribution prior of the model weight distribution. [Default=0.1]",
                        "default": 0.1
                    },
                    {
                        "id": "Eps_Prior_B",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--epsprior_b",
                            "shellQuote": false,
                            "position": 0
                        },
                        "doc": "Scale parameter of the gamma distribution prior of the model weight distribution. [default=10]",
                        "default": 10
                    },
                    {
                        "id": "Tol",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--tol",
                            "shellQuote": false,
                            "position": 0
                        },
                        "doc": "Threshold for the increase in model evidence when optimizing hidden factor values. Estimation completes for a hidden factor when the increase in model evidence exceeds this value [default=0.001].",
                        "default": 0.001
                    },
                    {
                        "id": "Max_Iteration",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--max_iter",
                            "shellQuote": false,
                            "position": 0
                        },
                        "doc": "Max number of iterations for updating values of each hidden factor [default=1000].",
                        "default": 1000
                    }
                ],
                "outputs": [
                    {
                        "id": "std_out",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "std.out"
                        }
                    },
                    {
                        "id": "peer_covariates",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*peer_covariates.txt"
                        }
                    },
                    {
                        "id": "peer_weights",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*peer_weights.txt"
                        }
                    },
                    {
                        "id": "peer_precisions",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*peer_precisions.txt"
                        }
                    },
                    {
                        "id": "peer_residuals",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*peer_residuals.txt"
                        }
                    }
                ],
                "doc": "```\nUsage: \nRscript run_peer.R [options] --omic_data file <omic_data_file> --output_prefix <output_prefix> --num_factors <num_factors>\n\nProbabilistic Estimation of Expression Residuals (PEER)\n\nRun PEER using the R interface. PEER is a method designed to estimate surrogate variables/latent factors/hidden factors that contribute to gene expression variability, but it can be applied to other data types as well. For more information, please refer to https://doi.org/10.1038/nprot.2011.457.\n\nOptions:\n        -h, --help\n                Show this help message and exit\n\n        --omic_data_file=OMIC_DATA_FILE\n                [REQUIRED] A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized.\n\n        --output_prefix=OUTPUT_PREFIX\n                [REQUIRED] File name prefix for output files. To specify an output directory as well, use --output_dir.\n\n        --num_factors=NUM_FACTORS\n                [REQUIRED] Number of hidden factors to estimate. PEER uses automatic relevance determination to choose a suitable effective number of factors, so this parameter needs only to be set to a sufficiently large value. Without prior information available, a general recommendation is to use 25% of the number of samples but no more than 100 factors.\n\n        --cov_file=COV_FILE\n                A tab delimited file with a .tab file extension containing a matrix of size M + 1 × C + 1, where M >= N and is the number of samples for which covariate data is provided. If this file is input, the set of samples used in the hidden factor estimation procedure will be the intersection of samples in the covariate matrix and omic data matrix. C is the number of known covariates to be included in association test regression models of downstream analyses. Examples of common covariates include sex, age, batch variables, and quality metrics. Categorical variables (e.g., batch number) have to be encoded as D - 1 indicator/binary variables, where D is the number of categories for a given categorical variable. For the indicator variables, a value of 1 signifies membership in the category and a value of 0 indicates otherwise. The first row and column must contain sample IDs and covariate IDs respectively [default=NULL].\n\n        --alphaprior_a=ALPHAPRIOR_A\n                Shape parameter of the gamma distribution prior of the model noise distribution [default=0.001].\n\nroot@80bc4504de2b:/opt# Rscript run_peer.R --help\nUsage: \nRscript run_peer.R [options] --omic_data file <omic_data_file> --output_prefix <output_prefix> --num_factors <num_factors>\n\nProbabilistic Estimation of Expression Residuals (PEER)\n\nRun PEER using the R interface. PEER is a method designed to estimate surrogate variables/latent factors/hidden factors that contribute to gene expression variability, but it can be applied to other data types as well. For more information, please refer to https://doi.org/10.1038/nprot.2011.457.\n\nOptions:\n        -h, --help\n                Show this help message and exit\n\n        --omic_data_file=OMIC_DATA_FILE\n                [REQUIRED] A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized.\n\n        --output_prefix=OUTPUT_PREFIX\n                [REQUIRED] File name prefix for output files. To specify an output directory as well, use --output_dir.\n\n        --num_factors=NUM_FACTORS\n                [REQUIRED] Number of hidden factors to estimate. PEER uses automatic relevance determination to choose a suitable effective number of factors, so this parameter needs only to be set to a sufficiently large value. Without prior information available, a general recommendation is to use 25% of the number of samples but no more than 100 factors.\n\n        --cov_file=COV_FILE\n                A tab delimited file with a .tab file extension containing a matrix of size M + 1 × C + 1, where M >= N and is the number of samples for which covariate data is provided. If this file is input, the set of samples used in the hidden factor estimation procedure will be the intersection of samples in the covariate matrix and omic data matrix. C is the number of known covariates to be included in association test regression models of downstream analyses. Examples of common covariates include sex, age, batch variables, and quality metrics. Categorical variables (e.g., batch number) have to be encoded as D - 1 indicator/binary variables, where D is the number of categories for a given categorical variable. For the indicator variables, a value of 1 signifies membership in the category and a value of 0 indicates otherwise. The first row and column must contain sample IDs and covariate IDs respectively [default=NULL].\n\n        --alphaprior_a=ALPHAPRIOR_A\n                Shape parameter of the gamma distribution prior of the model noise distribution [default=0.001].\n\n        --alphaprior_b=ALPHAPRIOR_B\n                Scale parameter of the gamma distribution prior of the model noise distribution. [default=0.01]\n\n        --epsprior_a=EPSPRIOR_A\n                Shape parameter of the gamma distribution prior of the model weight distribution. [default=0.1]\n\n        --epsprior_b=EPSPRIOR_B\n                Scale parameter of the gamma distribution prior of the model weight distribution. [default=10]\n\n        --tol=TOL\n                Threshold for the increase in model evidence when optimizing hidden factor values. Estimation completes for a hidden factor when the increase in model evidence exceeds this value [default=0.001].\n\n        --var_tol=VAR_TOL\n                Threshold for the variance of model residuals when optimizing hidden factor values. Estimation completes for a hidden factor when the variance of residuals is smaller than this value [default=1e-05].\n\n        --max_iter=MAX_ITER\n                Max number of iterations for updating values of each hidden factor [default=1000].\n\n        -o OUTPUT_DIR, --output_dir=OUTPUT_DIR\n                Directory in which to save outputs [default=.].\n\n        -v, --version\n                Print PEER version number.\n\n\n\n\n``",
                "label": "peer",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "davidroberson/peer:20.10.14"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "run_new_peer.R",
                                "entry": "# Modified from https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R\n# Original Author: Francois Aguet\n# Modified by: Bryan Quach <bquach@rti.org>\n\nlibrary(peer, quietly = T)\nlibrary(optparse, quietly = T)\npeer.version <- 1.3 #software version\n\nWriteTable <- function(data, filename, index.name) {\n    datafile <- file(filename, open = \"wt\")\n    on.exit(close(datafile))\n    header <- c(index.name, colnames(data))\n    writeLines(paste0(header, collapse = \"\\t\"), con = datafile, sep = \"\\n\")\n    write.table(data, datafile, sep = \"\\t\", col.names = F, quote = F)\n}\n\n# Generate usage doc and retrieve command line args\np <- OptionParser(usage = \"\\n%prog [options] --omic_data file <omic_data_file> --output_prefix <output_prefix> --num_factors <num_factors>\",\n    description = \"\\nProbabilistic Estimation of Expression Residuals (PEER)\\n\\nRun PEER using the R interface. PEER is a method designed to estimate surrogate variables/latent factors/hidden factors that contribute to gene expression variability, but it can be applied to other data types as well. For more information, please refer to https://doi.org/10.1038/nprot.2011.457.\",\n    prog = \"Rscript run_peer.R\")\np <- add_option(object = p, opt_str = c(\"--omic_data_file\"), default = NULL, type = \"character\",\n    help = \"[REQUIRED] A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized.\")\np <- add_option(object = p, opt_str = c(\"--output_prefix\"), default = NULL, type = \"character\",\n    help = \"[REQUIRED] File name prefix for output files. To specify an output directory as well, use --output_dir.\")\np <- add_option(object = p, opt_str = c(\"--num_factors\"), type = \"integer\", default = NULL,\n    help = \"[REQUIRED] Number of hidden factors to estimate. PEER uses automatic relevance determination to choose a suitable effective number of factors, so this parameter needs only to be set to a sufficiently large value. Without prior information available, a general recommendation is to use 25% of the number of samples but no more than 100 factors.\")\np <- add_option(object = p, opt_str = c(\"--cov_file\"), help = \"A tab delimited file with a .tab file extension containing a matrix of size M + 1 × C + 1, where M >= N and is the number of samples for which covariate data is provided. If this file is input, the set of samples used in the hidden factor estimation procedure will be the intersection of samples in the covariate matrix and omic data matrix. C is the number of known covariates to be included in association test regression models of downstream analyses. Examples of common covariates include sex, age, batch variables, and quality metrics. Categorical variables (e.g., batch number) have to be encoded as D - 1 indicator/binary variables, where D is the number of categories for a given categorical variable. For the indicator variables, a value of 1 signifies membership in the category and a value of 0 indicates otherwise. The first row and column must contain sample IDs and covariate IDs respectively [default=%default].\")\np <- add_option(object = p, opt_str = c(\"--alphaprior_a\"), type = \"double\", default = 0.001,\n    help = \"Shape parameter of the gamma distribution prior of the model noise distribution [default=%default].\")\np <- add_option(object = p, opt_str = c(\"--alphaprior_b\"), type = \"double\", default = 0.01,\n    help = \"Scale parameter of the gamma distribution prior of the model noise distribution. [default=%default]\")\np <- add_option(object = p, opt_str = c(\"--epsprior_a\"), type = \"double\", default = 0.1,\n    help = \"Shape parameter of the gamma distribution prior of the model weight distribution. [default=%default]\")\np <- add_option(object = p, opt_str = c(\"--epsprior_b\"), type = \"double\", default = 10,\n    help = \"Scale parameter of the gamma distribution prior of the model weight distribution. [default=%default]\")\np <- add_option(object = p, opt_str = c(\"--tol\"), type = \"double\", default = 0.001,\n    help = \"Threshold for the increase in model evidence when optimizing hidden factor values. Estimation completes for a hidden factor when the increase in model evidence exceeds this value [default=%default].\")\np <- add_option(object = p, opt_str = c(\"--var_tol\"), type = \"double\", default = 0.00001,\n    help = \"Threshold for the variance of model residuals when optimizing hidden factor values. Estimation completes for a hidden factor when the variance of residuals is smaller than this value [default=%default].\")\np <- add_option(object = p, opt_str = c(\"--max_iter\"), type = \"double\", default = 1000,\n    help = \"Max number of iterations for updating values of each hidden factor [default=%default].\")\np <- add_option(object = p, opt_str = c(\"--output_dir\", \"-o\"), default = \".\",\n    help = \"Directory in which to save outputs [default=%default].\")\np <- add_option(object = p, opt_str = c(\"--version\", \"-v\"), action = \"store_true\", default = F, \n    help = \"Print PEER version number.\")\nargv <- parse_args(p)\n\n# Quick execution for printing version number\nif(argv$version){\n    cat(paste0(\"PEER v\", peer.version))\n    quit(save = \"no\")\n}\n\n# Check if positional arguments were given \nif(is.null(argv$omic_data_file)){\n    stop(\"Error: Please provide a value for --omic_data_file\")\n}\nif(is.null(argv$output_prefix)){\n    stop(\"Error: Please provide a value for --output_prefix\")\n}\nif(is.null(argv$num_factors)){\n    stop(\"Error: Please provide a value for --num_factors\")\n}\n\n\n# Check validity of argument inputs\nif(!file.exists(argv$omic_data_file)){ \n    stop(paste0(\"Error: \", argv$omic_data_file, \n        \" not found. Check your file path and name.\")) \n}\nif(!grepl(x = argv$omic_data_file, pattern = \"\\\\.tab$\", perl = T)){ \n    stop(\"Error: --omic_data_file requires .tab extension.\")\n}\nif(!is.null(argv$cov_file) && !file.exists(argv$cov_file)){ \n    stop(paste0(\"Error: \", argv$cov_file, \n        \" not found. Check your file path and name.\")) \n}\nif(!is.null(argv$cov_file) &&\n    !grepl(x = argv$cov_file, pattern = \"\\\\.tab$\", perl = T)){ \n    stop(\"Error: --cov_file requires .tab extension.\")\n}\nif(is.na(argv$num_factors) | argv$num_factors <= 0 | \n   !is.finite(argv$num_factors) | argv$num_factors != as.integer(argv$num_factors)){ \n    stop(paste0(\"Error: Please provide a valid number for 'num_factors'. Use --help for more details.\"))\n}\nif(argv$max_iter <= 0 | !is.finite(argv$max_iter) | \n   argv$max_iter != as.integer(argv$max_iter)){ \n    stop(paste0(\"Error: Please provide a valid number for --max_iter. Use --help for more details.\"))\n}\nif(argv$tol <= 0 | !is.finite(argv$tol)){ \n    stop(paste0(\"Error: Please provide a valid value for --tol. Use --help for more details.\"))\n}\nif(argv$var_tol <= 0 | !is.finite(argv$var_tol)){ \n    stop(paste0(\"Error: Please provide a valid value for --var_tol. Use --help for more details.\"))\n}\nif(argv$alphaprior_a < 0 | !is.finite(argv$alphaprior_a)){ \n    stop(paste0(\"Error: Please provide a valid value for --alphaprior_a. Use --help for more details.\"))\n}\nif(argv$alphaprior_b < 0 | !is.finite(argv$alphaprior_b)){ \n    stop(paste0(\"Error: Please provide a valid value for --alphaprior_b. Use --help for more details.\"))\n}\nif(argv$epsprior_a < 0 | !is.finite(argv$epsprior_a)){ \n    stop(paste0(\"Error: Please provide a valid value for --epsprior_a. Use --help for more details.\"))\n}\nif(argv$epsprior_b < 0 | !is.finite(argv$epsprior_b)){ \n    stop(paste0(\"Error: Please provide a valid value for --epsprior_b. Use --help for more details.\"))\n}\n\n# Create output directory if needed\ndir.create(argv$output_dir, showWarnings = F)\n\n# Load omic data\ncat(paste0(\"Loading data from \", argv$omic_data_file, \" ...\"))\nomic.data <- read.table(argv$omic_data_file, sep = \"\\t\", header = T, \n    check.names = F, comment.char = \"\", row.names = 1)\nomic.data <- as.matrix(omic.data)\nn.samples <- nrow(omic.data)\nn.features <- ncol(omic.data)\ncat(\"Done.\\n\")\ncat(paste0(\"Loaded data matrix with \", n.samples, \" rows and \", \n    n.features, \" columns.\\n\"))\n\n# Load covariate data\ncov.data <- NULL\nif(!is.null(argv$cov_file)){\n    cat(paste0(\"Loading covariate data from \", argv$cov_file, \" ...\"))\n    cov.data <- read.table(argv$cov_file, sep = \"\\t\", header = T, as.is = T, \n        check.names = F, comment.char = \"\", row.names = 1)\n    cov.data <- as.matrix(cov.data)\n    n.vars <- ncol(cov.data)\n    cat(\"Done.\\n\")\n    cat(paste0(\"Loaded \", n.vars, \" covariates.\\n\"))\n    # Subset and match rows between covariate and omic data matrix\n    cov.subset <- cov.data[rownames(cov.data) %in% rownames(omic.data), , drop = F]\n    omic.subset <- omic.data[rownames(omic.data) %in% rownames(cov.subset), , drop = F]\n    match.order <- match(rownames(cov.subset), table = rownames(omic.subset))\n    cov.subset <- cov.subset[match.order, , drop = F]\n    if(nrow(omic.subset) < nrow(omic.data)){\n        cat(paste0(\"Data reduced to \", nrow(omic.subset), \n        \" samples after matching with covariate data.\\n\"))\n    }\n    omic.data <- omic.subset\n    cov.data <- cov.subset\n}\n\n# Set method parameters\ncat(paste0(\"Setting initialization parameters ...\"))\nmodel <- PEER()\ninvisible(PEER_setNk(model, argv$num_factors))\ninvisible(PEER_setPhenoMean(model, omic.data))\ninvisible(PEER_setPriorAlpha(model, argv$alphaprior_a, argv$alphaprior_b))\ninvisible(PEER_setPriorEps(model, argv$epsprior_a, argv$epsprior_b))\ninvisible(PEER_setTolerance(model, argv$tol))\ninvisible(PEER_setVarTolerance(model, argv$var_tol))\ninvisible(PEER_setNmax_iterations(model, argv$max_iter))\nif(!is.null(cov.data)){\n    invisible(PEER_setCovariates(model, cov.data))\n}\ncat(\"Done.\\n\")\n\n# Run inference routine\ncat(paste0(\"Beginning estimation of \", argv$num_factors, \" hidden factors.\\n\"))\ntime <- system.time(PEER_update(model))\ncat(\"Finished estimation procedure.\\n\")\ncat(paste0(\"Hidden factor estimation completed in \", round(time[3], 2) , \" seconds (\", round(time[3]/60, 2) ,\" minutes).\\n\"))\n\n# Retrieve results\nfactor.mat <- PEER_getX(model)  # samples x PEER factors\nweight.mat <- PEER_getW(model)  # omic features x PEER factors\nprecision.mat <- PEER_getAlpha(model)  # PEER factors x 1\nresid.mat <- t(PEER_getResiduals(model))  # omic features x samples\n\n# Add relevant row/column names\npeer.var.names <- paste0(\"peer.factor\", 1:ncol(factor.mat))\nrownames(factor.mat) <- rownames(omic.data)\ncolnames(factor.mat) <- peer.var.names\ncolnames(weight.mat) <- peer.var.names\nrownames(weight.mat) <- colnames(omic.data)\nrownames(precision.mat) <- peer.var.names\ncolnames(precision.mat) <- \"alpha\"\nprecision.mat <- as.data.frame(precision.mat)\nprecision.mat$relevance <- 1.0 / precision.mat$alpha\nrownames(resid.mat) <- colnames(omic.data)\ncolnames(resid.mat) <- rownames(omic.data)\n\n# Write results\ncat(\"Exporting results ... \")\nWriteTable(factor.mat, file.path(argv$output_dir, paste0(argv$output_prefix, \"_peer_covariates.txt\")), \"ID\")  \nWriteTable(weight.mat, file.path(argv$output_dir, paste0(argv$output_prefix, \"_peer_weights.txt\")), \"ID\")\nWriteTable(precision.mat, file.path(argv$output_dir, paste0(argv$output_prefix, \"_peer_precisions.txt\")), \"ID\")\nWriteTable(resid.mat, file.path(argv$output_dir, paste0(argv$output_prefix, \"_peer_residuals.txt\")), \"ID\")\ncat(\"Done.\\n\")",
                                "writable": false
                            }
                        ]
                    }
                ],
                "stdout": "std.out",
                "sbg:appVersion": [
                    "v1.1"
                ],
                "sbg:content_hash": "a25c6d4875eb1d525dd7f3a89371fb389d55b0971bfdf57ee9951b729caf9df63",
                "sbg:contributors": [
                    "rk.johnson",
                    "e.esquinca",
                    "dave"
                ],
                "sbg:createdBy": "dave",
                "sbg:createdOn": 1602689199,
                "sbg:id": "dave/peer-development/peer/27",
                "sbg:image_url": null,
                "sbg:latestRevision": 27,
                "sbg:modifiedBy": "e.esquinca",
                "sbg:modifiedOn": 1607574135,
                "sbg:project": "dave/peer-development",
                "sbg:projectName": "PEER Development",
                "sbg:publisher": "sbg",
                "sbg:revision": 27,
                "sbg:revisionNotes": "Transpose Peer Factor Results",
                "sbg:revisionsInfo": [
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602689199,
                        "sbg:revision": 0,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602689235,
                        "sbg:revision": 1,
                        "sbg:revisionNotes": "added help text"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602689438,
                        "sbg:revision": 2,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602689466,
                        "sbg:revision": 3,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602689579,
                        "sbg:revision": 4,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602689985,
                        "sbg:revision": 5,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602690267,
                        "sbg:revision": 6,
                        "sbg:revisionNotes": "std out"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602690736,
                        "sbg:revision": 7,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602690866,
                        "sbg:revision": 8,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602690914,
                        "sbg:revision": 9,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1602691010,
                        "sbg:revision": 10,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603125039,
                        "sbg:revision": 11,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603125227,
                        "sbg:revision": 12,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603125661,
                        "sbg:revision": 13,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603127020,
                        "sbg:revision": 14,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603128327,
                        "sbg:revision": 15,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603128415,
                        "sbg:revision": 16,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603911381,
                        "sbg:revision": 17,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603911497,
                        "sbg:revision": 18,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603911583,
                        "sbg:revision": 19,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603923045,
                        "sbg:revision": 20,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1603924351,
                        "sbg:revision": 21,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1604007308,
                        "sbg:revision": 22,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "rk.johnson",
                        "sbg:modifiedOn": 1604509981,
                        "sbg:revision": 23,
                        "sbg:revisionNotes": "Added output port for PEER covariates"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1604535891,
                        "sbg:revision": 24,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1604535931,
                        "sbg:revision": 25,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1604535981,
                        "sbg:revision": 26,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1607574135,
                        "sbg:revision": 27,
                        "sbg:revisionNotes": "Transpose Peer Factor Results"
                    }
                ],
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": []
            },
            "label": "peer",
            "sbg:x": -363.38470458984375,
            "sbg:y": 0.5600862503051758
        },
        {
            "id": "split_snp_annot_by_chromosome",
            "in": [
                {
                    "id": "snp_annotation",
                    "source": "snp_annotation"
                },
                {
                    "id": "output_prefrix",
                    "default": "snp_annot"
                }
            ],
            "out": [
                {
                    "id": "snp_annot_split_by_chr_files"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.1",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "rk.johnson/predictdb/split-snp-annot-by-chromosome/12",
                "baseCommand": [
                    "python split_snp_annot_by_chr.py"
                ],
                "inputs": [
                    {
                        "id": "snp_annotation",
                        "type": "File?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 1
                        },
                        "doc": "Input snp annotation file is expected to be a tab-delimited text file,\nwith a header row, with fields for chromosome, position, variantID,\nreference allele, alternative allele, rsid_label1, rsid_label2, and \nnumber of alternative alleles per site."
                    },
                    {
                        "id": "output_prefrix",
                        "type": "string?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 2
                        },
                        "doc": "The suffix 'chrN.txt' will be added\nto the prefix provided, where N is the chromosome number"
                    }
                ],
                "outputs": [
                    {
                        "id": "snp_annot_split_by_chr_files",
                        "doc": "The output files will be tab-delimited text files with chromosome,\nposition, variantID, reference allele, effect allele, and rsid.\nNOTE: the rsid number chosen is from rsidlabel2.",
                        "type": "File[]",
                        "outputBinding": {
                            "glob": "*.txt"
                        }
                    }
                ],
                "doc": "Script to split a SNP annotation file into multiple files by chromosome.\nFrom commandline, first argument is the snp annotation file, second is\nthe prefix for the output files.  The suffix 'chrN.txt' will be added\nto the prefix provided, where N is the chromosome number.\nIn splitting, script will only keep unambiguously stranded SNPs. I.e.,\nno INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G\nand vice-versa.\n\nInput snp annotation file is expected to be a tab-delimited text file,\nwith a header row, with fields for chromosome, position, variantID,\nreference allele, alternative allele, rsid_label1, rsid_label2, and \nnumber of alternative alleles per site. See file\nGTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt.gz\nfrom gtexportal.org for an example of such a file.\n\nThe output files will be tab-delimited text files with chromosome,\nposition, variantID, reference allele, effect allele, and rsid.\nNOTE: the rsid number chosen is from rsidlabel2.",
                "label": "Split SNP Annot By Chromosome",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "python:3.10.0a4-slim"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "split_snp_annot_by_chr.py",
                                "entry": "#! /usr/bin/env python3\n\nimport sys\n\n'''\nScript to split a SNP annotation file into multiple files by chromosome.\nFrom commandline, first argument is the snp annotation file, second is\nthe prefix for the output files.  The suffix 'chrN.txt' will be added\nto the prefix provided, where N is the chromosome number.\nIn splitting, script will only keep unambiguously stranded SNPs. I.e.,\nno INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G\nand vice-versa.\nInput snp annotation file is expected to be a tab-delimited text file,\nwith a header row, with fields for chromosome, position, variantID,\nreference allele, alternative allele, rsid_label1, rsid_label2, and \nnumber of alternative alleles per site. See file\nGTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt.gz\nfrom gtexportal.org for an example of such a file.\nThe output files will be tab-delimited text files with chromosome,\nposition, variantID, reference allele, effect allele, and rsid.\nNOTE: the rsid number chosen is from rsidlabel2.\n'''\n\nSNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}\nHEADER_FIELDS = ['chr','pos','varID','ref_vcf','alt_vcf','rsid']\n\ndef split_snp_annot(annot_file, out_prefix):\n    # Make output file names from prefix.\n    snps_by_chr_files= [out_prefix + '.chr' + str(i) + '.txt' for i in range(1,23)]\n    # Open connection to each output file\n    snp_by_chr = [open(f, 'w') for f in snps_by_chr_files]\n    # Write header in each file.\n    header = '\\t'.join(HEADER_FIELDS)+'\\n'\n    for f in snp_by_chr:\n        f.write(header)\n    with open(annot_file, 'r') as ann:\n        # Skip header from input file\n        ann.readline()\n        # Extract rows from input and write to body in appropriate output.\n        for line in ann:\n            attrs = line.split()\n            chr = attrs[0]\n            pos = attrs[1]\n            varID = attrs[2]\n            refAllele = attrs[3]\n            effectAllele = attrs[4]\n            rsid = attrs[5]\n            # Skip non-single letter polymorphisms\n            if len(refAllele) > 1 or len(effectAllele) > 1:\n                continue\n            # Skip ambiguous strands\n            if SNP_COMPLEMENT[refAllele] == effectAllele:\n                continue\n            if rsid == '.':\n                continue\n            index = int(chr) - 1\n            row = '\\t'.join([chr,pos,varID,refAllele,effectAllele,rsid])+'\\n'\n            snp_by_chr[index].write(row)\n    # Close connection to each output file.\n    for f in snp_by_chr:\n        f.close()\n\nif __name__ == '__main__':\n    annot_file = sys.argv[1]\n    out_prefix = sys.argv[2]\n    split_snp_annot(annot_file, out_prefix)",
                                "writable": false
                            }
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.py"
                    }
                ],
                "sbg:appVersion": [
                    "v1.1"
                ],
                "sbg:content_hash": "af4ad9b477a5714209b27c4c225acbecd68a56758d640c0f02fa22ff60cdb544b",
                "sbg:contributors": [
                    "e.esquinca"
                ],
                "sbg:createdBy": "e.esquinca",
                "sbg:createdOn": 1610154582,
                "sbg:id": "rk.johnson/predictdb/split-snp-annot-by-chromosome/12",
                "sbg:image_url": null,
                "sbg:latestRevision": 12,
                "sbg:modifiedBy": "e.esquinca",
                "sbg:modifiedOn": 1612209538,
                "sbg:project": "rk.johnson/predictdb",
                "sbg:projectName": "predictdb",
                "sbg:publisher": "sbg",
                "sbg:revision": 12,
                "sbg:revisionNotes": "changed back to 5",
                "sbg:revisionsInfo": [
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610154582,
                        "sbg:revision": 0,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610155165,
                        "sbg:revision": 1,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610999467,
                        "sbg:revision": 2,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611004492,
                        "sbg:revision": 3,
                        "sbg:revisionNotes": "update for our data"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611341601,
                        "sbg:revision": 4,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611344004,
                        "sbg:revision": 5,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611601137,
                        "sbg:revision": 6,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611601532,
                        "sbg:revision": 7,
                        "sbg:revisionNotes": "match names"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611851797,
                        "sbg:revision": 8,
                        "sbg:revisionNotes": "saved 5 for tutorial data run"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611852471,
                        "sbg:revision": 9,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612206342,
                        "sbg:revision": 10,
                        "sbg:revisionNotes": "rsID -> rsid"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612206826,
                        "sbg:revision": 11,
                        "sbg:revisionNotes": "needed to rerun tutorial data changed from 5 to 6"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612209538,
                        "sbg:revision": 12,
                        "sbg:revisionNotes": "changed back to 5"
                    }
                ],
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": []
            },
            "label": "Split SNP Annot By Chromosome",
            "sbg:x": 73.9994888305664,
            "sbg:y": -381.50750732421875
        },
        {
            "id": "split_chromosome",
            "in": [
                {
                    "id": "genotype_file",
                    "source": "genotype_file"
                },
                {
                    "id": "out_prefix",
                    "default": "genotype"
                }
            ],
            "out": [
                {
                    "id": "geno_split_by_chr_files"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.1",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "rk.johnson/predictdb/split-chromosome/14",
                "baseCommand": [
                    "python split_genotype_by_chr.py"
                ],
                "inputs": [
                    {
                        "id": "genotype_file",
                        "type": "File?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 1
                        },
                        "doc": "Input file is expected to be a tab-delimited text file including a\nheader row, where the header field is ID (for snp varID) and then a\nvariable number of fields with the sample ID numbers.  The first column\nhas the snp ID in the format of (chr_pos_refAll_effAll_build) and the\ndosages are encoded on a 0-2 scale representing the number or imputed\nnumber of the effect alleles the sample possesses."
                    },
                    {
                        "id": "out_prefix",
                        "type": "string?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 2
                        },
                        "doc": "String appended to out put files. The suffix 'chrN.txt' will be added\nto the prefix provided, where N is the chromosome number."
                    }
                ],
                "outputs": [
                    {
                        "id": "geno_split_by_chr_files",
                        "type": "File[]",
                        "outputBinding": {
                            "glob": "*.txt",
                            "outputEval": "$(inheritMetadata(self, inputs.text_file_input))"
                        }
                    }
                ],
                "doc": "Script to split a GTEx genotype file into multiple files by chromosome.\nFrom commandline, first argument is the genotype file, second is the\nprefix for the output files.  The suffix 'chrN.txt' will be added to the\nprefix provided, where N is the chromosome number.\n\nIn splitting, script will only keep unambiguously stranded SNPs. I.e.,\nno INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G\nand vice-versa.\n\nThe input file is expected to be a tab-delimited text file including a\nheader row, where the header field is ID (for snp varID) and then a\nvariable number of fields with the sample ID numbers.  The first column\nhas the snp ID in the format of (chr_pos_refAll_effAll_build) and the\ndosages are encoded on a 0-2 scale representing the number or imputed\nnumber of the effect alleles the sample possesses.",
                "label": "Split Genotype by Chromosome",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "python:3.10.0a4-slim"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "split_genotype_by_chr.py",
                                "entry": "#! /usr/bin/env python3\n\nimport os\nimport sys\n\n'''\nScript to split a GTEx genotype file into multiple files by chromosome.\nFrom commandline, first argument is the genotype file, second is the\nprefix for the output files.  The suffix 'chrN.txt' will be added to the\nprefix provided, where N is the chromosome number.\nIn splitting, script will only keep unambiguously stranded SNPs. I.e.,\nno INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G\nand vice-versa.\nThe input file is expected to be a tab-delimited text file including a\nheader row, where the header field is ID (for snp varID) and then a\nvariable number of fields with the sample ID numbers.  The first column\nhas the snp ID in the format of (chr_pos_refAll_effAll_build) and the\ndosages are encoded on a 0-2 scale representing the number or imputed\nnumber of the effect alleles the sample posseses.\n'''\n\nSNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}\n\ndef split_genotype(geno_file, out_prefix):\n    # Make output file names from prefix.\n    geno_by_chr_fns = [out_prefix + '.chr' + str(i) + '.txt' for i in range(1,23)]\n    # Open connection to each output file.\n    geno_by_chr = [open(f, 'w') for f in geno_by_chr_fns]\n\n    with open(geno_file, 'r') as geno:\n        # Write header in each file\n        header = geno.readline()\n        snps = set()\n        for f in geno_by_chr:\n            f.write(header)\n\n        for line in geno:\n            # First attribute of line is is chr_pos_refAllele_effAllele_build\n            # Extract this attribute and parse into list\n            varID_list = (line.split()[0].split('_'))\n            chr = varID_list[0]\n            refAllele = varID_list[2]\n            effectAllele = varID_list[3]\n            # Skip non_single letter polymorphisms\n            if len(refAllele) > 1 or len(effectAllele) > 1:\n                continue\n            # Skip ambiguous strands\n            if SNP_COMPLEMENT[refAllele] == effectAllele:\n                continue\n            varID = '_'.join(varID_list)\n            # Some snps have 2 rows for some reason. Attributes are nearly\n            # identical. Only keep the first one found.\n            if varID in snps:\n                continue\n            snps.add(varID)\n            # Write line to appropriate file\n            index = int(chr) - 1\n            geno_by_chr[index].write(line)\n\n    for f in geno_by_chr:\n        f.close()\n\nif __name__ == '__main__':\n    genotype_file = sys.argv[1]\n    out_prefix = sys.argv[2]\n    split_genotype(genotype_file, out_prefix)",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.py"
                    }
                ],
                "sbg:appVersion": [
                    "v1.1"
                ],
                "sbg:content_hash": "a4028e06034916fcea0a56bde1b41e03949812bb51c9718a4d3bb2e2c93b298b1",
                "sbg:contributors": [
                    "e.esquinca",
                    "dave"
                ],
                "sbg:createdBy": "e.esquinca",
                "sbg:createdOn": 1608313025,
                "sbg:id": "rk.johnson/predictdb/split-chromosome/14",
                "sbg:image_url": null,
                "sbg:latestRevision": 14,
                "sbg:modifiedBy": "e.esquinca",
                "sbg:modifiedOn": 1611000955,
                "sbg:project": "rk.johnson/predictdb",
                "sbg:projectName": "predictdb",
                "sbg:publisher": "sbg",
                "sbg:revision": 14,
                "sbg:revisionNotes": "",
                "sbg:revisionsInfo": [
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1608313025,
                        "sbg:revision": 0,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1608313253,
                        "sbg:revision": 1,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610122612,
                        "sbg:revision": 2,
                        "sbg:revisionNotes": "added input output ports"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610122648,
                        "sbg:revision": 3,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610123685,
                        "sbg:revision": 4,
                        "sbg:revisionNotes": "added sed"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610123960,
                        "sbg:revision": 5,
                        "sbg:revisionNotes": "sed 's/Id/varID/g'"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610124281,
                        "sbg:revision": 6,
                        "sbg:revisionNotes": "> added"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610124422,
                        "sbg:revision": 7,
                        "sbg:revisionNotes": "python:3.10.0a4-slim"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1610125209,
                        "sbg:revision": 8,
                        "sbg:revisionNotes": "removed sed"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610151487,
                        "sbg:revision": 9,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610153190,
                        "sbg:revision": 10,
                        "sbg:revisionNotes": "update name"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610154527,
                        "sbg:revision": 11,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610155359,
                        "sbg:revision": 12,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610392664,
                        "sbg:revision": 13,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1611000955,
                        "sbg:revision": 14,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": []
            },
            "label": "Split Genotype by Chromosome",
            "sbg:x": -57.45125961303711,
            "sbg:y": -223.16526794433594
        },
        {
            "id": "mlr",
            "in": [
                {
                    "id": "Peer_Covariates",
                    "source": "peer/peer_covariates"
                },
                {
                    "id": "Covariates",
                    "source": "Covariates"
                },
                {
                    "id": "Omic_Data",
                    "source": "Omics_Data_File"
                },
                {
                    "id": "Output_Prefix",
                    "source": "Output_File_Prefix"
                }
            ],
            "out": [
                {
                    "id": "Adjusted_Omics_Residuals"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.1",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "rk.johnson/predictdb/mlr/37",
                "baseCommand": [
                    "Rscript MLR.R"
                ],
                "inputs": [
                    {
                        "id": "Peer_Covariates",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--PEER_covariates",
                            "shellQuote": false,
                            "position": 3
                        },
                        "label": "PEER_Covariates",
                        "doc": "--PEER_Covariates\n(Also known as PEER Factors)\n\n[REQUIRED]: PEER - Probabilistic Estimation of Expression Residuals obtained from PEER tool. Expected to be a tab delimited file with N+1 by M+1, where N is the number of samples, and M the number of PEER factors"
                    },
                    {
                        "id": "Covariates",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--covariates_file",
                            "shellQuote": false,
                            "position": 2
                        },
                        "label": "Covariate_Data",
                        "doc": "--Covariates_File\n\n[REQUIRED]: A tab delimited file with a .tab file extension containing a matrix of size M + 1 × K + 1, where M >= N and is the number of samples for which covariate data is provided. ."
                    },
                    {
                        "id": "Omic_Data",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--omic_file",
                            "shellQuote": false,
                            "position": 1
                        },
                        "doc": "--Omic_Data_File\n\n[REQUIRED] A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.)."
                    },
                    {
                        "id": "Output_Prefix",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--output_prefix",
                            "shellQuote": false,
                            "position": 4
                        },
                        "doc": "--Output_Prefix\n\n[REQUIRED] File name prefix for output files."
                    }
                ],
                "outputs": [
                    {
                        "id": "Adjusted_Omics_Residuals",
                        "doc": "After the MLR is done on every column of the peer factors, the residuals will be store in the matrix of size m x n. The rows are the samples and the columns will be the residuals.\nThese are the adjusted omics data residuals.",
                        "type": "File",
                        "outputBinding": {
                            "glob": "*residuals.txt"
                        }
                    }
                ],
                "label": "MLR",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "r-base"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "MLR.R",
                                "entry": "# R Script used for the Multiple Linear Regression\n\n\n# Script will loop through all the columns of the transposed gene expression which correspond to each \n# gene and for each gene it runs linear regression on the PEER factors & covariates. Then\n# it sets the residuals to the new expression for that gene.\n\n\n# Install Dependencies\ninstall.packages(\"optparse\")\nlibrary(optparse, quietly = T)\n\n\n\n# Generate usage doc and retrieve command line arguments\np <- OptionParser(usage = \"\\n%prog [options] --omic_file <omic_file> --covariates_data_file <covariates_data_file> --PEER_covariates <PEER_covariates> --output_prefix <output_prefix>\",\n                  description = \"\\nScript will loop through all the columns of the transposed gene expression which correspond to each \n                        gene and for each gene it runs linear regression on the PEER factors & seperate covariates file is entered. Then\n                        it sets the residuals to the new expression for that gene\",\n                  prog = \"Rscript MLR.R\")\n\np <- add_option(object = p, opt_str = c(\"--omic_file\"), default = NULL, type = \"character\",\n                help = \"[REQUIRED] A tab delimited file with .tab file extension containing N + 1 rows and G + 1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized.\")\np <- add_option(object = p, opt_str = c(\"--covariates_file\"), default = NULL, type = \"character\",\n                help = \"A tab deliminated file with N+1 rows and K+1 columns, where N is the number of samples, and K is the desired covariates.\")\np <- add_option(object = p, opt_str = c(\"--PEER_covariates\"), type = \"character\", default = NULL,\n                help = \"[REQUIRED] PEER - Probabilistic Estimation of Expression Residuals obtained from PEER tool. Expected to be a tab deliminated file with N+1 by M+1, where N is the number of samples, and M the number of PEER factors \")\np <- add_option(object = p, opt_str = c(\"--output_prefix\"), default = NULL, type = \"character\",\n                help = \"[REQUIRED] File name prefix for output files.\")\n\n\nargv <- parse_args(p)\n\n\n# Check if positional arguments were given \nif(is.null(argv$omic_file)){\n  stop(\"Error: Please provide a value for --omic_file\")\n}\nif(is.null(argv$PEER_covariates)){\n  stop(\"Error: Please provide a value for --PEER_covariates\")\n}\nif(is.null(argv$output_prefix)){\n  stop(\"Error: Please provide a value for --output_prefix\")\n}\n\n# Read in data\ngene.exp <- read.table(argv$omic_file, sep = \"\\t\", header = T, \n                       check.names = F, comment.char = \"\", row.names = 1)\n\npeer.covariates <- read.table(argv$PEER_covariates, sep = \"\\t\", header = T, \n                           check.names = F, comment.char = \"\", row.names = 1)\n\n\ncovar.data = NULL\nif(!is.null(argv$covariates_file)){\n  covar.data <- read.table(argv$covariates_file, sep = \"\\t\", header = T, \n                           check.names = F, comment.char = \"\", row.names = 1)\n}\n\n\n# Make a copy of the gene.exp df and fill in with the residuals\nexpression <- gene.exp\n\n\n# Run MLR\n\n\n# If covariate data entered, merge peer covariates file with extra covariates file\n# Then run the MLR\n\nif(!is.null(argv$covariates_file)){\n  # Will merge by rownames\n  merged.covar <- merge(peer.covariates, covar.data, by = 0)\n  \n  # Make first column rownames and get rid of extra first column\n  rownames(merged.covar) <- merged.covar[,1]\n  merged.covar <- merged.covar[,-1]\n  \n  # Run the MLR\n  for (i in 1:length(colnames(gene.exp))) {\n    fit <- lm(gene.exp[,i] ~ as.matrix(merged.covar))\n    expression[,i] <- fit$residuals\n  }\n  \n}else{ \n\n  for (i in 1:length(colnames(gene.exp))) {\n  fit <- lm(gene.exp[,i] ~ as.matrix(peer.covariates))\n  expression[,i] <- fit$residuals\n}\n}\n \n\n\n# Write results\nwrite.table(expression, row.names = T, sep = \"\\t\", col.names = T, file = paste0(argv$output_prefix, \"_adjusted_omics_residuals.txt\"))\n         \n    \n            ",
                                "writable": false
                            }
                        ]
                    }
                ],
                "sbg:appVersion": [
                    "v1.1"
                ],
                "sbg:content_hash": "ac6dff55ea6f0a7e9332b7698514032bb6944d9b157ada845d603747d6012dade",
                "sbg:contributors": [
                    "e.esquinca"
                ],
                "sbg:createdBy": "e.esquinca",
                "sbg:createdOn": 1607551910,
                "sbg:id": "rk.johnson/predictdb/mlr/37",
                "sbg:image_url": null,
                "sbg:latestRevision": 37,
                "sbg:modifiedBy": "e.esquinca",
                "sbg:modifiedOn": 1612853548,
                "sbg:project": "rk.johnson/predictdb",
                "sbg:projectName": "predictdb",
                "sbg:publisher": "sbg",
                "sbg:revision": 37,
                "sbg:revisionNotes": "",
                "sbg:revisionsInfo": [
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1607551910,
                        "sbg:revision": 0,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1607552400,
                        "sbg:revision": 1,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1607552562,
                        "sbg:revision": 2,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1608318765,
                        "sbg:revision": 3,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1608318898,
                        "sbg:revision": 4,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610153154,
                        "sbg:revision": 5,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610653560,
                        "sbg:revision": 6,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610653637,
                        "sbg:revision": 7,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610653786,
                        "sbg:revision": 8,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610654249,
                        "sbg:revision": 9,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610654330,
                        "sbg:revision": 10,
                        "sbg:revisionNotes": "Updated Script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610655503,
                        "sbg:revision": 11,
                        "sbg:revisionNotes": "Update Script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610656726,
                        "sbg:revision": 12,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610657587,
                        "sbg:revision": 13,
                        "sbg:revisionNotes": "added out log file"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610657913,
                        "sbg:revision": 14,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610658907,
                        "sbg:revision": 15,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610663483,
                        "sbg:revision": 16,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610663915,
                        "sbg:revision": 17,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610664251,
                        "sbg:revision": 18,
                        "sbg:revisionNotes": "binding"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610751467,
                        "sbg:revision": 19,
                        "sbg:revisionNotes": "update script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610751592,
                        "sbg:revision": 20,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610752059,
                        "sbg:revision": 21,
                        "sbg:revisionNotes": "update script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610752826,
                        "sbg:revision": 22,
                        "sbg:revisionNotes": "Descriptions"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610753016,
                        "sbg:revision": 23,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610753180,
                        "sbg:revision": 24,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610993763,
                        "sbg:revision": 25,
                        "sbg:revisionNotes": "output script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610994304,
                        "sbg:revision": 26,
                        "sbg:revisionNotes": "write.table edit"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1610995238,
                        "sbg:revision": 27,
                        "sbg:revisionNotes": "output"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612832140,
                        "sbg:revision": 28,
                        "sbg:revisionNotes": "made covariates optional"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612842028,
                        "sbg:revision": 29,
                        "sbg:revisionNotes": "edited output"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612848474,
                        "sbg:revision": 30,
                        "sbg:revisionNotes": "edited names"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612849368,
                        "sbg:revision": 31,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612851486,
                        "sbg:revision": 32,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612851877,
                        "sbg:revision": 33,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612851919,
                        "sbg:revision": 34,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612851953,
                        "sbg:revision": 35,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612853089,
                        "sbg:revision": 36,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612853548,
                        "sbg:revision": 37,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": []
            },
            "label": "MLR",
            "sbg:x": -78.39913940429688,
            "sbg:y": 168.8870391845703
        },
        {
            "id": "model_bulding",
            "in": [
                {
                    "id": "population_name",
                    "source": "population_name"
                },
                {
                    "id": "tissue_name",
                    "source": "tissue_name"
                },
                {
                    "id": "gene_annotation",
                    "source": "gene_annotation"
                },
                {
                    "id": "snp_annotation",
                    "source": [
                        "split_snp_annot_by_chromosome/snp_annot_split_by_chr_files"
                    ]
                },
                {
                    "id": "genotype_file",
                    "source": [
                        "split_chromosome/geno_split_by_chr_files"
                    ]
                },
                {
                    "id": "adjusted_expression_file",
                    "source": "mlr/Adjusted_Omics_Residuals"
                },
                {
                    "id": "chr_list",
                    "source": [
                        "chr_list"
                    ]
                },
                {
                    "id": "model_prefix",
                    "source": "model_prefix"
                },
                {
                    "id": "seed",
                    "source": "seed"
                },
                {
                    "id": "RAM",
                    "source": "RAM"
                }
            ],
            "out": [
                {
                    "id": "weights"
                },
                {
                    "id": "covariances"
                },
                {
                    "id": "summary"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "rk.johnson/predictdb-development/model-creation/50",
                "baseCommand": [
                    "bash elnet.sh"
                ],
                "inputs": [
                    {
                        "id": "population_name",
                        "type": "string?",
                        "doc": "Population from which the samples originated from."
                    },
                    {
                        "id": "tissue_name",
                        "type": "string?",
                        "doc": "The tissue the expression was measured in"
                    },
                    {
                        "id": "gene_annotation",
                        "type": "File?"
                    },
                    {
                        "id": "snp_annotation",
                        "type": "File[]?"
                    },
                    {
                        "id": "genotype_file",
                        "type": "File[]?"
                    },
                    {
                        "id": "adjusted_expression_file",
                        "type": "File?",
                        "doc": "Expression was adjusted by performing a multivariate linear regression with all covariates, pulling the residual values, and then assigning the residuals to be the new expression values."
                    },
                    {
                        "id": "chr_list",
                        "type": "string[]?",
                        "doc": "Enter the chromosome number for which you want to perform the nested elastic models for every gene in that desired chromosome. Please enter at least 1 chromosome and each chromosome will be in its own line."
                    },
                    {
                        "id": "model_prefix",
                        "type": "string?",
                        "doc": "This prefix will be added to the beginning of all the summary, covariance, and weights files produced.",
                        "default": "\"Model_Training\""
                    },
                    {
                        "id": "seed",
                        "type": "int?"
                    },
                    {
                        "id": "MAF",
                        "type": "float?",
                        "default": 0.01
                    },
                    {
                        "id": "n_folds",
                        "type": "int?",
                        "default": 10
                    },
                    {
                        "id": "n_train_test_folds",
                        "type": "int?",
                        "default": 5
                    },
                    {
                        "id": "alpha",
                        "type": "float?",
                        "doc": "The alpha parameter used with the R package glmnet to train the elastic net model.",
                        "default": 0.5
                    },
                    {
                        "id": "RAM",
                        "type": "int?",
                        "doc": "In  MB"
                    }
                ],
                "outputs": [
                    {
                        "id": "weights",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "weights/*"
                        }
                    },
                    {
                        "id": "covariances",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "covariances/*"
                        }
                    },
                    {
                        "id": "summary",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "summary/*"
                        }
                    }
                ],
                "doc": "Program runs elastic net prediction models following PredictDB_Pipeline_GTEx_v7, \nas described here: https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7\n\n1. Functions defined\n2. Define inputs, directories, & call functions\n\nPart 3 is conducted in the Database Summary Tool\n\n3. Combine models for all chromosomes entered in one database file. Then the database will be filtered for use in PrediXcan\n\nNested Cross Validated Elastic-Net - In previous versions of PredictDB, we employed 10-fold cross-validated elastic-net to tune the parameter lambda, and then estimated the significance of the model. It recently became apparent that this was biasing the significance measures because we were using the same data to tune the parameter lambda and assess the performance. To correct for this problem, we used the following \"nested\" cross validation procedure:\n\nRandomly split the data into 5 folds.\n\nFor each fold:\n\na. Remove the fold from the data.\n\nb. Use the remaining data to train an elastic-net model using 10-fold cross-validation to tune the lambda parameter.\n\nc. With the trained model, predict on the hold out fold, and get various test statistics for how the model performs.\n\nCalculate the average and standard deviation of each of the significance statistics, where applicable. This should provide a reasonable estimate for how well the model will generalize to new data.\n\nTrain a new elastic-net model using all of the data. Again, use 10-fold cross validation to tune the lambda parameter. The non-zero weights from this final model are what are saved in the database, provided the model meets significance criteria.\n\nA model was determined to be \"significant\" if the average pearson correlation between predicted and observed during nested cross validation was greater than 0.1 (equivalent to R2 > 0.01) and the estimated p-value for this statistic was less than 0.05. See below for how the p-value was calculated.",
                "label": "Elastic Models",
                "requirements": [
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "$(inputs.RAM)"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sb.biodatacatalyst.nhlbi.nih.gov/dave/predictdb:v2021_02_13"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "gtex_v7_nested_cv_elnet_training_combined.R",
                                "entry": "#! /usr/bin/env Rscript\n\n# Program runs elastic net prediction models following PredictDB_Pipeline_GTEx_v7, \n# as described here: https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7\n\n# 1. Functions defined\n# 2. Define inputs, directories, call functions\n# 3. Combine models for all chromosomes, create and filter database for use in PrediXcan\n\nsuppressMessages(library(dplyr))\nsuppressMessages(library(glmnet))\nsuppressMessages((library(reshape2)))\nsuppressMessages(library(methods))\nsuppressMessages(library(RSQLite))\nsuppressMessages(library(data.table))\n\n\"%&%\" <- function(a,b) paste(a,b, sep='')\n\n\n##################################\n# 1. Define all functions for prediction models\n\nget_filtered_snp_annot <- function(snp_annot_file) {\n  snp_annot <- read.table(snp_annot_file, header = T, stringsAsFactors = F) %>%\n    filter(!((ref_vcf == 'A' & alt_vcf == 'T') |\n               (ref_vcf == 'T' & alt_vcf == 'A') |\n               (ref_vcf == 'C' & alt_vcf == 'G') |\n               (ref_vcf == 'G' & alt_vcf == 'C')) &\n             !(is.na(rsid))) %>%\n    distinct(varID, .keep_all = TRUE)\n  snp_annot\n}\n\n\nget_maf_filtered_genotype <- function(genotype_file,  maf, samples) {\n  gt_df <- read.table(genotype_file, header = T, stringsAsFactors = F, row.names = 1)\n  gt_df <- gt_df[,(colnames(gt_df) %in% samples )] %>% t() %>% as.data.frame()\n  effect_allele_freqs <- colMeans(gt_df) / 2\n  gt_df <- gt_df[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf))]\n  gt_df\n}\n\nget_gene_annotation <- function(gene_annot_file, chrom, gene_types=c('protein_coding', 'pseudogene', 'lincRNA')){\n  gene_df <- read.table(gene_annot_file, header = TRUE, stringsAsFactors = FALSE) %>%\n    filter((chr == chrom) & gene_type %in% gene_types)\n  gene_df\n}\n\nget_gene_type <- function(gene_annot, gene) {\n  filter(gene_annot, gene_id == gene)$gene_type\n}\n\n# Got rid of t() done twice which should cancel out. \nget_gene_expression <- function(expression_file, gene_annot) {\n  expr_df <- as.data.frame((read.table(expression_file, header = T, stringsAsFactors = F, row.names = 1)))\n  #expr_df <- expr_df %>% t() %>% as.data.frame()\n  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))\n  expr_df\n}\n\nget_gene_coords <- function(gene_annot, gene) {\n  row <- gene_annot[which(gene_annot$gene_id == gene),]\n  c(row$start, row$end)\n}\n\nget_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {\n  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(rsid)) & (pos <= (coords[2] + cis_window)))\n  if (nrow(snp_info) == 0)\n    return(NA)\n  #Check if the varID exist in the data\n  if (TRUE %in% (snp_info$varID %in% names(gt_df))) {\n    cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))\n  } else {\n    return(NA) # the varID doesn't exist in the gt_df dataset\n  }\n  column_labels <- colnames(cis_gt)\n  row_labels <- rownames(cis_gt)\n  # Convert cis_gt to a matrix for glmnet\n  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt))\n  colnames(cis_gt) <- column_labels\n  rownames(cis_gt) <- row_labels\n  cis_gt\n}\n\n#get_covariates <- function(covariates_file, samples) {\n#  cov_df <- read.table(covariates_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1)\n#  cov_df <- cov_df[(rownames(cov_df) %in% samples),] %>% as.data.frame() # %&% t() \n  # We have peer covariates coming out as samples as rows from the tool so deleted the extra t()\n  # and adjusted reading in rownames instead\n#  cov_df\n#}\n\ngenerate_fold_ids <- function(n_samples, n_folds=10) {\n  n <- ceiling(n_samples / n_folds)\n  fold_ids <- rep(1:n_folds, n)\n  sample(fold_ids[1:n_samples])\n}\n\n# adjust_for_covariates <- function(expression_vec, cov_df) {\n#  combined_df <- cbind(expression_vec, cov_df)\n#  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals\n#  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)\n#  expr_resid\n# }\n\ncalc_R2 <- function(y, y_pred) {\n  tss <- sum(y**2)\n  rss <- sum((y - y_pred)**2)\n  1 - rss/tss\n}\n\ncalc_corr <- function(y, y_pred) {\n  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))\n}\n\nnested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha, samples) {\n  # Gets performance estimates for k-fold cross-validated elastic-net models.\n  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,\n  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is\n  # cross validated. Then get performance measures for how the model predicts on the hold-out\n  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis\n  # is there is no correlation between prediction and observed.\n  #\n  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values\n  # are combined using Fisher's method.\n  \n  # for testing line by line only\n  # x = cis_gt\n  # y = adj_expression\n  # n_k_folds = n_folds\n  \n  R2_folds <- rep(0, n_train_test_folds)\n  corr_folds <- rep(0, n_train_test_folds)\n  zscore_folds <- rep(0, n_train_test_folds)\n  pval_folds <- rep(0, n_train_test_folds)\n  # Outer-loop split into training and test set.\n  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)\n  for (test_fold in 1:n_train_test_folds) {\n    train_idxs <- which(train_test_fold_ids != test_fold)\n    test_idxs <- which(train_test_fold_ids == test_fold)\n    x_train <- x[(rownames(x) %in% samples[train_idxs]), ]\n    y_train <- y[(rownames(y) %in% rownames(x_train))]\n    x_test <- x[(rownames(x) %in% samples[test_idxs]), ]\n    y_test <- y[(rownames(y) %in% rownames(x_test))]\n    # Inner-loop - split up training set for cross-validation to choose lambda.\n    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)\n    y_pred <- tryCatch({\n      # Fit model with training data.\n        # Parallel\n        library(doMC)\n        registerDoMC(cores = 16)\n      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids,\n        parallel = TRUE)\n      # Predict test data using model that had minimal mean-squared error in cross validation.\n      predict(fit, x_test, s = 'lambda.min')},\n      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)\n      error = function(cond) rep(mean(y_train), length(y_test)))\n    R2_folds[test_fold] <- calc_R2(y_test, y_pred)\n    # Get p-value for correlation test between predicted y and actual y.\n    # If there was no model, y_pred will have var=0, so cor.test will yield NA.\n    # In that case, give a random number from uniform distribution, which is what would\n    # usually happen under the null.\n    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)\n    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation\n    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))\n  }\n  R2_avg <- mean(R2_folds)\n  R2_sd <- sd(R2_folds)\n  rho_avg <- mean(corr_folds)\n  rho_se <- sd(corr_folds)\n  rho_avg_squared <- rho_avg**2\n  # Stouffer's method for combining z scores.\n  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)\n  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)\n  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method\n  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)\n  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)\n}\n\ndo_covariance <- function(gene_id, cis_gt, rsids, varIDs) {\n  model_gt <- cis_gt[,varIDs, drop=FALSE]\n  colnames(model_gt) <- rsids\n  geno_cov <- cov(model_gt)\n  geno_cov[lower.tri(geno_cov)] <- NA\n  cov_df <- reshape2::melt(geno_cov, varnames = c(\"rsid1\", \"rsid2\"), na.rm = TRUE) %>%\n    mutate(gene=gene_id) %>%\n    select(GENE=gene, RSID1=rsid1, RSID2=rsid2, VALUE=value) %>%\n    arrange(GENE, RSID1, RSID2)\n  cov_df\n}\n\n# Refine eventually: will want to make the fields defined here to be optional input parameters (maf, n_folds, etc.), and behave similarly to seed where there is a default input unless user-defined.\nmain <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file,\n                  chrom, prefix, summary_path, weights_path, covariances_path, maf=0.01, n_folds=10, n_train_test_folds=5,\n                 seed=NA, cis_window=1e6, alpha=0.5, null_testing=FALSE) {\n  gene_annot <- get_gene_annotation(gene_annot_file, chrom)\n  expr_df <- get_gene_expression(expression_file, gene_annot)\n  samples <- rownames(expr_df)\n  n_samples <- length(samples)\n  genes <- colnames(expr_df)\n  n_genes <- length(expr_df)\n  snp_annot <- get_filtered_snp_annot(snp_annot_file)\n  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples)\n  #covariates_df <- get_covariates(covariates_file, samples)\n  \n  #Update: order all subject-level data frames to have sampleids in the same order as the expr_df\n  gt_df = gt_df[match(samples, rownames(gt_df)),]\n  #covariates_df = covariates_df[match(samples, rownames(covariates_df)),]\n  \n  # Set seed----\n  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)\n  set.seed(seed)\n  \n  # Prepare output data----\n  model_summary_file <- summary_path %&% prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'\n  model_summary_cols <- c('chrom','gene_id', 'gene_name', 'gene_type', 'alpha', 'n_snps_in_window', 'n_snps_in_model', 'lambda_min_mse',\n                          'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',\n                          'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',\n                          'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')\n  write(model_summary_cols, file = model_summary_file, ncol = 25, sep = '\\t')\n  \n  weights_file <- weights_path %&% prefix %&% '_chr' %&% chrom %&% '_weights.txt'\n  weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')\n  write(weights_col, file = weights_file, ncol = 6, sep = '\\t')\n  \n  tiss_chr_summ_f <- summary_path %&% prefix %&% '_chr' %&% chrom %&% '_summary.txt'\n  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')\n  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)\n  colnames(tiss_chr_summ) <- tiss_chr_summ_col\n  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\\t')\n  \n  covariance_file <- covariances_path %&% prefix %&% '_chr' %&% chrom %&% '_covariances.txt'\n  covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')\n  write(covariance_col, file = covariance_file, ncol = 4, sep = '\\t')\n  \n  # #for testing line by line only\n  # i = 1\n  \n  # Attempt to build model for each gene----\n  for (i in 1:n_genes) {\n    cat(i, \"/\", n_genes, \"\\n\")\n    gene <- genes[i]\n    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]\n    gene_type <- get_gene_type(gene_annot, gene)\n    coords <- get_gene_coords(gene_annot, gene)\n    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)\n    if (all(is.na(cis_gt))) {\n      # No snps within window for gene.\n      model_summary <- c(chrom,gene, gene_name, gene_type, alpha, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)\n      write(model_summary, file = model_summary_file, append = TRUE, ncol = 25, sep = '\\t')\n      next\n    }\n    model_summary <- c(chrom,gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)\n    if (ncol(cis_gt) >= 2) {\n      # expression_vec <- expr_df[,i]\n      # adj_expression <- adjust_for_covariates(expression_vec, covariates_df) #cautious using the adjust_for_covariates function because it assumes covariates and expression have same sample id order/sorting\n      \n      adj_expression1 <- expr_df[,i, drop = F] # use this instead of the adjust for covariates so we don't adjust twice\n      \n      # will try to center and scale residuals as that was done in adjust for covariates function\n      adj_expression <- scale(adj_expression1, center = TRUE, scale = TRUE)\n      \n      adj_expression <- as.matrix(adj_expression[(rownames(adj_expression) %in% rownames(cis_gt)),])\n      \n      if (null_testing)\n        adj_expression <- sample(adj_expression)\n      perf_measures <- nested_cv_elastic_net_perf(cis_gt, adj_expression, n_samples, n_train_test_folds, n_folds, alpha, samples)\n      R2_avg <- perf_measures$R2_avg\n      R2_sd <- perf_measures$R2_sd\n      pval_est <- perf_measures$pval_est\n      rho_avg <- perf_measures$rho_avg\n      rho_se <- perf_measures$rho_se\n      rho_zscore <- perf_measures$rho_zscore\n      rho_avg_squared <- perf_measures$rho_avg_squared\n      zscore_pval <- perf_measures$zscore_pval\n      # Fit on all data\n       # Parallel\n        library(doMC)\n        registerDoMC(cores = 16)\n      cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)\n      fit <- tryCatch(cv.glmnet(cis_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE,  parallel = TRUE),\n                      error = function(cond) {message('Error'); message(geterrmessage()); list()})\n      if (length(fit) > 0) {\n        cv_R2_folds <- rep(0, n_folds)\n        cv_corr_folds <- rep(0, n_folds)\n        cv_zscore_folds <- rep(0, n_folds)\n        cv_pval_folds <- rep(0, n_folds)\n        best_lam_ind <- which.min(fit$cvm)\n        for (j in 1:n_folds) {\n          fold_idxs <- which(cv_fold_ids == j)\n          adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]\n          cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)\n          cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)\n          cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation\n          cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))\n        }\n        cv_R2_avg <- mean(cv_R2_folds)\n        cv_R2_sd <- sd(cv_R2_folds)\n        adj_expr_pred <- predict(fit, as.matrix(cis_gt), s = 'lambda.min')\n        training_R2 <- calc_R2(adj_expression, adj_expr_pred)\n        \n        cv_rho_avg <- mean(cv_corr_folds)\n        cv_rho_se <- sd(cv_corr_folds)\n        cv_rho_avg_squared <- cv_rho_avg**2\n        # Stouffer's method for combining z scores.\n        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)\n        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)\n        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)\n        \n        if (fit$nzero[best_lam_ind] > 0) {\n          \n          weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]\n          weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]\n          weighted_snps_info <- snp_annot %>% filter(varID %in% weighted_snps) %>% select(rsid, varID, ref_vcf, alt_vcf)\n          weighted_snps_info$gene <- gene\n          weighted_snps_info <- weighted_snps_info %>%\n            merge(data.frame(weights = weights, varID=weighted_snps), by = 'varID') %>%\n            select(gene, rsid, varID, ref_vcf, alt_vcf, weights)\n          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\\t')\n          covariance_df <- do_covariance(gene, cis_gt, weighted_snps_info$rsid, weighted_snps_info$varID)\n          write.table(covariance_df, file = covariance_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = \"\\t\")\n          model_summary <- c(chrom, gene, gene_name, gene_type, alpha, ncol(cis_gt), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,\n                             rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)\n        } else {\n          model_summary <- c(chrom, gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,\n                             cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,\n                             cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)\n        }\n      } else {\n        model_summary <- c(chrom, gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, R2_avg, R2_sd, NA, NA, NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,\n                           NA, NA, NA, NA, NA, NA)\n      }\n    }\n    write(model_summary, file = model_summary_file, append = TRUE, ncol = 25, sep = '\\t')\n  }\n}\n\n#############################################\n# 2. Define inputs and calls to run functions for each chromosome\n\n# define which chromosomes to run entered by user\nsource(\"cwl_inputs.R\")\n\ncat(chrom)\n \nsnp_annot_file <- \"snp_annot.chr\" %&% chrom %&% \".txt\"\ngenotype_file <- \"genotype.chr\" %&% chrom %&% \".txt\"\n\nsummary_path <- \"summary/\"\ncovariances_path <- \"covariances/\"\nweights_path <- \"weights/\"\n\n\nmain(snp_annot_file, gene_annot_file, genotype_file, expression_file, seed=seed, maf=as.numeric(maf), n_folds=as.numeric(n_folds), n_train_test_folds=as.numeric(n_train_test_folds), alpha=as.numeric(alpha), chrom=as.numeric(chrom), prefix, summary_path, weights_path, covariances_path, cis_window=1e6, null_testing=FALSE)\n\n\n\n# Check inputs\n'maf' = maf\n'n_folds' = n_folds\n'n_train_test_folds' = n_train_test_folds\n'alpha' = alpha\n'seed' = seed\n",
                                "writable": false
                            },
                            {
                                "entryname": "elnet.sh",
                                "entry": "\nmkdir summary\nmkdir weights\nmkdir covariances\n\nRscript gtex_v7_nested_cv_elnet_training_combined.R",
                                "writable": false
                            },
                            {
                                "entryname": "cwl_inputs.R",
                                "entry": "\nchrom = $(inputs.chr_list)\n\ngene_annot_file = \"$(inputs.gene_annotation.path)\"\nexpression_file = \"$(inputs.adjusted_expression_file.path)\"\n\nexpression_path = \"$(inputs.adjusted_expression_file.path)\"\n\n\npop_name = \"$(inputs.population_name)\"\ntiss_name = \"$(inputs.tissue_name)\"\nprefix = \"$(inputs.model_prefix)\"\nseed = \"$(inputs.seed)\"\nmaf = \"$(inputs.MAF)\"\nn_folds = \"$(inputs.n_folds)\"\nn_train_test_folds = \"$(inputs.n_train_test_folds)\"\nalpha =  \"$(inputs.alpha)\"\n",
                                "writable": false
                            },
                            "$(inputs.snp_annotation)",
                            "$(inputs.genotype_file)"
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.R"
                    },
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.sh"
                    },
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.Rda"
                    },
                    {
                        "class": "sbg:SaveLogs",
                        "value": "standard.out"
                    }
                ],
                "stdout": "standard.out",
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:content_hash": "aad4bc0e7fb2c2a61018ba625032577ab38ad53175d960ebba651ac3f19a103f3",
                "sbg:contributors": [
                    "e.esquinca"
                ],
                "sbg:createdBy": "e.esquinca",
                "sbg:createdOn": 1613684407,
                "sbg:id": "rk.johnson/predictdb-development/model-creation/50",
                "sbg:image_url": null,
                "sbg:latestRevision": 50,
                "sbg:modifiedBy": "e.esquinca",
                "sbg:modifiedOn": 1628791479,
                "sbg:project": "rk.johnson/predictdb-development",
                "sbg:projectName": "Predictdb Development",
                "sbg:publisher": "sbg",
                "sbg:revision": 50,
                "sbg:revisionNotes": "",
                "sbg:revisionsInfo": [
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613684407,
                        "sbg:revision": 0,
                        "sbg:revisionNotes": "Copy of rk.johnson/predictdb/model-creation/47"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613684552,
                        "sbg:revision": 1,
                        "sbg:revisionNotes": "edited script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613684971,
                        "sbg:revision": 2,
                        "sbg:revisionNotes": "adjusted expression file input in cwl"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613685129,
                        "sbg:revision": 3,
                        "sbg:revisionNotes": "add code to part 3"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613685196,
                        "sbg:revision": 4,
                        "sbg:revisionNotes": "change name to know the difference"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613687567,
                        "sbg:revision": 5,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613689092,
                        "sbg:revision": 6,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614021874,
                        "sbg:revision": 7,
                        "sbg:revisionNotes": "took out adjusting twice for covariates"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614040923,
                        "sbg:revision": 8,
                        "sbg:revisionNotes": "fixed seed to 2421 for comparison"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614042657,
                        "sbg:revision": 9,
                        "sbg:revisionNotes": "same edit as last time"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614101174,
                        "sbg:revision": 10,
                        "sbg:revisionNotes": "removed lapply"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614111498,
                        "sbg:revision": 11,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614111556,
                        "sbg:revision": 12,
                        "sbg:revisionNotes": "readded double adjust line 235"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614484814,
                        "sbg:revision": 13,
                        "sbg:revisionNotes": "remove adjust for covariates"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614488470,
                        "sbg:revision": 14,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614489468,
                        "sbg:revision": 15,
                        "sbg:revisionNotes": "cleaned up script"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614489471,
                        "sbg:revision": 16,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614489558,
                        "sbg:revision": 17,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614527381,
                        "sbg:revision": 18,
                        "sbg:revisionNotes": "added back in adjust for covariates"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614640295,
                        "sbg:revision": 19,
                        "sbg:revisionNotes": "removed adjust for covariates"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614705235,
                        "sbg:revision": 20,
                        "sbg:revisionNotes": "scale and center"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614799550,
                        "sbg:revision": 21,
                        "sbg:revisionNotes": "updated script so sample ID's match"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614829387,
                        "sbg:revision": 22,
                        "sbg:revisionNotes": "set same seed"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614833511,
                        "sbg:revision": 23,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614833563,
                        "sbg:revision": 24,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617211818,
                        "sbg:revision": 25,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617649134,
                        "sbg:revision": 26,
                        "sbg:revisionNotes": "added parallel code"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617649340,
                        "sbg:revision": 27,
                        "sbg:revisionNotes": "update parallel code"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617649553,
                        "sbg:revision": 28,
                        "sbg:revisionNotes": "added extra ports"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617650205,
                        "sbg:revision": 29,
                        "sbg:revisionNotes": "added extra ports"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617650690,
                        "sbg:revision": 30,
                        "sbg:revisionNotes": "Added some descriptions"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617650857,
                        "sbg:revision": 31,
                        "sbg:revisionNotes": "Updated app info"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617651437,
                        "sbg:revision": 32,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617651542,
                        "sbg:revision": 33,
                        "sbg:revisionNotes": "delete unnecessary"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617821262,
                        "sbg:revision": 34,
                        "sbg:revisionNotes": "added chrom column for analysis"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617821440,
                        "sbg:revision": 35,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617834672,
                        "sbg:revision": 36,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617899297,
                        "sbg:revision": 37,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618432345,
                        "sbg:revision": 38,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618631011,
                        "sbg:revision": 39,
                        "sbg:revisionNotes": "Remove Covariates file"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618631210,
                        "sbg:revision": 40,
                        "sbg:revisionNotes": "Edited Readme"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618663885,
                        "sbg:revision": 41,
                        "sbg:revisionNotes": "remove covariates file path"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618803112,
                        "sbg:revision": 42,
                        "sbg:revisionNotes": "remove covariates file in main"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618803456,
                        "sbg:revision": 43,
                        "sbg:revisionNotes": "remove all covariate code"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618804194,
                        "sbg:revision": 44,
                        "sbg:revisionNotes": "edit inputs"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618810016,
                        "sbg:revision": 45,
                        "sbg:revisionNotes": "edit inputs"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618861365,
                        "sbg:revision": 46,
                        "sbg:revisionNotes": "remove coriate code"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618873969,
                        "sbg:revision": 47,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1618874013,
                        "sbg:revision": 48,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1628790945,
                        "sbg:revision": 49,
                        "sbg:revisionNotes": "added descriptions"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1628791479,
                        "sbg:revision": 50,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": []
            },
            "label": "Elastic Models",
            "scatter": [
                "chr_list"
            ],
            "sbg:x": 441.701904296875,
            "sbg:y": -163.73800659179688
        },
        {
            "id": "database_summary",
            "in": [
                {
                    "id": "population_name",
                    "source": "population_name"
                },
                {
                    "id": "tissue_name",
                    "source": "tissue_name"
                },
                {
                    "id": "summary",
                    "source": [
                        "model_bulding/summary"
                    ]
                },
                {
                    "id": "covariances",
                    "source": [
                        "model_bulding/covariances"
                    ]
                },
                {
                    "id": "weights",
                    "source": [
                        "model_bulding/weights"
                    ]
                }
            ],
            "out": [
                {
                    "id": "db_output"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "rk.johnson/predictdb/database-summary/33",
                "baseCommand": [
                    "bash summary.sh"
                ],
                "inputs": [
                    {
                        "id": "population_name",
                        "type": "string?",
                        "doc": "Population from which the samples originated from."
                    },
                    {
                        "id": "tissue_name",
                        "type": "string?",
                        "doc": "The tissue the expression was measured in"
                    },
                    {
                        "id": "summary",
                        "type": "File[]?",
                        "doc": "Summary files from the Nested Elastic Net Models split by chromosome that output from the tool."
                    },
                    {
                        "id": "covariances",
                        "type": "File[]?",
                        "doc": "Covariance files from the Nested Elastic Net Models split by chromosome that output from the tool."
                    },
                    {
                        "id": "weights",
                        "type": "File[]?",
                        "doc": "Weight files from the Nested Elastic Net Models split by chromosome that output from the tool."
                    }
                ],
                "outputs": [
                    {
                        "id": "db_output",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "dbs/*.db"
                        }
                    }
                ],
                "doc": "After the Nested Elastic Net Models are run for each chromosome, we take the summary, covariance, and weights file. \nCombine models for all chromosomes entered in one database file. Then the database will be filtered for use in PrediXcan",
                "label": "Database_Summary",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sb.biodatacatalyst.nhlbi.nih.gov/dave/predictdb:v2021_02_13"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "summary.R",
                                "entry": "#############################################\nsuppressMessages(library(dplyr))\nsuppressMessages(library(glmnet))\nsuppressMessages((library(reshape2)))\nsuppressMessages(library(methods))\nsuppressMessages(library(RSQLite))\nsuppressMessages(library(data.table))\n\n\"%&%\" <- function(a,b) paste(a,b, sep='')\n\n# 3. Combine summaries of all chromosomes, and create database of results\n\nsummary_path <- \"summary/\"\ncovariances_path <- \"covariances/\"\nweights_path <- \"weights/\"\n\nsource(\"cwl_inputs.R\")\n\n\n\ndriver <- dbDriver('SQLite')\n\nfilenames = list.files(path = summary_path, pattern = \"*_model_summaries.txt\", full.names = TRUE)\nmodel_summaries <- rbindlist(lapply(filenames, fread))\n\nfilenames2 =list.files(path = summary_path, pattern = \"*_summary.txt\", full.names = TRUE)\ntiss_summary = rbindlist(lapply(filenames2, fread))\n\nn_samples = unique(tiss_summary$n_samples)\n\nmodel_summaries <- rename(model_summaries, gene = gene_id)\nconn <- dbConnect(drv = driver, 'dbs/'%&%pop_name%&%'_'%&%tiss_name%&%'_ncvelnet_'%&%n_samples%&% '.db')\ndbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)\ndbExecute(conn, \"CREATE INDEX gene_model_summary ON model_summaries (gene)\")\n\n# Weights Table -----\nfilenames3 = list.files(path = weights_path, pattern = \"*_weights.txt\", full.names = TRUE)\nweights = rbindlist(lapply(filenames3, fread))\nweights <- rename(weights, gene = gene_id)\ndbWriteTable(conn, 'weights', weights, overwrite = TRUE)\ndbExecute(conn, \"CREATE INDEX weights_rsid ON weights (rsid)\")\ndbExecute(conn, \"CREATE INDEX weights_gene ON weights (gene)\")\ndbExecute(conn, \"CREATE INDEX weights_rsid_gene ON weights (rsid, gene)\")\n\n# Sample_info Table ----\nsample_info <- data.frame(n_samples = n_samples, population = pop_name, tissue = tiss_name)\ndbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)\n\n# Construction Table ----\nconstruction <- tiss_summary %>%\n  select(chrom, cv_seed) %>%\n  rename(chromosome = chrom)\ndbWriteTable(conn, 'construction', construction, overwrite = TRUE)\ndbDisconnect(conn)\n\n\n\n# Filter the databases to get significant values\nunfiltered_db <- 'dbs/' %&% pop_name %&% '_' %&% tiss_name %&% '_ncvelnet_' %&% n_samples %&%  '.db'\nfiltered_db <-  'dbs/' %&% pop_name %&% '_' %&% tiss_name %&% '_ncvelnet_' %&% n_samples %&% '_filtered_signif' %&% '.db'\n\nin_conn <- dbConnect(driver, unfiltered_db)\nout_conn <- dbConnect(driver, filtered_db)\n\nmodel_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg_squared > 0.01')\nmodel_summaries <- model_summaries %>% \n  rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)\nmodel_summaries$pred.perf.qval <- NA\ndbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)\n\nconstruction <- dbGetQuery(in_conn, 'select * from construction')\ndbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)\n\nsample_info <- dbGetQuery(in_conn, 'select * from sample_info')\ndbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)\n\nweights <- dbGetQuery(in_conn, 'select * from weights')\nweights <- weights %>%\n    filter(gene %in% model_summaries$gene) %>%\n    rename(eff_allele = alt, ref_allele = ref, weight = beta)\ndbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)\n\ndbExecute(out_conn, \"CREATE INDEX weights_rsid ON weights (rsid)\")\ndbExecute(out_conn, \"CREATE INDEX weights_gene ON weights (gene)\")\ndbExecute(out_conn, \"CREATE INDEX weights_rsid_gene ON weights (rsid, gene)\")\ndbExecute(out_conn, \"CREATE INDEX gene_model_summary ON extra (gene)\")\n\ndbDisconnect(in_conn)\ndbDisconnect(out_conn)\n",
                                "writable": false
                            },
                            {
                                "entryname": "cwl_inputs.R",
                                "entry": "pop_name = \"$(inputs.population_name)\"\ntiss_name = \"$(inputs.tissue_name)\"",
                                "writable": false
                            },
                            {
                                "entryname": "summary.sh",
                                "entry": "mkdir dbs\nmkdir weights\nmkdir summary\nmkdir covariances\n\nmv *weights.txt weights\nmv *model_summaries.txt summary\nmv *_summary.txt summary\nmv *covariances.txt covariances\n\nRscript summary.R",
                                "writable": false
                            },
                            "$(inputs.summary)",
                            "$(inputs.covariances)",
                            "$(inputs.weights)"
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.R"
                    },
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.Rda"
                    },
                    {
                        "class": "sbg:SaveLogs",
                        "value": "*.sh"
                    },
                    {
                        "class": "sbg:SaveLogs",
                        "value": "standard.out"
                    }
                ],
                "stdout": "standard.out",
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:content_hash": "af62c9ec347cc1e56438178451672aa86dcaa7d281216ac9132e00b6d13876ddb",
                "sbg:contributors": [
                    "dave",
                    "e.esquinca"
                ],
                "sbg:createdBy": "e.esquinca",
                "sbg:createdOn": 1612889428,
                "sbg:id": "rk.johnson/predictdb/database-summary/33",
                "sbg:image_url": null,
                "sbg:latestRevision": 33,
                "sbg:modifiedBy": "e.esquinca",
                "sbg:modifiedOn": 1628791365,
                "sbg:project": "rk.johnson/predictdb",
                "sbg:projectName": "predictdb",
                "sbg:publisher": "sbg",
                "sbg:revision": 33,
                "sbg:revisionNotes": "",
                "sbg:revisionsInfo": [
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612889428,
                        "sbg:revision": 0,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1612889505,
                        "sbg:revision": 1,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613688240,
                        "sbg:revision": 2,
                        "sbg:revisionNotes": "added inputs"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613688788,
                        "sbg:revision": 3,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613688926,
                        "sbg:revision": 4,
                        "sbg:revisionNotes": "added directories from previous tool"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1613752519,
                        "sbg:revision": 5,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1613752825,
                        "sbg:revision": 6,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1613753925,
                        "sbg:revision": 7,
                        "sbg:revisionNotes": "mv"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613760110,
                        "sbg:revision": 8,
                        "sbg:revisionNotes": "added library"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1613762514,
                        "sbg:revision": 9,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1613769966,
                        "sbg:revision": 10,
                        "sbg:revisionNotes": "removed unnecessary args files from command line."
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1613769968,
                        "sbg:revision": 11,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1613770124,
                        "sbg:revision": 12,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614109689,
                        "sbg:revision": 13,
                        "sbg:revisionNotes": "added dbs/"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614109953,
                        "sbg:revision": 14,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614187886,
                        "sbg:revision": 15,
                        "sbg:revisionNotes": "commented out code"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614188620,
                        "sbg:revision": 16,
                        "sbg:revisionNotes": "save logs"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614188724,
                        "sbg:revision": 17,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614189634,
                        "sbg:revision": 18,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614190553,
                        "sbg:revision": 19,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614192972,
                        "sbg:revision": 20,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614194099,
                        "sbg:revision": 21,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614194811,
                        "sbg:revision": 22,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614195226,
                        "sbg:revision": 23,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614195711,
                        "sbg:revision": 24,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1614262773,
                        "sbg:revision": 25,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1614277370,
                        "sbg:revision": 26,
                        "sbg:revisionNotes": "fixing dbs output"
                    },
                    {
                        "sbg:modifiedBy": "dave",
                        "sbg:modifiedOn": 1614277382,
                        "sbg:revision": 27,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1617821777,
                        "sbg:revision": 28,
                        "sbg:revisionNotes": "changed t0 (rho_avg_squared) > 0.01"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1619632116,
                        "sbg:revision": 29,
                        "sbg:revisionNotes": "rename databases"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1619647551,
                        "sbg:revision": 30,
                        "sbg:revisionNotes": "update name"
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1621981733,
                        "sbg:revision": 31,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1622131779,
                        "sbg:revision": 32,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:modifiedBy": "e.esquinca",
                        "sbg:modifiedOn": 1628791365,
                        "sbg:revision": 33,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": []
            },
            "label": "Database_Summary",
            "sbg:x": 923.9014892578125,
            "sbg:y": -188.00721740722656
        }
    ],
    "requirements": [
        {
            "class": "ScatterFeatureRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "sbg:image_url": "https://platform.sb.biodatacatalyst.nhlbi.nih.gov/ns/brood/images/dave/predixscan/predictdb-main/1.png",
    "sbg:projectName": "prediXscan",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedBy": "dave",
            "sbg:modifiedOn": 1645031328,
            "sbg:revisionNotes": null
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedBy": "dave",
            "sbg:modifiedOn": 1645031372,
            "sbg:revisionNotes": "changed specific parameters"
        }
    ],
    "sbg:appVersion": [
        "v1.2",
        "v1.1"
    ],
    "sbg:id": "dave/predixscan/predictdb-main/1",
    "sbg:revision": 1,
    "sbg:revisionNotes": "changed specific parameters",
    "sbg:modifiedOn": 1645031372,
    "sbg:modifiedBy": "dave",
    "sbg:createdOn": 1645031328,
    "sbg:createdBy": "dave",
    "sbg:project": "dave/predixscan",
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [],
    "sbg:contributors": [
        "dave"
    ],
    "sbg:latestRevision": 1,
    "sbg:publisher": "sbg",
    "sbg:content_hash": "ad2e0faffd7831337f3964dfba52df178e8688896f434ae968aa2ebfc45a5d7bd",
    "sbg:workflowLanguage": "CWL"
}