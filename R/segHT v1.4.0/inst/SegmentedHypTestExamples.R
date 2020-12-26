# See seghtQuickStart.Rmd for detailed explanations


#' @importFrom utils install.packages


segHT_manual_setup <- function()
{
  rm(list = ls())

  # Install packages if necessary
  if (!require(pracma)) install.packages('pracma')
  if (!require(stats)) install.packages('stats')
  if (!require(SuppDists)) install.packages('SuppDists')

  library(pracma)
  library(SuppDists)
  library(stats)


  source("R/SegmentedHypTestResult.R")
  source("R/SegmentedResearcher.R")
  source("R/StatProcedureBase.R")
  source("R/StatProcedureOneSampleT.R")
  source("R/StatProcedureTwoSampleT.R")
  source("R/StatProcedurePearsonR.R")
  source("R/StatProcedureOneSampleZ.R")
  source("R/StatProcedureTwoSampleZ.R")
  source("R/TrueEffect.R")
  source("R/SegmentedHypTestEngine.R")
  source("R/SegmentedHypTestUtilities.R")
}



# segHT Utility Functions
#
# alpha_weak: Accepts maximum number of segments, alpha total and alpha strong. Returns the required alpha weak.
#
# n_for_power: Accepts a full set of study parameters (see example below) and a desired power level. Returns the number of subjects per segment required to achieve the target power.
#
# segmented_hyp_test_outcomes: Accepts a full set of study parameters (see example below). Returns an overall summary of the expected performance of the segmented hypothesis testing protocol including probabilities of true and false positives and negatives, probabilities of termination at each segment, expected number of subjects tested, etc.
#
# search_kmax: Accepts a full set of study parameters **except for maximum number of segments** and returns outcome summaries for each possible number of segments between 2 and 10. These results help the researcher to select the optimal number of segments in conjunction with practical constraints.
#
# See detailed code examples of these functions below.

################ alpha_weak ################
alpha_weak_demo <- function()
{
  # Initialise value parameters
  max_n_segments <- 3
  alpha_total <- 0.05
  alpha_strong <- 0.025

  # call function
  required_alpha_weak <- alpha_weak(max_n_segments, alpha_total, alpha_strong)

  # display output
  print(required_alpha_weak)
}


################ n_for_power ################
n_for_power_demo <- function()
{
  # Identify hypothesis test to be used. Values are 1t, 2t, r, 1z, 2z
  stat_procedure_name <- "1t"

  # Initialise value parameters. NB: the value of alpha_weak will be computed for you automatically
  max_n_segments <- 3
  alpha_total <- 0.05
  alpha_strong <- 0.025

  # State the expected effect size
  effect_size <- 0.5

  # State the desired power
  target_power <- 0.8

  # Call the function
  required_n_per_segment <- n_for_power(target_power,
                                        stat_procedure_name,
                                        effect_size,
                                        max_n_segments,
                                        alpha_total,
                                        alpha_strong)


  # Display results
  print(required_n_per_segment)
}


################ segmented_hyp_test_outcomes ################
segmented_hyp_test_outcomes_demo <- function()
{
  #########################################################
  # With default base rate effect probability = 1

  # Identify hypothesis test to be used. Values are 1t, 2t, r, 1z, 2z
  stat_procedure_name <- "1t"

  # Initialise value parameters. NB: the value of alpha_weak will be computed for you automatically.
  max_n_segments <- 3
  n_per_segment <- 20
  alpha_total <- 0.05
  alpha_strong <- 0.025

  # State the expected effect size.
  effect_sizes <- 0.2

  # Call the function.
  # When not provided, the final base_rate parameter defaults to 1
  seght_outcomes <- segmented_hyp_test_outcomes(max_n_segments,
                                                n_per_segment,
                                                alpha_total,
                                                alpha_strong,
                                                stat_procedure_name,
                                                effect_sizes)


  # Look at your results
  print("Outcomes Example 1:")
  print(seght_outcomes)


  #########################################################
  # With specified base rate probability

  # Identify hypothesis test to be used. Values are 1t, 2t, r, 1z, 2z
  stat_procedure_name <- "1t"

  # Initialise value parameters. NB: the value of alpha_weak will be computed for you
  max_n_segments <- 3
  n_per_segment <- 20
  alpha_total <- 0.05
  alpha_strong <- 0.025

  # State the expected effect size and base rate probability.
  effect_sizes <- 0.2
  base_rates <- 0.6

  # Call the method, supplying a value for the base_rate argument
  seght_outcomes <- segmented_hyp_test_outcomes(max_n_segments,
                                                n_per_segment,
                                                alpha_total,
                                                alpha_strong,
                                                stat_procedure_name,
                                                effect_sizes,
                                                base_rates)


  # Look at your results
  print("Outcomes Example 2:")
  print(seght_outcomes)
}


################ search_kmax ################
search_kmax_demo <- function()
{
  # Identify hypothesis test to be used. Values are 1t, 2t, r, 1z, 2z
  stat_procedure_name <- "2t"

  # Initialise value parameters. NB: the value of alpha_weak will be computed for you automatically.
  alpha_total <- 0.05
  alpha_strong <- 0.025

  # State the expected effect size and base rate probability. When the base rate is less than 1, the remaining
  # 1-base rate is assumed to be the probability of cases with effect size = 0 (that is, H0 is true).
  effect_size <- 0.2
  current_base_rate <- 0.6

  # State the desired power level
  target_power <- 0.80

  # Call the method
  kmax_summary <- search_kmax(alpha_total,
                              alpha_strong,
                              stat_procedure_name,
                              target_power,
                              effect_size,
                              base_rate = current_base_rate)

  # Display your results
  print(kmax_summary)
}

