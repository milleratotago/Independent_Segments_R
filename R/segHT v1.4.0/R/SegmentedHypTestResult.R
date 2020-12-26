# Class SegmentedHypTestResult

# This class holds information about the expected outcomes for a segmented hypothesis testing experiment for given
# SegmentedHypTestResult (kmax, alphaTotal, alphaWeak and alphaStrong) & StatType (e.g. one sample t-test or correlation r).

# Can be produced either for a single true effect size or for a probability-weighted average of
# multiple true effect sizes which may include and effect size of 0 (i.e. H0 is true) with associated baserate probability.

# Note that some values are vectors holding a result for each of the possible ending segments 1 to k.

#############################################################
# Class properties - all default to 0 at creation
#############################################################
# max_n_segments:         maximum number of segments that can be run in the described experiment before termination
# pr_true_pos:            hit - reject Ho when it is false (there is an effect)
# pr_false_pos:           false alarm - reject Ho when it is true (there is no effect)
# pr_true_neg:            correct failure to reject
# pr_false_neg:           miss - failure to reject when there is really an effect
# pr_reject_total:        probability of rejecting (correctly or incorrectly) at any of the k segments
# exp_n_subj:             expected (mean across inf iterations) number of subjects used when testing for this effect.
# exp_n_subj_sqr:         expected (mean across inf iterations) square of the number of subjects used. Useful for variance computations.
# exp_effect_size:        expected (mean across inf iterations) observed effect size (not generally equal to the true effect size)
# pr_reject_by_segment    vector giving prob of rejecting at each segment
# pr_ftr_by_segment       vector giving probability of stopping with FTR at each segment
# pr_continue_by_segment  vector giving probability of continuing at each segment
#
# n_simulations           0 for computed results, or number of simulations for simulated results (not used in this project)

# The public-facing constructor requires only max_n_segments, and sets all other fields to 0
# but if you want to pack an instance with a set of prepared result values, you can assign to l0 the
# scalar input args.

#####################################################################
# Creation methods - following Wickham
#####################################################################
new_SegmentedHypTestResult <- function(max_n_segments,
                                       pr_true_pos = 0,
                                       pr_false_pos = 0,
                                       pr_true_neg = 0,
                                       pr_false_neg = 0,
                                       pr_reject_total = 0,
                                       exp_n_subj = 0,
                                       exp_n_subj_sqr = 0,
                                       exp_effect_size = 0,
                                       n_simulations = 0)
{

  new_segmented_hyp_test_result <- list(max_n_segments = max_n_segments,
                                        pr_true_pos = pr_true_pos,
                                        pr_false_pos = pr_false_pos,
                                        pr_true_neg = pr_true_neg,
                                        pr_false_neg = pr_false_neg,
                                        pr_reject_total = pr_reject_total,
                                        exp_n_subj = exp_n_subj,
                                        exp_n_subj_sqr = exp_n_subj_sqr,
                                        exp_effect_size = exp_effect_size,
                                        n_simulations = n_simulations,
                                        pr_reject_by_segment = numeric(max_n_segments),
                                        pr_ftr_by_segment = numeric(max_n_segments),
                                        pr_continue_by_segment = numeric(max_n_segments))



  attr(new_segmented_hyp_test_result, "class") <- "SegmentedHypTestResult"

  # Return instance
  return(new_segmented_hyp_test_result)
}


#----------------------------------------------------------------------
validate_SegmentedHypTestResult <- function(segemented_hyp_test_result)
{
  # Perform additional checking here, if wanted.
}


#---------------------------------------------------------------------------
#-----------------
# Roxygen comments

#' SegmentedHypTestResult Constructor
#'
#' A SegmentedHypTestResult instance holds information about the expected
#' outcomes of a segmented hypothesis testing study for given maximum number of
#' segments, alpha total, alpha strong, and statistical procedure.These objects
#' are typically created with no parameter values supplied so all fields default
#' to 'empty', and initisalised during computation to hold the simulated or
#' numerically-derived expected outcomes for the scenario. See
#' segHT::run_scenario.SegmentedHypTestEngine for an example.
#'
#' @param max_n_segments Maximum number of segments that can be run in the
#'   described experiment before termination
#' @param pr_true_pos Overall probability of a hit - reject Ho when it is false
#'   (there is an effect)
#' @param pr_false_pos Overall probability of a false alarm - reject Ho when it is
#'   true (there is no effect)
#' @param pr_true_neg Overall probability of a correct failure to reject
#' @param pr_false_neg Overall probability of a miss - failure to reject when there is really
#'   an effect
#' @param pr_reject_total Overall probability of rejecting (correctly or
#'   incorrectly) at any of the k segments
#' @param exp_n_subj Expected (mean across infinite iterations) number of
#'   subjects used when testing for this effect.
#' @param exp_n_subj_sqr Expected (mean across infinite iterations) square of
#'   the number of subjects used. Useful for variance computations.
#' @param exp_effect_size Expected (mean across infinite iterations) observed
#'   effect size (not generally equal to the true effect size)
#' @param n_simulations 0 for computed results, or the number of
#'   simulations for simulated results (to allow extension to simulation)
#' @return SegmentedHypTestResult instance
#'
#' @export
SegmentedHypTestResult <- function(max_n_segments,
                                    pr_true_pos = 0,
                                    pr_false_pos = 0,
                                    pr_true_neg = 0,
                                    pr_false_neg = 0,
                                    pr_reject_total = 0,
                                    exp_n_subj = 0,
                                    exp_n_subj_sqr = 0,
                                    exp_effect_size = 0,
                                    n_simulations = 0)
{

  # Create instance
  instance = new_SegmentedHypTestResult(max_n_segments,
                                        pr_true_pos,
                                        pr_false_pos,
                                        pr_true_neg,
                                        pr_false_neg,
                                        pr_reject_total,
                                        exp_n_subj,
                                        exp_n_subj_sqr,
                                        exp_effect_size,
                                        n_simulations)

  # Validate instance. Execution halts here if things go wrong
  validate_SegmentedHypTestResult(instance)

  # Return validated instance
  return(instance)
}


#####################################################################
# Public methods
#####################################################################
#========================================================================================================
# print.SegmentendHypTestResult: System will dispatch to this method on print(instance)
#========================================================================================================
#' @export
print.SegmentedHypTestResult <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called

  # Print only the scalars, since this is primarily for debugging anyway. If you include the
  # vectors, R replicates the scalars to fill the data frame and it is confusing.
  # The scalars are the first 10 fields.

  n_scalar_fields <- 10
  df <- data.frame(data_values[1:n_scalar_fields])
  print(df)
}








