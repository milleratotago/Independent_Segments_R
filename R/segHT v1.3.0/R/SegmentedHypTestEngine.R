# Class SegmentedHypTestEngine

# Runs (computes the proababilistic results of) complete segmented hypothesis testing scenarios

#########################################################################################################
# Creation methods - following Wickham
#########################################################################################################
#========================================================================================================
# new_SegmentedHypTestEngine
#========================================================================================================
new_SegmentedHypTestEngine <- function()
{

  # No fields needed at creation
  new_SegmentedHypTestEngine <- list()


  attr(new_SegmentedHypTestEngine, "class") <- "SegmentedHypTestEngine"

  # Return instance
  return(new_SegmentedHypTestEngine)

} # new_SegmentedHypTestEngine
#========================================================================================================

#========================================================================================================
# validate_SegmentedHypTestEngine
#========================================================================================================
# Checks additional properties of elements, for example, that values are in allowed ranges -- business rules
validate_SegmentedHypTestEngine <- function(segmented_hyp_test_engine)
{
  # More elaborate checking in here if needed. Stop if things go wrong.

} # end validate_SegmentedHypTestEngine
#========================================================================================================

#========================================================================================================
# SegmentedHypTestEngine - Public ctor
#========================================================================================================
#-----------------
# Roxygen comments

#' SegmentedHypTestEngine
#'
#' \code{SegmentedHypTestEngine}  is the public-facing SegmentedHypTestEngine
#' constructor. It accepts no parameters. This object exposes the primary methods
#' for numerically determing the expected outcomes of a segmented hypothesis testing
#' study. See run_scenario.SegmentedHypTestEngine for an example.
#' @return Validated SegmentedHypTestEngine instance
#'
#' @export
SegmentedHypTestEngine <- function()
{
  # Create instance
  instance = new_SegmentedHypTestEngine()

  # Validate instance. Execution halts here if things go wrong
  validate_SegmentedHypTestEngine(instance)

  # Return validated instance
  return(instance)

} # end SegmentedHypTestEngine ctor
#========================================================================================================


#########################################################################################################
# Public methods
#########################################################################################################
#' @export
#========================================================================================================
# print.SegmentedHypTestEngine: System will dispatch to this method on print(instance)
#========================================================================================================
print.SegmentedHypTestEngine <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}


#========================================================================================================
# run_scenario: Top-level wrapper for running a complete scenario.
#
# Accepts researcher, statistical procedure and list (must be a list) of TrueEffect instances. Each TrueEffect
# hold size and baserate probability. A TrueEffect with size 0 (i.e. H0 is true) may be included.
# Probabilities across all effect sizes must sum to 1. A single SegmentedResult object is returned, which is the
# average across all effect sizes, weighted by associated probability.
#========================================================================================================
#-----------------
# Roxygen comments

#' run_scenario.SegmentedHypTestEngine
#'
#' \code{run_scenario.SegmentedHypTestEngine} numerically determines expected
#' outcomes for a complete study scenario. It wraps various internal methods of
#' class SegmentedHypTestEngine (developers can cf. source code).
#' @param segmented_hyp_test_engine A SegmentedHypTestEngine instance for method
#'   dispatch
#' @param segmented_researcher Correctly initialised instance of
#'   SegmentedResearcher
#' @param stat_procedure Correctly initialised instance of any class descended
#'   from StatProcedureBase, e.g. OneSampleT, PearsonR, etc.
#' @param effects_list A list of TrueEffect instances. This must be a list, not
#'   simply a vector. Cast with list() if necessary. Each TrueEffect instance
#'   holds fields effect_size and effect_size_probability. A TrueEffect with
#'   effect_size 0 (i.e. H0 is true) may be included. Probabilities across all effect
#'   sizes must sum to 1 or an exception is thrown.
#' @return A SegmentedHypTestResult instance which holds the average expected outcomes
#' across all effect sizes, weighted by associated probability.
#'
#' @export
run_scenario.SegmentedHypTestEngine <- function(segmented_hyp_test_engine,
                                                segmented_researcher,
                                                stat_procedure,
                                                effects_list)
{
  # Check that effect probs sum to 1. If they do not, throw an exception

  # Grab all the probability values from the true_effects list
  probs <- sapply(effects_list, function(te) te[["effect_size_probability"]])

  if (sum(probs) != 1)
  {
    stop("The sum of all effect size probabilities must be 1")
  }

  # Determine how many effects were passed in
  n_effects <- length(effects_list)

  # Make skeleton to hold results
  individual_results <- vector(mode="list", length = n_effects)

  # Run the experimental scenario (i.e. compute outcome probabilities) for each effect size and accrue results
  for(i in 1:n_effects)
  {
    individual_results[[i]] <- run_any_subjects_scenario(segmented_researcher, stat_procedure, effects_list[[i]])
  }

  # Pass results and probs to weighted averager. Returns a SegmentedHypTestResult instance
  result <- weighted_avg_hyp_test_result(individual_results, probs)

  # Return weighted average SegmentedHypTestResult instance
  return(result)
}

#####################################################################################
# average_power(): By computational definition, the average power across a set of
# experiments. Needs to run from the raw data, as SegHTResults do not carry their
# original true effect size(s), so needs access to the scenario methods.
# Called by user-facing methods in SegmentedHypTester.R
#####################################################################################
#' @export
average_power.SegmentedHypTestEngine <- function(segmented_hyp_test_engine,
                                                 max_n_segments,
                                                 n_per_segment,
                                                 alpha_total,
                                                 alpha_strong,
                                                 stat_procedure,
                                                 effects_list)
{

  # We exclude effect sizes of zero (i.e. when H0 is true) from this power calculation
  # Drop any element of effects_list with effect_size == 0
  real_effects_list <- effects_list[sapply(effects_list, function(te) te[["effect_size"]] > 0)]


  # Skeleton to hold all results, one for each provided effect size greater than 0
  n_real_effects <- length(real_effects_list)
  results_list <- list(n_real_effects)


  # For each provided true effect, run the experimental scenario, producing a SegmentedHypTestResult

  # Same researcher information used throughout
  segmented_researcher <- SegmentedResearcher(max_n_segments,n_per_segment, alpha_total, alpha_strong)

  # We call run_any_subject_scenario here to run the partial analysis with total effect prob < 1
  # Here we just want this single result, with or without baserate, not the weighted average
  # Run for each effect and accrue the results
  for (i in 1:length(real_effects_list))
  {
      curr_effect <- real_effects_list[[i]]
      result <- run_any_subjects_scenario(segmented_researcher, stat_procedure, curr_effect)
      results_list[[i]] <- result
  }

  # Get the average power. Computations by definition (cf. Miller, 2020)

  # We multiply each baserate by its associated probability of a true positive, and sum these.
  # Average power is that sum divided by the total baserate where effect size is greater than 0.

  pr_true_pos <- sapply(results_list, function(result) result[["pr_true_pos"]])[1:n_real_effects]
  effect_probs <- sapply(real_effects_list, function(te) te[["effect_size_probability"]])

  # Perform the computation as described
  power_numerator <- sum(pr_true_pos * effect_probs)
  power_denominator <- sum(effect_probs)
  avg_power <- power_numerator/power_denominator

  return(avg_power)

} # end average_power

#########################################################################################################
# Internal Methods
#########################################################################################################
#========================================================================================================
# outcome_probabilities
#========================================================================================================

# This method computes the probability of each outcome (reject, fail to reject, or continue) at each
# of kmax stages for a given experimental scenario. See Miller & Ulrich, 2020, esp. Figures 2A and 2B for further
# computational details.

# Note that when Ho is true, "pr_beat_strong" (causing E to reject and stop) and "pr_beat_weak" (1 - pr_beat_weak causing
# E to fail to reject and stop) are exactly equal to alpha strong and alpha weak. However, when Ho is false these
# values depend on the effect size and statistical test being used.

outcome_probabilities <- function(max_n_segments,pr_beat_strong,pr_beat_weak)
{
  # The method returns a list with the following fields:

  #   pr_reject_total: Total probability of rejecting at any step during the experiment
  #   pr_reject_by_segment: Probability of stopping after step i & rejecting.
  #   pr_ftr_by_segmenmt: Probability of stopping after step i & failing to reject.
  #   pr_continue_by_segment: Probability of continuing (neither reject or ftr) at step i
  #   expected_n_segments: The mean of the number of steps taken.
  #   variance_n_segments: The variance of the number of steps taken.


  # Prepare skeleton lists to hold outputs
  pr_reject_by_segment <- numeric(max_n_segments)
  pr_ftr_by_segment <- numeric(max_n_segments)
  pr_continue_by_segment <- numeric(max_n_segments)

  # continue if p is neither greater than weak or less than strong
  pr_continue <- pr_beat_weak - pr_beat_strong

  # Compute the probabilistic outcomes for each of max_n_segments steps
  for (curr_segment in 1:max_n_segments)
  {
    # To reject at this segment, you have continued curr_segment-1 previous times,
    # and now you beat (i.e. are less than or equal to) alpha strong.
    # Take the product of the probabilities of these events
    pr_reject_by_segment[curr_segment] <- pr_beat_strong * pr_continue^(curr_segment-1)

    # To fail to reject at this segment, you have continued curr_segment-1 previous times,
    # and now you DO NOT BEAT (i.e. are not less than or equal to; are greater than) pr_beat_weak
    pr_ftr_by_segment[curr_segment] <- (1-pr_beat_weak) * pr_continue^(curr_segment-1)

    # To continue at this segment, you have continued at all previous segments, and this one as well
    pr_continue_by_segment[curr_segment] <- pr_continue^curr_segment
  }

  # The probability of continuing at the last segments is, by definition 0. Adjust that here.
  pr_continue_by_segment[max_n_segments] <- 0;

  excess_prob <- 1 - sum(pr_reject_by_segment) - sum(pr_ftr_by_segment)
  pr_reject_by_segment[max_n_segments] <- pr_reject_by_segment[max_n_segments] + excess_prob

  # Total probabilitiy of rejecting is the sum pr(reject) across all segments
  pr_reject_total <- sum(pr_reject_by_segment)


  # By definition
  expected_n_segments <- (1 - pr_continue^max_n_segments)/(1-pr_continue)

  # This vector holds the probability of terminating at each segment, by either decision
  pr_termination_by_segment <- pr_reject_by_segment + pr_ftr_by_segment

  # By definition
  segment_ord_squared <- (1:max_n_segments)^2
  expected_n_segments_squared <- sum(pr_termination_by_segment * segment_ord_squared)
  variance_n_segments <- expected_n_segments_squared - expected_n_segments^2

  result_list <- list(pr_reject_total = pr_reject_total,
                      pr_reject_by_segment = pr_reject_by_segment,
                      pr_ftr_by_segment = pr_ftr_by_segment,
                      pr_continue_by_segment = pr_continue_by_segment,
                      expected_n_segments = expected_n_segments,
                      variance_n_segments = variance_n_segments)

  return(result_list)

} # end outcome_probabilities

#========================================================================================================
# run_int_subjects_scenario: Run scenario with an integer number of subjects
#========================================================================================================

# This method generates a description of the probabalistic results of a single segmented
# hypothesis testing scenario for a given "researcher" (alphaweak,alphastrong and kmax),
# statistical test and one or more effect sizes whose total probability sums to 1. It returns
# an instance of SegmentedHypTestResult, which is effectivley a table of the useful descriptive
# metrics (see source file).

# This method requires that the provided SegmentedResearcher's number of subjects is an integer. See wrapper
# run_any_subjects_scenario for handling of non-integer cases that can be encountered, for example, during simulations.

run_int_subjects_scenario  <- function(segmented_researcher, stat_procedure, true_effect)
{

  # Create a skeleton instance to hold the results. Number of segments is provided to accomodate
  # the output vectors which hold a value for each segment (e.g. pr(reject) at segment i)
  results <- SegmentedHypTestResult(segmented_researcher$max_n_segments) # all other properties default to 0


  # A SegmentedResearcher object specified the number of subjects to be used for each segment
  # totalled across all groups. Thus the actual size per group depends on the number of groups under
  # the statistical procedure being used. For example, a one-sample t places all subjects in a segment
  # into a single group while a two-sample t puts half in each group. The statistical procedure classes
  # expose method divide_total_sample to perform this calculation

  n_per_group <- divide_total_sample(stat_procedure, segmented_researcher$n_per_segment)

  # Query the statistical procedure instance to compute the probabilities of "beating" (i.e. getting a
  # smaller p observed value than) alpha_strong for the effect size(s). This is logically equivalent to power
  # when Ho is false or the Type I error rate when Ho is true, those outcomes associated with rejecting Ho

  # These values are required to compute the prob of rejecting, failing to reject, or continuing at each
  # segment. When H0 is true (effect size = 0), they are equal to alpha_strong and alpha_weak, respectively.
  # When Ho is false they depend on the statistical procedure, the effect size and number of subjects per group.

  pr_beat_strong <- compute_power(stat_procedure,
                                  segmented_researcher$alpha_strong,
                                  true_effect$effect_size,
                                  n_per_group)

  # Perform the equivalent computation for alpha_weak. Failure to beat alpha_weak triggers failure to reject.

  pr_beat_weak <- compute_power(stat_procedure,
                                  segmented_researcher$alpha_weak,
                                  true_effect$effect_size,
                                  n_per_group)

  # Use these values to compute the various outcome probabilities (for reject, ftr and continue) for each
  # segment.

  outcomes <- outcome_probabilities(segmented_researcher$max_n_segments, pr_beat_strong, pr_beat_weak)

  # The outcome_probabilities method returns a list with the fields as shown. We store some of
  # these in the skeleton results object.

    #   pr_reject_total: Total probability of rejecting at any step during the experiment
    #   pr_reject_by_segment: Probability of stopping after step i & rejecting.
    #   pr_ftr_by_segmenmt: Probability of stopping after step i & failing to reject.
    #   pr_continue_by_segment: Probability of continuing (neither reject or ftr) at step i
    #   expected_n_segments: The mean of the number of steps taken.
    #   variance_n_segments: The variance of the number of steps taken.

  results$pr_reject_total <- outcomes$pr_reject_total
  results$pr_reject_by_segment <- outcomes$pr_reject_by_segment
  results$pr_ftr_by_segment <- outcomes$pr_ftr_by_segment
  results$pr_continue_by_segment <- outcomes$pr_continue_by_segment

  # Set up values of alpha-beta matrix and store in result
  if (true_effect$effect_size == 0) {
    alpha_beta <- c(0, results$pr_reject_total, 1 - results$pr_reject_total, 0)
  } else {
    alpha_beta <- c(results$pr_reject_total, 0, 0, 1 - results$pr_reject_total)
  }

  results$pr_true_pos <- alpha_beta[1] # R vectors are 1-indexed
  results$pr_false_pos <- alpha_beta[2]
  results$pr_true_neg <- alpha_beta[3]
  results$pr_false_neg <- alpha_beta[4]

  # With all these important values computed, the expected numbers of subjects and expected effect size
  # for cases where Ho is rejected can be computed, by the magic of maths, as follows (pers comm, Miller, 2020):

  exp_n_subj <- outcomes$expected_n_segments * segmented_researcher$n_per_segment
  exp_n_segments_sqr <- outcomes$variance_n_segments + outcomes$expected_n_segments^2
  exp_n_subj_sqr <- exp_n_segments_sqr * segmented_researcher$n_per_segment^2

  # Assign the relevant values to the results object
  results$exp_n_subj <- exp_n_subj
  results$exp_n_subj_sqr <- exp_n_subj_sqr

  # Expected size of significant effect, conditional on rejection. Useful for simulations onle.

  # Find the expected effect size (mean of effect size across infinitely many iterations) conditional on
  # a significant result, for a given researcher object (kmax, nsubjects, alphastrong and alpha weak), and
  # statistical procedure. Simulations using this computation have shown that the expected effect size is
  # generally larger than the true effect size. That is, looking only at those experiments where Ho is
  # rejectd (as is typical in the publication process) leads to a general overestimation of true effect size.

  # Adjust sample size for number of groups per experiment (e.g. one-sample has 1, two-sample has 2)
  n_per_group <- divide_total_sample(stat_procedure, segmented_researcher$n_per_segment)

  # Compute probability-weighted average of expected significant effect sizes across all possible strong cut-off stopping segments.
  # Segments preceding the strong cut-off stop have p's bounded between alpha weak and alpha strong.


  # Compute the three expected effect sizes -- against alpha weak, between alpha weak and strong, and against alpha strong
  exp_effect_size_weak <- expected_significant_effect(stat_procedure, segmented_researcher$alpha_weak, n_per_group, true_effect$effect_size)
  exp_effect_size_between <- expected_significant_effect_bounded(stat_procedure, segmented_researcher$alpha_strong, segmented_researcher$alpha_weak, n_per_group, true_effect$effect_size)
  exp_effect_size_strong <- expected_significant_effect(stat_procedure, segmented_researcher$alpha_strong, n_per_group, true_effect$effect_size)

  # The weak result in the baseline nonsegmented case
  if ((segmented_researcher$alpha_strong == 0) || (segmented_researcher$max_n_segments == 1)){
    exp_effect_size_overall <- exp_effect_size_weak
  } else {
    # Iteratively determine the expected effect size for rejection at each of the kmax segments

    # Prepare skeleton to hold results
    exp_effect_size_by_segment <- numeric(segmented_researcher$max_n_segments)

    for(curr_segment in 1:segmented_researcher$max_n_segments - 1)
    {
      exp_effect_size_by_segment[curr_segment] <- (exp_effect_size_strong + ((curr_segment-1) * exp_effect_size_between))/curr_segment
    }
    exp_effect_size_by_segment[segmented_researcher$max_n_segments] <- (exp_effect_size_weak + ((segmented_researcher$max_n_segments -1)* exp_effect_size_between))/segmented_researcher$max_n_segments

    # Compute the overall expected significant effect size as the weighted average of each segment's expected effect size and probability
    # These probabilities have been computed above.
    exp_effect_size_overall <- sum(exp_effect_size_by_segment * outcomes$pr_reject_by_segment)/ outcomes$pr_reject_total
  } # end else


  # Store the value in the results data structre
  results$exp_effect_size <- exp_effect_size_overall

  return(results)

} # end run_int_subjects_scenario

#========================================================================================================
# run_any_subjects_scenarioL Run scenario with integer or real number of subjects. Wrapper for
# run_int_subjects_scenario
#========================================================================================================
# This function handles cases in which Researcher.NperSegment is a real number (i.e., need not be an integer). In that
# case, the scenario is simulated once with the floor and once with the ceiling of NperSegment, and the weighted average
# of those two results is returned. Weighting based on fractional part of original NPerSegment. For example, if nPerSegment
# is 5.14, simulate with 5 and 6, then take their weighted average with weights of .86 and .14 respectively.

run_any_subjects_scenario <- function(segmented_researcher, stat_procedure, true_effect)
{
  # Get number of subjects to use
  n_per_segment <- segmented_researcher$n_per_segment

  # If n is an integer, simply run the integer scenarion
  if (n_per_segment %% 1 == 0)
    result <- run_int_subjects_scenario(segmented_researcher, stat_procedure, true_effect)
  else
  {
    # If n has a fractional part, do the adjustment described above

    # Get the integer bounds
    lower_n <- floor(n_per_segment)
    upper_n <- ceiling(n_per_segment)

    # Assign the fractional part and 1- fractional part as probabilities for the weighted averaging
    pr_upper_n = n_per_segment - lower_n;
    pr_lower_n = 1 - pr_upper_n;

    # We will run the scenario twice -- once with n = lower and once with n = upper. In each case, we
    # need to pass the researcher into the run method, so we need to be able to change its n_per_segment
    # value. To avoid risk of destroying it in some way, we'll make a copy here.

    researcher_data_values <- unclass(segmented_researcher)
    # SegmentedResearcher <- function(max_n_segments, n_per_segment, alpha_total, alpha_strong = 0, alpha_weak = NULL)
    temp_segmented_researcher <- SegmentedResearcher(researcher_data_values$max_n_segments,
                                                     researcher_data_values$n_per_segment,
                                                     researcher_data_values$alpha_total,
                                                     researcher_data_values$alpha_strong)

    # run the lower n
    temp_segmented_researcher$n_per_segment <- lower_n
    lower_result <- run_int_subjects_scenario(temp_segmented_researcher, stat_procedure, true_effect)

    # run the upper n
    temp_segmented_researcher$n_per_segment <- upper_n
    upper_result <- run_int_subjects_scenario(temp_segmented_researcher, stat_procedure, true_effect)

    # Compute the weighted average of the two results. See method definition below.
    result <- weighted_avg_hyp_test_result(list(lower_result, upper_result), list(pr_lower_n, pr_upper_n))

  } # end else
  return(result)
} # end run_any_subjects_scenario

#========================================================================================================
# Accepts a list of SegmentedHypTestResults and a parallel list of probabilities.
# Returns a SegmentedHypTestResult instance whose values are the weighted (by prob)
# average of the input objects.
#========================================================================================================
weighted_avg_hyp_test_result <- function(results_list, probs_list)
{
  if (length(results_list) != length(probs_list)){
    stop("Results and probabilities lists must be the same length")
  }

  # We can pull representative values from any element of results_list. We use the 1st one, as per convention

  # Get the number of segments
  max_n_segments <- results_list[[1]]$max_n_segments

  # Get the names of all the fields, for associative indexing
  field_list <- names(results_list[[1]])

  # Fields 1 to 10 ields are scalars. we can process them in a single loop.
  # The remaining fields are vectors, which are processed separately

  n_scalars = 10
  scalar_list <- field_list[1:n_scalars]


  # Prepare an empty instance to hold the summary values. Missing args all default to 0
  weighted_avg_results <- SegmentedHypTestResult(max_n_segments)

  # iterate over the scalar fields, taking the weighted mean of each value, and inserting into the
  # prepared results object
  for (scalar_field_name in scalar_list)
  {

    # Pull the values for this scalar field into a list
    field_values <- sapply(results_list, function(result) result[scalar_field_name])

    # stats::weighted.mean wants vectors, not lists, so we need to unlist here.
    weighted_avg <- stats::weighted.mean(unlist(field_values), unlist(probs_list))

    # Store in skeleton
    weighted_avg_results[scalar_field_name] <- weighted_avg
  }

  # Iterate over vector fields
  # Each result has three vector fields for the  probabilities of rejecting, failing to reject and
  # continuing. Each of these vectors has one element for each segment.
  # Here we want to take a weighted average of those values, using the probability of each result.
  # So, for example, assume that we have 5 results. Let kmax = 3, so each of the five result objects
  # has a 3-element vector of pr(reject) for segments 1, 2 and 3, and also has a probability of occurrence.
  # We take the pr(reject vector) of result 1 and multiple each element vector-wise by the probability of
  # occurence of result 1. We then take result 2 and multiply its vector elements by its probability, and
  # WE SUM THE VALUES FROM 1 AND 2 VECTOR-WISE, obtaining a new vector, still of length three. We continue
  # over the remaining results, accruing the sum for each. The resulting 3-element vector which contains the
  # sum of five weighted values is the result to return.

  n_results <- length(results_list) # How many results have been passed in

  probs_vector <- unlist(probs_list)
  for (i in 1:n_results)
  {
    wght_reject <- results_list[[i]]$pr_reject_by_segment * probs_vector[i]
    weighted_avg_results$pr_reject_by_segment =  weighted_avg_results$pr_reject_by_segment + wght_reject

    wght_ftr <- results_list[[i]]$pr_ftr_by_segment * probs_vector[i]
    weighted_avg_results$pr_ftr_by_segment = weighted_avg_results$pr_ftr_by_segment + wght_ftr

    wght_cont <- results_list[[i]]$pr_continue_by_segment * probs_vector[i]
    weighted_avg_results$pr_continue_by_segment = weighted_avg_results$pr_continue_by_segment + wght_cont
  }

  return(weighted_avg_results)

} # end weighted_avg_hyp_test_result


#########################################################################################################
# Generics for public-facing methods
#########################################################################################################
#-----------------
# Roxygen comments

#' run_scenario.SegmentedHypTestEngine
#'
#' \code{run_scenario.SegmentedHypTestEngine} numerically determines expected
#' outcomes for a complete study scenario. It wraps various internal methods of
#' class SegmentedHypTestEngine (developers can cf. source code).
#' @param segmented_hyp_test_engine A SegmentedHypTestEngine instance for method
#'   dispatch
#' @param segmented_researcher Correctly initialised instance of
#'   SegmentedResearcher
#' @param stat_procedure Correctly initialised instance of any class descended
#'   from StatProcedureBase, e.g. OneSampleT, PearsonR, etc.
#' @param effects_list A list of TrueEffect instances. This must be a list, not
#'   simply a vector. Cast with list() if necessary. Each TrueEffect instance
#'   holds fields effect_size and effect_size_probability. A TrueEffect with
#'   effect_size 0 (i.e. H0 is true) may be included. Probabilities across all effect
#'   sizes must sum to 1 or an exception is thrown.
#' @return A SegmentedHypTestResult instance which holds the average expected outcomes
#' across all effect sizes, weighted by associated probability.
#'
#' @export
run_scenario <- function(segmented_hyp_test_engine,
                         segmented_researcher,
                         stat_procedure,
                         effects_list)
{
  UseMethod("run_scenario", segmented_hyp_test_engine)
}


#-----------------
# Roxygen comments

#' average_power.SegmentedHypTestEngine
#'
#' \code{average_power.SegmentedHypTestEngine} computes the average power across
#' a set of studies as defined by the effects_list argument. Needs to run from
#' the raw data, not SegHTResults, as SegHTResults do not carry their original
#' true effect size(s). Called by user-facing methods in SegmentedHypTester.R.
#' See source code for logic.
#' @param segmented_hyp_test_engine A SegmentedHypTestEngine instance for method
#'   dispatch
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param n_per_segment Number of subjects to be tested in each segment
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   experiment
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, reject and stop.
#' @param stat_procedure Character code for hypothesis testing procedure to
#'   be used in the experiment. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param effects_list List of TrueEffect instances

#' @return Average power across all studies based on effects in effect_list
#'
#' @export
average_power <- function(segmented_hyp_test_engine,
                          max_n_segments,
                          n_per_segment,
                          alpha_total,
                          alpha_strong,
                          stat_procedure,
                          effects_list)
{
  UseMethod("average_power", segmented_hyp_test_engine)
}
