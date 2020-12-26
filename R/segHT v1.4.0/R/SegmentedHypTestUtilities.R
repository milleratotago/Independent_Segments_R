# For maximum ease of access, these end-user facing methods are not contained
# in a class, so they need no generics or class instance for dispatch.
# Top level methods for non-programming users.

#' @importFrom pracma fzero
#' @importFrom pracma fminbnd


# Ignore this 13-03-20
placeholder <- function()
{
}

###########################################################################################
# alpha_weak: User-facing wrapper for alpha_weak_computation
###########################################################################################
#-----------------
# Roxygen comments

#' alpha_weak
#'
#' \code{alpha_weak}  returns the value of alpha_weak which
#' will produce the correct overall Type 1 error probability (i.e. alpha_total)
#' for a segmented hypothesis testing design, given values for maximum
#' number of segments, alpha total and alpha strong. This function wraps the
#' alpha_weak_computation method of the SegmentedResearcher class. If alpha
#' strong has been specified as zero (boundary case) alpha weak has an
#' analytical solution alpha_total^(1/max_n_segments). Otherwise alpha weak is
#' found via numerical search using function pracma::fzero. See Miller and
#' Ulrich, 2020 for discussion.
#'
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   study
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, stop and reject Ho.
#'
#' @return Numeric. Value of alpha weak.
#'
#' @export
alpha_weak <- function(max_n_segments, alpha_total, alpha_strong)
{
  # We create a SegmentedResearcher to call its alpha_weak method.
  segmented_researcher <- SegmentedResearcher()
  alpha_weak <- alpha_weak_computation(segmented_researcher, max_n_segments, alpha_total, alpha_strong)

  return(alpha_weak)
}

###########################################################################################
# n_for_power: User-facing wrapper for n_for_power_computation
###########################################################################################

#-----------------
# Roxygen comments

#' n_for_power
#'
#' \code{n_for_power}  returns the number of subjects
#' required to achieve a target statistical power for a segmented
#' hypothesis testing design, given the other design descriptors (see Arguments,
#' below). This method wraps n_for_power_computation to prevent the end user from
#' having to create his/her own component class instances. The result is
#' determined by numerical search using pracma::fminbnd. The function minimised
#' computes the difference between the power obtained by a given value of n and
#' the target power. The power obtained for any given value of n is computed
#' numerically. See n_for_power_computation in this file and
#' seght::run_scenario.SegmentedHypTestEngine for computation logic.
#'
#' @param target_power Desired power level (probability of detecting a false H0
#'   when present) across all segments of the study
#' @param stat_procedure_name Character code for hypothesis testing procedure to
#'   be used in the study. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param effect_size Size of the true effect when it is present.
#'   For tests of correlation, this is the true correlation rho.
#'   For tests of means (1t, 2t, 1z, 2z), this Cohen's d.
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   study
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha_strong, stop and reject Ho.
#' @param alpha_weak For decision rule: If p-observed for a segment is greater
#'   than alpha weak, stop and fail to reject Ho. This argument is not normally provided,
#'   in which case the correct value is computed in the method. If the argument is
#'   provided, it is used without being checked (for speed) so it is the user's
#'   responsibility to ensure that it is the correct value.
#'
#' @return Returns the number
#'   of subjects per segment that yields the target power. The result is not
#'   necessarily a whole number, in which case you must use the next
#'   smaller or larger integer and get a little lower or higher power.
#' @export
n_for_power <- function(target_power,
                        stat_procedure_name,
                        effect_size,
                        max_n_segments = 3,
                        alpha_total = 0.5,
                        alpha_strong = 0.025,
                        alpha_weak = NA)
{

  # Assign a generous starting guess for n_per_segment. Computation code will raise
  # this if necessary to gain the required power level.
  n_per_segment <- 50

  # Build necessary classes for computation
  segmented_researcher <- SegmentedResearcher(max_n_segments, n_per_segment, alpha_total, alpha_strong, alpha_weak)
  stat_procedure <- stat_procedure_factory(stat_procedure_name)

  # Call the computation method
  # This method returns a three-item vector [n per segment, exact power, computed alpha weak]
  n_for_power_outcomes <- n_for_power_computation(target_power,
                                        segmented_researcher,
                                        stat_procedure,
                                        effect_size)

  n_for_power <- n_for_power_outcomes[1] # First item of the three returned by no_for_power_computation

  # return the result.
  return(n_for_power)

} # n_for_power wrapper function


###########################################################################################
# segmented_hyp_test_outcomes: User-facing wrapper for SegmentedHypTestEngine scenario methods.
#
# User supplies scalar descriptors for an experimental scenario (kmax, total alpha, alpha strong,
# statistical test, total n per segment, true effect size(s) and base rate(s), and receives
# a summary of the expected experimental outcomes (alpha weak, average power, expected number of subjects used,
# pr(stop & reject) and pr(stop and fail to reject) at each segment.
###########################################################################################
#-----------------
# Roxygen comments

#' segmented_hyp_test_outcomes
#'
#' \code{segmented_hyp_test_outcomes}  wraps
#' seght::run_scenario.SegmentedHypTestEngine to numerically determine the
#' probabilities of all outcomes (reject, fail to reject, or continue) for each
#' potential segment of a segmented hypothesis testing design, along with a
#' variety of useful summary descriptors of expected study performance.
#'
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param n_per_segment Number of subjects to be tested in each segment
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   experiment
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, stop and reject Ho.
#' @param stat_procedure_name Character code for hypothesis testing procedure to
#'   be used in the experiment. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param effect_sizes Value, or vector of values, representing the true effect when it is present.
#'   For tests of correlation, this is the true correlation rho.
#'   For tests of means (1t, 2t, 1z, 2z), this Cohen's d.
#'   Any effect size of zero is logically Ho = true.
#' @param base_rates Base rate or vector of base rates giving the probability of each effect size.
#'   Each base_rate value matches the corresponding element of the effect_sizes argument, above.
#'   Each base rate must be between 0 and 1, and their total should be less than or equal to 1.
#'   If the total of base_rates is less than 1, the
#'   remaining probability is allocated to effect size 0 (i.e. H0 is true)
#'
#' @return Returns a list (collection of named items) containing all of the method
#'   input values (for reference), plus: \itemize{ \item exp_n_sub: Expected (mean
#'   across infinite iterations) number of subjects used in the study \item
#'   avg_power: Power level of the study on average across all cases where Ho is false
#'   (weighted average according to base rates).
#'   \item pr_reject_by_segment: Probablility of rejecting H0 at each segment
#'   \item pr_ftr_by_segment: Probability of failing to reject H0 at each
#'   segment }
#'
#' @export
segmented_hyp_test_outcomes <- function(max_n_segments,
                                        n_per_segment,
                                        alpha_total,
                                        alpha_strong,
                                        stat_procedure_name,
                                        effect_sizes,
                                        base_rates = 1)
{

  # We are aiming for a call to run_scenario.SegmentedHypTestEngine. We prepare our arguments first.
  segmented_researcher <- SegmentedResearcher(max_n_segments, n_per_segment, alpha_total, alpha_strong)
  stat_procedure <- stat_procedure_factory(stat_procedure_name)

  # run_scenario takes a list of true effects whose probabilities sum to 1. If the user has passed
  # in a baserate or vector of baserates that do not sum to 1, we add an element with
  # effect size 0 (i.e. H0 true) with base rate as required to get the total to 1

  if (length(base_rates) != length(effect_sizes)) {
    stop("must have equal numbers of effect_sizes and base_rates")
  }

  total_effect_prob <- sum(base_rates)

  if (total_effect_prob > 1) {
    stop("sum of base_rates cannot exceed 1.0")
  }

  if (total_effect_prob < 1) {
    ho_base_rate <- 1 - total_effect_prob
    ho_effect_size <- 0
    effect_sizes <- c(effect_sizes, ho_effect_size)
    base_rates <- c(base_rates, ho_base_rate)
  }

  # Create TrueEffect objects for each effect_size/base_rate pair
  effects_list <- list(length(effect_sizes))
  for (i in 1:length(effect_sizes))
  {
    te <- TrueEffect(effect_sizes[i], base_rates[i])
    effects_list[[i]] <- te
  }

  # Make the call
  engine <- SegmentedHypTestEngine()
  experiment_results <- run_scenario(engine, segmented_researcher, stat_procedure, effects_list)

  # Compute the average power
  avg_power <-average_power(engine,
                            max_n_segments,
                            n_per_segment,
                            alpha_total,
                            alpha_strong,
                            stat_procedure,
                            effects_list)


  # Pull other necessary values to build the return list
  alpha_weak <- segmented_researcher$alpha_weak # This value initialised in the SegmentedResearcher ctor

  exp_n_subj <- experiment_results$exp_n_subj
  exp_n_subj_sqr = experiment_results$exp_n_subj_sqr
  var_exp_n_subj <- exp_n_subj_sqr - (exp_n_subj^2)
  sd_exp_n_subj <- sqrt(var_exp_n_subj)

  pr_reject_by_segment = experiment_results$pr_reject_by_segment
  pr_ftr_by_segment = experiment_results$pr_ftr_by_segment

  # Build it
  output_value <- list(max_n_segments = max_n_segments,
                       n_per_segment = n_per_segment,
                       alpha_total = alpha_total,
                       alpha_strong = alpha_strong,
                       alpha_weak = alpha_weak,
                       stat_procedure = stat_procedure_name,
                       effects_list = effects_list,
                       exp_n_subj = exp_n_subj,
                       sd_exp_n_subj = sd_exp_n_subj,
                       avg_power = avg_power,
                       pr_reject_by_segment = pr_reject_by_segment,
                       pr_ftr_by_segment = pr_ftr_by_segment)

  # Return it
  return(output_value)
}


####################################################################################
# search_kmax: Accepts total alpha, alpha strong, a test name, an effect size with optional
# effect probability (baserate; defaults to 1 if not provided) and a target power level. Returns
# a dataframe with one row for each value of kmax between 2 and 10 (a practical maximum).
# Columns indicate the appropriate alpha weak, the n per segment required to achieve the target power,
# and the expected number of subjects and expected number of segments for the study
####################################################################################

#-----------------
# Roxygen comments

#' search_kmax
#'
#' \code{search_kmax} computes a set  of performance
#' descriptors (see Values, below) for a specified study scenario for each
#' value of kmax (maximum number of segments) between 2 and 10. Returned values
#' are based on the n_per_segment required to achieve the provided target statistical power.
#'
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   study
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, stop and reject Ho.
#' @param stat_procedure_name Character code for hypothesis testing procedure to
#'   be used in the experiment. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param target_power Desired power level (probability of detecting a false H0
#'   when present) across all segments of the study
#' @param effect_size Size of the true effect when it is present.
#'   For tests of correlation, this is the true correlation rho.
#'   For tests of means (1t, 2t, 1z, 2z), this Cohen's d.
#' @param base_rate Base rate of the true effect_size. Must be between 0 and 1.
#'   If this value is less than 1, the remaining probability is allocated to
#'   effect size 0 (i.e. H0 is true)
#' @return Returns a data frame with one row for each value of kmax between 2
#'   and 10, and columns: \itemize{ \item kmax: Maximum number of segments in
#'   the study \item alpha_weak: For decision rule: If p-observed for a segment
#'   is greater than alpha weak, stop and fail to reject Ho. \item n per segment:
#'   Number of subjects per segment needed to achieve the target power. \item expected n subjects: Expected (mean)
#'   number of subjects required across infinite iterations of the experiment
#'   \item expected n segments: Expected (mean) number of segments that will be
#'   performed prior to termination (with any outcome) across infinite
#'   iterations of the experiment. }
#'
#' @export
search_kmax <- function(alpha_total,
                        alpha_strong,
                        stat_procedure_name,
                        target_power,
                        effect_size,
                        base_rate = 1)
{
  # Create engine instance to run scenarios
  engine <- SegmentedHypTestEngine()

  # Prepare skeleton vectors to hold outputs
  kmax_holder <- c()
  alpha_weak_holder <- c()
  n_per_segment_holder <- c()
  exp_n_subj_holder <- c()
  exp_n_segments_holder <- c()

  # Prepare required args
  stat_procedure <- stat_procedure_factory(stat_procedure_name)

  n_segment_upper_limit = 10 # By definition
  n_subjects_default = 50

  # Loop over potential values of kmax

  # For each potential value of kmax, we call n_for_power_computation (internal method of this class)
  # to find an appropriate n_per_segment to get the desired power.
  # We then run the scenario with that value to get expected n subjects and expected n segments.

  # fzero throws if max_n_segments is < 2
  for (max_n_segments in 2:n_segment_upper_limit)
  {
    # SegmentedResearcher$alpha_weak is initialised in the ctor
    segmented_researcher <- SegmentedResearcher(max_n_segments,
                                                n_subjects_default,
                                                alpha_total,
                                                alpha_strong)

    # n_for_power returns a vector [n, power]
    n_per_segment_output <- n_for_power_computation(target_power,
                                                    segmented_researcher,
                                                    stat_procedure,
                                                    effect_size)
    n_per_segment <- n_per_segment_output[1]

    # We wish to call run_scenario with this n_per_segment.
    # run_scenario takes a list of true effects whose probabilities sum to 1. If the user has passed
    # in a baserate, we create the appropriate TrueEffect instance for effect size 0 (i.e. H0 true)

    real_effect <- TrueEffect(effect_size, base_rate)
    if (base_rate == 1) {
      effects_list <- list(real_effect)
    } else {
      null_effect_probability <- 1 - base_rate
      null_effect <- TrueEffect(0, null_effect_probability)
      effects_list <- list(real_effect, null_effect)
    }

    # Remake the researcher with the n_per_segment computed by find_n_for_power
    segmented_researcher <- SegmentedResearcher(max_n_segments,
                                                n_per_segment,
                                                alpha_total,
                                                alpha_strong)

    experiment_results <- run_scenario(engine, segmented_researcher, stat_procedure, effects_list)

    # Expected number of segments is (expected total number of subjects)/(number of subjects per segment)
    exp_n_subjects <- experiment_results$exp_n_subj
    exp_n_segments <- exp_n_subjects/n_per_segment

    # Gather up the output values
    kmax_holder <- c(kmax_holder, max_n_segments)
    alpha_weak_holder <- c(alpha_weak_holder, segmented_researcher$alpha_weak)
    n_per_segment_holder <- c(n_per_segment_holder, n_per_segment)
    exp_n_subj_holder <- c(exp_n_subj_holder, exp_n_subjects )
    exp_n_segments_holder <- c(exp_n_segments_holder, exp_n_segments)

  } # end for all values of kmax

  # Build the output dataframe
  output_df <- data.frame(kmax_holder,
                          alpha_weak_holder,
                          n_per_segment_holder,
                          exp_n_subj_holder,
                          exp_n_segments_holder)

  colnames(output_df) <- c("kmax", "alpha_weak", "n_per_segment", "exp_n_subjects", "exp_n_segments")

  # Return it
  return(output_df)
}



#####################
# Internal methods
#####################
#========================================================================================================
# stat_procedure_factory: Accepts string from {'1t', '2t', '1z', '2z', 'r'}. Returns initialised objects
# instance from StatProcedure family.
#========================================================================================================
stat_procedure_factory <- function(stat_procedure_name)
{
  stat_names <- c("1t", "2t", "1z", "2z", "r")

  if (!any(stat_names == stat_procedure_name))
  {
    stop("stat_procedure_name must be one of {'1t', '2t', '1z', '2z', 'r'} ")
  }

  switch(stat_procedure_name,
         "1t" = stat_procedure <- OneSampleT(),
         "2t" = stat_procedure <- TwoSampleT(),
         "1z" = stat_procedure <- OneSampleZ(),
         "2z" = stat_procedure <- TwoSampleZ(),
         "r"  = stat_procedure <- PearsonR())

  return(stat_procedure)
}


#####################################################################################
# Find the N that is needed _PER SEGMENT_ to have the indicated TargetPower for the
# specified researcher, stat type and true effect
#####################################################################################
n_for_power_computation <- function(target_power, segmented_researcher, stat_procedure, effect_size, sample_max = 50)
{

  # Copy input researcher for use in computation
  researcher_data_values <- unclass(segmented_researcher)
  temp_segmented_researcher <- SegmentedResearcher(researcher_data_values$max_n_segments,
                                                   researcher_data_values$n_per_segment,
                                                   researcher_data_values$alpha_total,
                                                   researcher_data_values$alpha_strong)

  # Never a base rate in this computation, by definition
  true_effect <- TrueEffect(effect_size, effect_size_probability = 1)

  engine <- SegmentedHypTestEngine()


  # Find minimum sample size for this statistical procedure. Maximum is an input argument
  # that defaults to 50 (pers. comm. Miller, 2020)
  sample_min <- stat_procedure$min_sample_n

  # Make sure that the desired power is possible between sample min and max. If not, slide the
  # bounds up in sample_max size steps, until it is.
  found_bounds <- FALSE

  while (!found_bounds)
  {
    temp_segmented_researcher$n_per_segment <- sample_max
    temp_results <- run_scenario(engine, temp_segmented_researcher, stat_procedure, list(true_effect))
    if (temp_results$pr_reject_total < target_power) {
      sample_min <- sample_max
      sample_max <- sample_max * 2
    } else {
      found_bounds <- TRUE
    }
  } # not found_bounds

  # With your bounds in hand, use numerical search to generate the target sample size

  # Define the function for minimisation of distance between pr(reject) and target power.
  # Variables other than n must be in scope an initialised
  power_error_fn <- function(n)
  {
    temp_segmented_researcher$n_per_segment <- n
    result <- run_scenario(engine, temp_segmented_researcher, stat_procedure, list(true_effect))
    error <- abs(result$pr_reject_total - target_power)
  }

  # Here is the search. It generates messages about its algorithms, which we do not want.
  # We use capture.output to suppress these.
  capture.output(best_n <- pracma::fminbnd(power_error_fn, sample_min, sample_max, tol = 1e-2))

  # Compute the precise power for your guess
  temp_segmented_researcher$n_per_segment <- best_n$xmin
  final_result <- run_scenario(engine, temp_segmented_researcher, stat_procedure, list(true_effect))

  output_value <- c(best_n$xmin, final_result$pr_reject_total, temp_segmented_researcher$alpha_weak)
  return(output_value)

} # find N for power computation
