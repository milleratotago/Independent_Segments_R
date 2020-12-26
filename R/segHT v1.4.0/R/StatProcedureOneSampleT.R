# Class OneSampleT descends from StatProcedureBase. See StatProcedureBase.R for additional comments

# All computations ported from the original matlab source code (Miller, 2020)

#' @importFrom stats qt
#' @importFrom stats pt
#' @importFrom stats dt
#' @importFrom pracma integral
#' @importFrom utils capture.output

#################################################################################
# Creation methods - following Wickaman
#################################################################################
#=================================================================================
new_OneSampleT <- function(min_sample_n, display_abbrev, display_name)
{
  new_one_sample_t <- list(min_sample_n = min_sample_n,
                           display_abbrev = display_abbrev,
                           display_name = display_name)


  attr(new_one_sample_t, "class") <- c("OneSampleT", "StatProcedureBase")

  return(new_one_sample_t)
}

#=================================================================================
validate_OneSampleT <- function(one_sample_t)
{
  # Add validation here. Use stop(msg) to exit if things go wrong
}

#=================================================================================
#-----------------
# Roxygen comments

#' OneSampleT Constructor
#'
#' OneSampleT is a child of StatProcedureBase, representing a one-sample t-test
#'
#' \code{OneSampleT}  is the public-facing OneSampleT constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see one_sample_t_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated OneSampleT instance
#'
#' @export
OneSampleT <- function(min_sample_n = 1, display_abbrev = "1t", display_name = '1-sample t-test')
{
  instance <- new_OneSampleT(min_sample_n, display_abbrev, display_name)
  validate_OneSampleT(instance)
  return(instance)
}

#################################
# Instance methods
#################################
#========================================================================================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL): Compute power for a single
# alpha level, effect size and sample n. Power is the probability (cumulative density) of t critical in the t
# distribution for the specified df and effect_size. This function first computes that t critical, and noncentrality,
# (a measure of the shift of the t distribution away from H0) then uses R.pt to fetch the cumulative probability.
#========================================================================================================
#' @export
compute_power.OneSampleT <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  # degrees of freedom for a one sample t-test is n-1
  df <- sample_n - 1

  # t-critical is the value which cuts an alpha-sized proportion off the tail of the
  # t distribution when Ho is true. Here we compute the size of the distribution below
  # t-critical (i.e. those values which lead to failure to reject), by definition, 1 - alpha.
  # This is the cumulative probability of (i.e. proportion of scores less than) t critical in
  # the t-distribution when Ho is true

  t_crit_cumul_prob <- 1 - alpha_one_tailed

  # We now use that cumulative probability to retrieve the actual value of t-critical
  # using the R.qt function, which performs this job for the t distribution. Because the ncp
  # parameter if not set in the call to qt, we are using the t distribution when Ho is true.

  # Alternatively, you can just pass the critical value into the function.
  if (is.null(critical_value)){
    t_crit <- qt(t_crit_cumul_prob, df)
  } else {
    t_crit <- critical_value
  }

  # This will dispatch to the object instance's noncentrality function
  noncentrality_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  # The pt function accepts a value of t and returns the cumulative probability (i.e. the chances
  # that value less than this occurs). The df and noncen params adjust for sample and effect size.
  # Power (pr of rejecting when H0 is false) is 1 - that value.
  power <- 1 - pt(t_crit, df, ncp = noncentrality_parameter)

  return(power)

} # end compute_power


#========================================================================================================
# p_level(one_sample_t, observed_stat, sample_n): Computes the p-level for a supplied t-observed by
# calling R.pt, which returns cumulative density (probability) for given values.
#========================================================================================================
#' @export
p_level.OneSampleT <- function(stat_procedure, observed_stat, sample_n)
{
  print("In p_level.OneSampleT")
  # by definition for one sample t
  df <- sample_n - 1

  # the observed sample exceeds this proportion of t when Ho is true
  cumul_prob_obs <- pt(observed_stat, df)

  # by definition of p-value
  p <- 1 - cumul_prob_obs

  return(p)

} # end p_level

#========================================================================================================
# critical_value(one_sample_t, alpha_one_tailed, sample_n): Returns the t-critical value for
# a given alpha (one-tailed) and sample n. Fetches from R.qt, which returns t values for given
# cumulative probabilities.
#========================================================================================================
#' @export
critical_value.OneSampleT <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  # To obtain your critical value, find the t at 1-alpha percentile in the Ho distribution
  # e.g. if alpha = 0.05, you want the t that cuts off the lower 95% of the Ho distribution.

  # degrees of freedom for a one sample t-test is n-1
  df <- sample_n - 1

  t_crit_cumul_prob <- 1 - alpha_one_tailed

  t_critical <- qt(t_crit_cumul_prob, df)

  return(t_critical)

} # end critical_value

#========================================================================================================
# noncentrality(one_sample_t, sample_n, effect_size): Parameter to describe how far the t-distribution is shifted
# from the standard (i.e. that which exists when effect size = 0, Ho is true)
#========================================================================================================
#' @export
noncentrality.OneSampleT <- function(stat_procedure, sample_n, effect_size)
{
  # By definition of noncentrality for one sample t
  noncen <- sqrt(sample_n) * effect_size
  return(noncen)

} # end noncentrality

#========================================================================================================
# effect_size_from_stat(one_sample_t, stat_value, sample_n): By definition for one sample-t, eta = tobs/sqrt(n)
#========================================================================================================
#' @export
effect_size_from_stat.OneSampleT <- function(stat_procedure, stat_value, sample_n)
{
  # By definition of effect size for one sample t
  effect_size <- stat_value / sqrt(sample_n);
  return(effect_size)

} # end effect_size_from_stat

#========================================================================================================
# expected_significant_effect: For simulations to demonstrate bias in published observed effect sizes.
# Determines the expected (mean) effect size, if one considers only those
# cases where the null hypothesis is rejected (i.e. the observed stat is greater than the critical). For each
# test, this is computed by integrating across stat * stat_density for all values above the critical
# to obtain an expected (mean) value for that portion of the distribution, than dividing by the total area
# under the curve to obtain the correct proportional value. See child code.
#========================================================================================================
#' @export
expected_significant_effect.OneSampleT <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  # Find the conditional expectation of the observed effect size, conditional on p<Alpha1Tailed.

  # Prepare the parameters for the computation
  t_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- t_critical
  upper_bound <- Inf

  df <- sample_n - 1
  total_area_above_crit <- 1 - pt(t_critical, df, noncen_parameter)


  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  t_prob_times_value <- function(t)
  {
    suppressWarnings(
      dt(t, df = df, ncp = noncen_parameter) * t
    )
  }

  # Take the integral
  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(t_prob_times_value, lower_bound, upper_bound))

  # Normalise to the proportion of the curve above critical
  expected_t <- expected_value_numerator/total_area_above_crit

  # Determine the effect size from the resulting expected mean t
  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_t, sample_n)

  return(expected_effect_size)

} # end expected_significant_effect
#========================================================================================================
# expected_significant_effect_bounded. As above, but for a region for the curve between the two critical
# values determined by alpha strong and alpha weak.
#========================================================================================================
#' @export
expected_significant_effect_bounded.OneSampleT <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{

  t_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  t_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- t_critical_weak
  upper_bound <- t_critical_strong

  df <- sample_n - 1
  total_area_above_strong <- 1 - pt(t_critical_strong, df, noncen_parameter)
  total_area_above_weak <- 1 - pt(t_critical_weak, df, noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong


  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  t_prob_times_value <- function(t)
  {
    suppressWarnings(
      dt(t, df = df, ncp = noncen_parameter) * t
    )
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(t_prob_times_value, lower_bound, upper_bound))
  expected_t <- expected_value_numerator/total_area_between

  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_t, sample_n)

  return(expected_effect_size)


} # end expected_significant_effect

#========================================================================================================
# divide_total_sample: Experiment descriptions (SegmentedResearcher class instances) hold n *per segment*,
# summed across all groups. Each test needs to compute the number of subjects *per group*. In a one-sample
# t-test, all subjects are in a single group so this value is simply n_per_segment. See TwoSampleT for example
# where computation is required.
#========================================================================================================
#' @export
divide_total_sample.OneSampleT <- function(stat_procedure, n_per_segment)
{
  return(n_per_segment)
} # end divide_total_sample
