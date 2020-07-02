# Class TwoSampleT descends from StatProcedureBase

# See class OneSampleT (StatProcedureOneSampleT.R) for detailed explanation of method logic

# All the computations ported from original source code (Miller, 2020)

# NB: All sample n's provided to these methods must be a vector of length 2. Passing
# a scalar or vector of length one throws an exception

#' @importFrom stats qt
#' @importFrom stats pt
#' @importFrom stats dt
#' @importFrom pracma integral
#' @importFrom utils capture.output

########################################
# Creation methods - following Wickham
########################################
new_TwoSampleT <- function(min_sample_n, display_abbrev, display_name)
{
  new_two_sample_t <- list(min_sample_n = min_sample_n,
                           display_abbrev = display_abbrev,
                           display_name = display_name)


  attr(new_two_sample_t, "class") <- c("TwoSampleT", "StatProcedureBase")

  return(new_two_sample_t)
}

#====================================================================
validate_TwoSampleT <- function(two_sample_t)
{
  # Do any validation here. Use stop(msg) to exit if things go wrong
}

#====================================================================
#-----------------
# Roxygen comments

#' TwoSampleT Constructor
#'
#' TwoSampleT is a child of StatProcedureBase, representing a two sample t-test
#'
#' \code{TwoSampleT}  is the public-facing TwoSampleT constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated TwoSampleT instance
#'
#' @export
TwoSampleT <- function(min_sample_n = 4, display_abbrev = "2t", display_name = '2-sample t-test')
{
  instance <- new_TwoSampleT(min_sample_n, display_abbrev, display_name)
  validate_TwoSampleT(instance)
  return(instance)
}
#################################
# Instance methods
#################################
#========================================================================================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
#========================================================================================================
#-----------------
# Roxygen comments

#' compute_power.TwoSampleT
#'
#' Implementation of compute_power method for class TwoSampleT
#'
#' \code{TwoSampleT} computes the statistical power of a described two sample
#' t-test as the probability density of t-critical in the t-distribution for a given effect size
#'
#' @param stat_procedure TwoSampleT instance for method dispatch
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param effect_size Presumed true effect size on Cohen's d scale
#' @param sample_n Group sizes. Must be a numeric vector of length 2
#' @param critical_value Optional to speed computation. If omitted, this value is computed in the method
#'
#' @return Computed statistical power (probability of rejecting H0 when the effect is present)
#'
#' @export
compute_power.TwoSampleT <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")

  # degrees of freedom for a two sample t-test is n1 + n2 - 2
  df <- sample_n[1] + sample_n[2] - 2

  t_crit_cumul_prob <- 1 - alpha_one_tailed

  if (is.null(critical_value))
  {
    t_crit <- qt(t_crit_cumul_prob, df)
  }
  else
  {
    t_crit <- critical_value
  }

  noncentrality_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  power <- 1 - pt(t_crit, df, ncp = noncentrality_parameter)

  return(power)
} # end compute_power.TwoSampleT


#========================================================================================================
# p_level(two_sample_t, observed_stat, sample_n)
#========================================================================================================
#' @export
p_level.TwoSampleT <- function(stat_procedure, observed_stat, sample_n)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")

  # degrees of freedom for a two sample t-test is n1 + n2 - 2
  df <- sample_n[1] + sample_n[2] - 2

  # the observed sample exceeds this proportion of t when Ho is true
  cumul_prob_obs <- pt(observed_stat, df)

  # by definition of p-value
  p <- 1 - cumul_prob_obs

  return(p)

} # end p_level

#========================================================================================================
# critical_value(two_sample_t, alpha_one_tailed, sample_n)
#========================================================================================================
#' @export
critical_value.TwoSampleT <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  #print(stat_procedure)

  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")

  # degrees of freedom for a two sample t-test is n1 + n2 - 2
  df <- sample_n[1] + sample_n[2] - 2

  t_crit_cumul_prob <- 1 - alpha_one_tailed

  t_critical <- qt(t_crit_cumul_prob, df)

  return(t_critical)

} # end critical_value

#========================================================================================================
# noncentrality(sample_n, effect_size)
#========================================================================================================
#' @export
noncentrality.TwoSampleT <- function(stat_procedure, sample_n, effect_size)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # By definition
  noncen <- sqrt((sample_n[1] * sample_n[2]) / (sample_n[1] + sample_n[2])) * effect_size
  return(noncen)

} # end noncentrality

#========================================================================================================
# effect_size_from_stat(stat_value, sample_n)
#========================================================================================================
#' @export
effect_size_from_stat.TwoSampleT <- function(stat_procedure, stat_value, sample_n)
{

  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # By definition
  effect_size <- stat_value * sqrt( (sample_n[1] + sample_n[2])/ (sample_n[1] * sample_n[2]))
  return(effect_size)

} # end effect_size_from_stat

#========================================================================================================
# expected_significant_effect
#========================================================================================================
#' @export
expected_significant_effect.TwoSampleT <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{

  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")

    t_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
    noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

    lower_bound <- t_critical
    upper_bound <- Inf

    df <- sample_n[1] + sample_n[2] - 2
    total_area_above_crit <- 1 - pt(t_critical, df, noncen_parameter)


    # The function to integrate over. The variables df and noncen_parameter must be
    # in scope and initialised at the time of the call to integral(args) below
    t_prob_times_value <- function(t)
    {
      dt(t, df = df, ncp = noncen_parameter) * t
    }

    # Integral produces a message about algorithm used. We don't want it, so we wrap the
    # call in capture.output
    capture.output(expected_value_numerator <- pracma::integral(t_prob_times_value, lower_bound, upper_bound))
    expected_t <- expected_value_numerator/total_area_above_crit

    expected_effect_size <- effect_size_from_stat(stat_procedure, expected_t, sample_n)

    return(expected_effect_size)
} # end expected_significant_effect

#========================================================================================================
#expected_significant_effect_bounded
#========================================================================================================
#' @export
expected_significant_effect_bounded.TwoSampleT <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  t_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  t_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- t_critical_weak
  upper_bound <- t_critical_strong

  df <- sample_n[1] + sample_n[2] - 2
  total_area_above_strong <- 1 - pt(t_critical_strong, df, noncen_parameter)
  total_area_above_weak <- 1 - pt(t_critical_weak, df, noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong


  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  t_prob_times_value <- function(t)
  {
    dt(t, df = df, ncp = noncen_parameter) * t
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(t_prob_times_value, lower_bound, upper_bound))
  expected_t <- expected_value_numerator/total_area_between

  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_t, sample_n)

  return(expected_effect_size)


} # end expected_significant_effect_bounded

#========================================================================================================
# divide_total_sample
#========================================================================================================
#' @export
divide_total_sample.TwoSampleT <- function(stat_procedure, n_per_segment)
{

  # split total subject into equal groups (or as close as possible, if total n is odd)
  group_sizes <- numeric(2)
  group_sizes[1] <- floor(n_per_segment/2)
  group_sizes[2]  <- n_per_segment - group_sizes[1]

  return(group_sizes)
} # end divid_total_sample

