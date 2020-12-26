# Class OneSampleZ descends from StatProcedureBase

# See class OneSampleT for detailed explanations of method logic

# All the computations ported from original source code (Miller, 2020)

# Note that the sampling distribution for z, when Ho is true, is a normal with mean = 0 and s = 1.
# Noncentrality shifts the distribution mean depending on n and the true effect size

# Use
# dnorm(x, mean = 0, sd = 1, log = FALSE)
# pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
# qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats dnorm
#' @importFrom pracma integral
#' @importFrom utils capture.output

########################################
# Creation methods - following Wickham
########################################
#====================================================================
new_OneSampleZ <- function(min_sample_n, display_abbrev, display_name)
{
  new_one_sample_z <- list(min_sample_n = min_sample_n,
                           display_abbrev = display_abbrev,
                           display_name = display_name)


  attr(new_one_sample_z, "class") <- c("OneSampleZ", "StatProcedureBase")

  return(new_one_sample_z)
}
#====================================================================
validate_OneSampleZ <- function(one_sample_z)
{
  # Do any validation here. Use stop(msg) to exit if things go wrong
}
#====================================================================
#-----------------
# Roxygen comments

#' OneSampleZ Constructor
#'
#' OneSampleZ is a child of StatProcedureBase, representing a one-sample z-test
#'
#' \code{OneSampleZ}  is the public-facing OneSampleZ constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated OneSampleZ instance
#'
#' @export
OneSampleZ <- function(min_sample_n = 2, display_abbrev = "1z", display_name = '1-sample z-test')
{
  instance <- new_OneSampleZ(min_sample_n, display_abbrev, display_name)
  validate_OneSampleZ(instance)
  return(instance)
}


#################################
# Instance methods
#################################
#========================================================================================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
#========================================================================================================
#' @export
compute_power.OneSampleZ <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  z_crit_cumul_prob = 1 - alpha_one_tailed

  if (is.null(critical_value))
  {
    z_crit = qnorm(z_crit_cumul_prob)
  }
  else
  {
    z_crit = critical_value
  }

  noncentrality_parameter = noncentrality(stat_procedure, sample_n, effect_size)

  # When Ho is false, z-critical remains the same, but the true z-distribution shifts upward by the effect
  # size and sample size as reflected by the noncentrality parameter. The area that causes rejection is that
  # above z-critical in this new distribution, i.e. that above z-critical - the noncen parameter, reflecting the
  # shift.
  power = 1 - pnorm(z_crit - noncentrality_parameter)

  return(power)

} # end compute_power


#========================================================================================================
# p_level(one_sample_z, observed_stat, sample_n)
#========================================================================================================
#' @export
p_level.OneSampleZ <- function(stat_procedure, observed_stat, sample_n)
{

  # the observed sample exceeds this proportion of the normal when Ho is true
  cumul_prob_obs <- pnorm(observed_stat)

  # by definition of p-value
  p <- 1 - cumul_prob_obs

  return(p)

} # end p_level

#========================================================================================================
# critical_value(one_sample_z, alpha_one_tailed, sample_n)
#========================================================================================================
#' @export
critical_value.OneSampleZ <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  #print(stat_procedure)


  # To obtain your critical value, find the z at 1-alpha percentile in the Ho: n(0,1) distribution
  # e.g. if alpha = 0.05, you want the z that cuts off the lower 95% of the Ho distribution.

  z_crit_cumul_prob = 1 - alpha_one_tailed

  z_critical = qnorm(z_crit_cumul_prob)

  return(z_critical)

} # end critical_value

#========================================================================================================
# noncentrality(sample_n, effect_size)
#========================================================================================================
#' @export
noncentrality.OneSampleZ <- function(stat_procedure, sample_n, effect_size)
{
  # By definition
  noncen = sqrt(sample_n) * effect_size
  return(noncen)

} # end noncentrality

#========================================================================================================
# effect_size_from_stat(stat_value, sample_n)
#========================================================================================================
#' @export
effect_size_from_stat.OneSampleZ <- function(stat_procedure, stat_value, sample_n)
{
  # By definition
  effect_size = stat_value / sqrt(sample_n);
  return(effect_size)

} # end effect_size_from_stat

#========================================================================================================
# expected_significant_effect - See StatProcedureOneSampleT.R for explanation of logic
#========================================================================================================
#' @export
expected_significant_effect.OneSampleZ <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{

  # Fetch the appropriate critical value and noncentrality parameter as defined for this statistic
  z_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  # We will integrate over the portion of the noncentral sampling distribution greater than z-critical
  lower_bound <- z_critical
  upper_bound <- Inf

  # Required to make a proportional adjustment in the computation of expected effect size
  total_area_above_crit <- 1 - pnorm(z_critical - noncen_parameter)

  # The function to integrate over. The variable noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  z_prob_times_value <- function(z)
  {
    dnorm(z - noncen_parameter) * z
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(z_prob_times_value, lower_bound, upper_bound))
  expected_z <- expected_value_numerator/total_area_above_crit

  expected_effect_size = effect_size_from_stat(stat_procedure, expected_z, sample_n)

  return(expected_effect_size)

} # end expected_significant_effect

#========================================================================================================
# expected_significant_effect_bounded - See StatProcedureOneSampleT.R for explanation of logic
#========================================================================================================
#' @export
expected_significant_effect_bounded.OneSampleZ <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  z_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  z_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- z_critical_weak
  upper_bound <- z_critical_strong

  total_area_above_strong <- 1 - pnorm(z_critical_strong - noncen_parameter)
  total_area_above_weak <- 1 - pnorm(z_critical_weak - noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong

  z_prob_times_value <- function(z)
  {
    dnorm(z - noncen_parameter) * z
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(z_prob_times_value, lower_bound, upper_bound))
  expected_z <- expected_value_numerator/total_area_between

  expected_effect_size = effect_size_from_stat(stat_procedure, expected_z, sample_n)

  return(expected_effect_size)

} # end expected_significant_effect_bounded


#========================================================================================================
# divide_total_sample
#========================================================================================================
#' @export
divide_total_sample.OneSampleZ <- function(stat_procedure, n_per_segment)
{
  # In a one-sample z-test, all subjects are in a single group
  return(n_per_segment)
} # end divid_total_sample
