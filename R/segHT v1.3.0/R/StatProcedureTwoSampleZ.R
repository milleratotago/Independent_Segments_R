# Class TwoSampleZ descends from StatProcedureBase

# See class OneSampleT (StatProcedureOneSampleT.R) for detailed explanation of method logic

# All computations ported from the original source code

# NB: All sample n's provided to these methods must be a vector of length 2. Passing
# a scalar or vector of length one throws an exception

# Note that the sampling distribution for z, when Ho is true, is a normal with mean = 0 and s = 1.
# Noncentrality shifts the distribution mean depending on n and the true effect size

#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom pracma integral
#' @importFrom utils capture.output

########################################
# Creation methods
########################################

#=============================================================================
new_TwoSampleZ <- function(min_sample_n, display_abbrev, display_name)
{
  new_two_sample_z <- list(min_sample_n = min_sample_n,
                           display_abbrev = display_abbrev,
                           display_name = display_name)


  attr(new_two_sample_z, "class") <- c("TwoSampleZ", "StatProcedureBase")

  return(new_two_sample_z)
}

#=============================================================================
validate_TwoSampleZ <- function(two_sample_z)
{
  # Do any validation here. Use stop(msg) to exit if things go wrong
}

#=============================================================================
#-----------------
# Roxygen comments

#' TwoSampleZ Constructor
#'
#' TwoSampleZ is a child of StatProcedureBase, representing a two sample z-test
#'
#' \code{TwoSampleZ}  is the public-facing TwoSampleZ constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated TwoSampleZ instance
#'
#' @export
TwoSampleZ <- function(min_sample_n = 2, display_abbrev = "2z", display_name = '2-sample z-test')
{
  instance <- new_TwoSampleZ(min_sample_n, display_abbrev, display_name)
  validate_TwoSampleZ(instance)
  return(instance)
}


#################################
# Instance methods
#################################
#========================================================================================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
#========================================================================================================
#' @export
compute_power.TwoSampleZ <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{

  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  z_crit_cumul_prob <- 1 - alpha_one_tailed

  if (is.null(critical_value))
  {
    z_crit <- qnorm(z_crit_cumul_prob)
  }
  else
  {
    z_crit <- critical_value
  }


  noncentrality_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  # When Ho is false, z-critical remains the same, but the true z-distribution shifts upward by the effect
  # size and sample size as reflected by the noncentrality parameter. The area that causes rejection is that
  # above z-critical in this new distribution, i.e. that above z-critical - the noncen parameter, reflecting the
  # shift.
  power <- 1 - pnorm(z_crit - noncentrality_parameter)

  return(power)

} # end compute_power

#========================================================================================================
# p_level(two_sample_z, observed_stat, sample_n)
#========================================================================================================
#' @export
p_level.TwoSampleZ <- function(stat_procedure, observed_stat, sample_n)
{

  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # the observed sample exceeds this proportion of the normal when Ho is true
  cumul_prob_obs <- pnorm(observed_stat)

  # by definition of p-value
  p <- 1 - cumul_prob_obs

  return(p)

} # end p_level

#========================================================================================================
# critical_value(two_sample_z, alpha_one_tailed, sample_n): Returns the z-critical value for
# a given alpha (one-tailed) and sample n. Fetches from R.qnorm
#========================================================================================================
#' @export
critical_value.TwoSampleZ <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  #print(stat_procedure)


  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # To obtain your critical value, find the z at 1-alpha percentile in the Ho (n(0,1) distribution
  # e.g. if alpha = 0.05, you want the z that cuts off the lower 95% of the Ho distribution.

  z_crit_cumul_prob <- 1 - alpha_one_tailed

  z_critical <- qnorm(z_crit_cumul_prob)

  return(z_critical)

} # end critical_value

#========================================================================================================
# noncentrality(sample_n, effect_size)
#========================================================================================================
#' @export
noncentrality.TwoSampleZ <- function(stat_procedure, sample_n, effect_size)
{
  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # By definition
  sum_ns <- sample_n[1] + sample_n[2]
  product_ns <- sample_n[1] * sample_n[2]
  noncen <- sqrt(product_ns/sum_ns) * effect_size
  return(noncen)

} # end noncentrality

#========================================================================================================
# effect_size_from_stat(stat_value, sample_n)
#========================================================================================================
#' @export
effect_size_from_stat.TwoSampleZ <- function(stat_procedure, stat_value, sample_n)
{

  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # By definition
  sum_ns <- sample_n[1] + sample_n[2]
  product_ns <- sample_n[1] * sample_n[2]
  effect_size <- stat_value * sqrt(sum_ns/product_ns);
  return(effect_size)

} # end effect_size_from_stat

#========================================================================================================
# expected_significant_effect
#========================================================================================================
#' @export
expected_significant_effect.TwoSampleZ <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  # Fetch the appropriate critical value and noncentrality parameter as defined for this statistic
  z_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  # We will integrate over the portion of the noncentral sampling distribution greater than z-critical
  lower_bound <- z_critical
  upper_bound <- Inf

  # Required to make a proportional adjustment in the computation of expected effect size
  total_area_above_crit <- 1 - pnorm(z_critical - noncen_parameter)


  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  z_prob_times_value <- function(z)
  {
    dnorm(z - noncen_parameter) * z
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(z_prob_times_value, lower_bound, upper_bound))
  expected_z <- expected_value_numerator/total_area_above_crit

  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_z, sample_n)

  return(expected_effect_size)

} # end expected_significant_effect


#========================================================================================================
# expected_significant_effect_bounded
#========================================================================================================
#' @export
expected_significant_effect_bounded.TwoSampleZ <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{

  # Check that sample_n is a vector of length two
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample z-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")


  z_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  z_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- z_critical_weak
  upper_bound <- z_critical_strong

  total_area_above_strong <- 1 - pnorm(z_critical_strong - noncen_parameter)
  total_area_above_weak <- 1 - pnorm(z_critical_weak - noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong


  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  z_prob_times_value <- function(z)
  {
    dnorm(z - noncen_parameter) * z
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(z_prob_times_value, lower_bound, upper_bound))
  expected_z <- expected_value_numerator/total_area_between

  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_z, sample_n)

  return(expected_effect_size)

} # end expected_significant_effect_bounded

#========================================================================================================
# divide_total_sample
#========================================================================================================
#' @export
divide_total_sample.TwoSampleZ <- function(stat_procedure, n_per_segment)
{
  # split total subject into equal groups (or as close as possible, if total n is odd)
  group_sizes <- numeric(2)
  group_sizes[1] <- floor(n_per_segment/2)
  group_sizes[2]  <- n_per_segment - group_sizes[1]

  return(group_sizes)
} # end divid_total_sample


