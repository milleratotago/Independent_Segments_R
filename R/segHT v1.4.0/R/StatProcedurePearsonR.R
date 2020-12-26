# Class PearsonRStatType descends from StatProcedureBase

# See class OneSampleT (StatProcedureOneSampleT.R) for detailed explanation of method logic.

# Small numerical differences (beyond 3 significant digits) between matlab and R implementations
# because of slight variations in the computation of the r distribution.

# All the computations ported from original source code (Miller, 2020)

#' @importFrom SuppDists pPearson
#' @importFrom SuppDists qPearson
#' @importFrom SuppDists dPearson
#' @importFrom pracma integral
#' @importFrom utils capture.output

########################################
# Creation methods - following Wickham
########################################

#====================================================================
new_PearsonR <- function(min_sample_n, display_abbrev, display_name)
{
  new_pearson_r <- list(min_sample_n = min_sample_n,
                        display_abbrev = display_abbrev,
                        display_name = display_name)


  attr(new_pearson_r, "class") <- c("PearsonR", "StatProcedureBase")

  return(new_pearson_r)
}

#====================================================================
validate_PearsonR <- function(pearson_r)
{
  # Do any validation here. Use stop(msg) to exit if things go wrong
}

#====================================================================
#-----------------
# Roxygen comments

#' PearsonR Constructor
#'
#' PearsonR is a child of StatProcedureBase, representing a Pearson Product Moment correlation
#'
#' \code{PearsonR}  is the public-facing PearsonR constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated PearsonR instance
#'
#' @export
PearsonR <- function(min_sample_n = 5, display_abbrev = "r", display_name = 'Pearson r')
{
  instance <- new_PearsonR(min_sample_n, display_abbrev, display_name)
  validate_PearsonR(instance)
  return(instance)
}

#################################
# Instance methods
#################################
#========================================================================================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
#========================================================================================================
#' @export
compute_power.PearsonR <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  r_crit_cumul_prob <- 1 - alpha_one_tailed

  if (is.null(critical_value))
  {
    r_crit <- SuppDists::qPearson(r_crit_cumul_prob, sample_n)
  }
  else
  {
    r_crit <- critical_value
  }

  # For the Pearson, "effect size" is simply the true value of Rho

  # The SuppDists::pPearson function accepts a value of r and returns the cumulative probability (i.e. the chances
  # that a value less than this occurs).
  # Power (pr of rejecting) is 1 - that value.
  power <- 1 - SuppDists::pPearson(r_crit, sample_n, rho = effect_size)

  return(power)

} # end compute_power.PearsonR


#========================================================================================================
# p_level(pearson_r, observed_stat, sample_n)
#========================================================================================================
#' @export
p_level.PearsonR <- function(stat_procedure, observed_stat, sample_n)
{

  # the observed sample exceeds this proportion of the sampling distribution of r when Ho is true
  cumul_prob_obs <- SuppDists::pPearson(observed_stat, sample_n)

  # by definition of p-value
  p <- 1 - cumul_prob_obs

  return(p)

} # end p_level

#========================================================================================================
# critical_value(pearson_r, alpha_one_tailed, sample_n)
#========================================================================================================
#' @export
critical_value.PearsonR <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  #print(stat_procedure)

  r_crit_cumul_prob <- 1 - alpha_one_tailed

  r_critical <- SuppDists::qPearson(r_crit_cumul_prob, sample_n)

  return(r_critical)

} # end critical_value

#========================================================================================================
# noncentrality(sample_n, effect_size): Parameter to describe how far the r-distribution is shifted
# from the standard (i.e. that which exists when effect size = 0, Ho is true). For pearson r, this is simply
# the true correlation rho, so it doesn't really need its own method, but this is included for interface consistency
# across the class family.
#========================================================================================================
#' @export
noncentrality.PearsonR <- function(stat_procedure, sample_n, effect_size)
{
  # By definition
  noncen <- effect_size
  return(noncen)

} # end noncentrality

#========================================================================================================
# effect_size_from_stat(stat_value, sample_n):  For pearson r, this is simply the true correlation rho,
# so it doesn't really need its own method, but this is included for interface consistency across the class family.
#========================================================================================================
#' @export
effect_size_from_stat.PearsonR <- function(stat_procedure, stat_value, sample_n)
{
  return(stat_value)

} # end effect_size_from_stat

#========================================================================================================
# expected_significant_effect: Equivalent to the similar methods for t-tests, except that the r distribution
# methods SuppDists::pPearson, SuppDists::qPearson, etc. accept the full sample n, not the precomputed degrees of freedom.
#========================================================================================================
#' @export
expected_significant_effect.PearsonR <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  r_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- r_critical
  upper_bound <- 1 # Maximum possible value of pearson r

  total_area_above_crit <- 1 - SuppDists::pPearson(r_critical, sample_n, noncen_parameter)


  # The function to integrate over. The variables N and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below.

  N <- sample_n
  r_prob_times_value <- function(r)
  {
    SuppDists::dPearson(r, N = N, rho = noncen_parameter) * r
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(r_prob_times_value, lower_bound, upper_bound))
  expected_r <- expected_value_numerator/total_area_above_crit

  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_r, sample_n)

  return(expected_effect_size)

} # end expected_significant_effect

#========================================================================================================
#expected_significant_effect_bounded
#========================================================================================================
#' @export
expected_significant_effect_bounded.PearsonR <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{

  r_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  r_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)

  lower_bound <- r_critical_weak
  upper_bound <- r_critical_strong

  total_area_above_strong <- 1 - SuppDists::pPearson(r_critical_strong, sample_n, noncen_parameter)
  total_area_above_weak <- 1 - SuppDists::pPearson(r_critical_weak, sample_n, noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong

  # The function to integrate over. The variables N and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below.

  N <- sample_n
  r_prob_times_value <- function(r)
  {
    SuppDists::dPearson(r, N = N, rho = noncen_parameter) * r
  }

  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(r_prob_times_value, lower_bound, upper_bound))
  expected_r <- expected_value_numerator/total_area_between

  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_r, sample_n)

  return(expected_effect_size)


} # end expected_significant_effect_bounded

#========================================================================================================
# divide_total_sample - Not really required for r, but included for consistency
#========================================================================================================
#' @export
divide_total_sample.PearsonR <- function(stat_procedure, n_per_segment)
{
  # For simple r, all subjects are in a single group
  return(n_per_segment)
} # end divid_total_sample
