# Class StatProcedureBase

# Base class for objects representing the statistical test being used in a research
# scenario. Children are 1 and 2 group t-tests, correlation and 1 and 2 group z-test
# following Miller & Ulrich 2020


#########################################################################################################
# Class properties: Values by definition for each child test, set in ctor
#########################################################################################################

# min_sample_n
# display_abbrev
# display_name


#########################################################################################################
# Creation methods (following Wickman)
#########################################################################################################

#========================================================================================================
# new_StatProcedureBase: Children will, in typical cases, provide hard-coded values for all fields.
# Values for display_abbrev must match SegmentedHypTestEngine$stat_procedure_factory switch statement.
#========================================================================================================
# Low-level constructor - following Wickman
new_StatProcedureBase <- function(min_sample_n, display_abbrev, display_name)
{

  new_stat_procedure_base <- list(min_sample_n = min_sample_n,
                                  display_abbrev = display_abbrev,
                                  display_name = display_name)


  attr(new_stat_procedure_base, "class") <- "StatProcedureBase"

  # Return instance
  return(new_stat_procedure_base)

} # new_StatProcedureBase
#========================================================================================================


#========================================================================================================
# validate_StatProcedureBase
#========================================================================================================
# Checks additional properties of elements, for example, that values are in allowed ranges -- business rules
validate_StatProcedureBase <- function(stat_procedure_base)
{
  if (stat_procedure_base$min_sample_n < 1) {
    stop("Minimum sample size must be greater than 0")
  }

} # end validate_StatProcedureBase
#========================================================================================================


#========================================================================================================
# StatProcedureBase - Public ctor
#========================================================================================================
#-----------------
# Roxygen comments

#' StatProcedureBase Placeholder Ctor
#'
#' StatProcedureBase is the abstract base class for all stat procedure classes.
#'
#' \code{StatProcedureBase}  should not be called directly. Descend all child classes
#' from StatProcedureBase by assigning the class attribute in the child constructor.
#' For example, the low-level constructor for class OneSampleT contains the statement:
#'
#' \code{attr(new_one_sample_t, "class") <- c("OneSampleT", "StatProcedureBase")}
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#'
#' @export
StatProcedureBase <- function(min_sample_n, display_abbrev, display_name)
{
  # Create instance
  instance = new_StatProcedureBase(min_sample_n, display_abbrev, display_name)

  # Validate instance. Execution halts here if things go wrong
    validate_StatProcedureBase(instance)

  # Return validated instance
  return(instance)

} # end StatProcedureBase ctor
#========================================================================================================


#########################################################################################################
# Class methods
#########################################################################################################

#========================================================================================================
# print.StatProcedureBase: System will dispatch to this method on print(instance)
#========================================================================================================
#' @export
print.StatProcedureBase <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}


# Remaining family methods are all abstract in the base type
#========================================================================================================
# compute_power: Given a directional alpha, anticipated effect size and sample n, return power for the child test.
# These have analytical solutions for each test. See children for details
#========================================================================================================
compute_power.StatProcedureBase <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  # Abstract
} # end compute_power


#========================================================================================================
# p_level: Given an observed stat value and sample n, compute the p-level for the child test.
# These are computed by selecting from the appropriate sampling distribution
#========================================================================================================
p_level.StatProcedureBase <- function(stat_procedure, observed_stat, sample_n)
{
  # Abstract
} # end p_level

#========================================================================================================
# critical_value: Given alpha and sample n, compute the critical value for the child test.
# Determined by accessing the cumulative density function for the appropriate sampling distribution
#========================================================================================================
critical_value.StatProcedureBase <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  # Abstract
} # end critical_value

#========================================================================================================
# noncentrality: Provides a measure of the extent to which the null hypothesis is false, determining the
# true distribution from which the observed statistic is drawn. These have analytical solutions for each
# test. See the child code.
#========================================================================================================
noncentrality.StatProcedureBase <- function(stat_procedure, sample_n, effect_size)
{
  # Abstract
} # end noncentrality

#========================================================================================================
# effect_size_from_stat: Measure of effect size derived from observed statistic value and n. These
# have analytical solutiosn. See the child code.
#========================================================================================================
effect_size_from_stat.StatProcedureBase <- function(stat_procedure, stat_value, sample_n)
{
  # Abstract
} # end effect_size_from_stat

#========================================================================================================
# expected_significant_effect: Determines the expected (mean) effect size, if one considers only those
# cases where the null hypothesis is reject (i.e. the observed stat is greater than the critical). For each
# test, this is computed by integrating across stat * stat_density for all values above the critical
# to obtain an expected (mean) value for that portion of the distribution, than dividing by the total area
# under the curve to obtain the correct proportional value. See child code.
#========================================================================================================
expected_significant_effect.StatProcedureBase <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  # Abstract
} # end expected_significant_effect


#========================================================================================================
# expected_significant_effect_bounded: As for expected_significant_effect, above, but integrated only
# over a region determined by two bounding alpha levels.
#========================================================================================================
expected_significant_effect_bounded.StatProcedureBase <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  # Abstract
} # end expected_significant_effect


#========================================================================================================
# divide_total_sample: Convert from n_per_segment (which is always the total subjects used per segment)
# to n_per_group, depending on the number of groups the child test applies to.
#========================================================================================================
divide_total_sample.StatProcedureBase <- function(stat_procedure, n_per_segment)
{
  # Abstract
} # end divide_total_sample



#######################################################################################################
# Generics: Following the S3 method dispatch model. This causes each child instance to invoke its own
# polymorphic method implementation. See Wickman and others for details of the function binding process.
#######################################################################################################

#========================================================================================================
# Generics -- inherited by child classes
#========================================================================================================
#-----------------
# Roxygen comments

#' compute_power
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{compute_power} computes the statistical power of a statistical test
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#' @param critical_value Optional to speed computation. If omitted, this value is computed in the method
#'
#' @return Computed statistical power (probability of rejecting H0 when the effect is present)
#'
#' @export
compute_power <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  UseMethod("compute_power", stat_procedure)
}

#-----------------
# Roxygen comments

#' p_level
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{p_level} computes the p_level for an observed statistical values and sample size
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param observed_stat Observed test value of interest
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Computed p-level
#'
#' @export
p_level <- function(stat_procedure, observed_stat, sample_n)
{
  UseMethod("p_level", stat_procedure)
}

#-----------------
# Roxygen comments

#' critical_value
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{critical_value} computes the critical value of a test for given alpha and n
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Critical value of test for alpha and n
#'
#' @export
critical_value <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  UseMethod("critical_value", stat_procedure)
}

#-----------------
# Roxygen comments

#' noncentrality
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{noncentrality} computes the noncentrality parameter for a given n and effect size
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Computed noncentrality parameter
#'
#' @export
noncentrality <- function(stat_procedure, sample_n, effect_size)
{
  UseMethod("noncentrality", stat_procedure)
}

#-----------------
# Roxygen comments

#' effect_size_from_stat
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#' children
#'
#' \code{effect_size_from_stat} computes the effect size for a test given alpha,
#' n, and an observed value of the test statistic
#'
#' @param stat_procedure Instance for method dispatch. Member of
#'   StatProcedureBase family
#' @param stat_value Value of test statistic
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Computed effect size
#'
#' @export
effect_size_from_stat <- function(stat_procedure,stat_value, sample_n)
{
  UseMethod("effect_size_from_stat", stat_procedure)
}

##-----------------
# Roxygen comments

#' expected_significant_effect
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#'children
#'
#' \code{expected_significant_effect} is used for simulations to demonstrate
#' bias in published observed effect sizes. Determines the expected (mean) effect
#' size, if one considers only those cases where the null hypothesis is rejected
#' (i.e. the observed stat is greater than the critical). For each test, this is
#' computed by integrating across stat * stat_density for all values above the
#' critical to obtain an expected (mean) value for that portion of the
#' distribution, than dividing by the total area under the curve to obtain the
#' correct proportional value. See child code.
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Expected (mean) effect size
#'
#' @export
expected_significant_effect <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  UseMethod("expected_significant_effect", stat_procedure)
}

##-----------------
## Roxygen comments

#' expected_significant_effect_bounded
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#'children
#'
#' \code{expected_significant_effect_bounded} as for expected_significant_effect but for
#' a specified portion of the sampling distribution
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_strong Upper bound in sampling distribution
#' @param alpha_weak Lower bound in sampling distribution
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Expected (mean) effect size
#'
#' @export
expected_significant_effect_bounded <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  UseMethod("expected_significant_effect_bounded", stat_procedure)
}



## Roxygen comments

#' divide_total_sample
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#'children
#'
#' \code{divide_total_sample} returns n per group given total n
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param n_per_segment Total n per segment
#'
#' @return N per group for the associated statistical procedure object
#'
#' @export
divide_total_sample <- function(stat_procedure, n_per_segment)
{
  UseMethod("divide_total_sample", stat_procedure)
}

