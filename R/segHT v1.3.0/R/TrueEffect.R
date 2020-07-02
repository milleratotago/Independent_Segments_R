# Class TrueEffect: Holds effect size (TrueEffect$effect_size) and probablity (TrueEffect$effect_size_probability),
# information needed to describe the true state of the world.
#
# An effect size of 0 implies that the null hypothesis (H0) is true.
#
# An experimental scenario may have one or more true effects, each with a probability of occurrence,
# (more generally called a "base rate"). The most typical situations are:
#
#    1) a single effect with probability 1, meaning we consider only those cases where the effect is present
#
#    2) a single effect with probability (base rate) < 1. The researcher expects the effect to be present in this
#    proportion of studies across an infinite number of studies in the area. The researcher expects there to be no effect (i.e. H0 to be true) in the remaining 1- base rate
#    proportion of studies.


#########################################################################################################
# Class properties
#########################################################################################################
#   effect_size:              Default = 0
#   effect_size_probability:  Default = 1
#########################################################################################################


#########################################################################################################
# Creation methods (following Wickham)
#########################################################################################################
#========================================================================================================
# new_TrueEffect - low-level ctor.
#========================================================================================================
new_TrueEffect <- function(effect_size = 0, effect_size_probability = 1)
{
  new_true_effect <- list(effect_size = effect_size,
                          effect_size_probability = effect_size_probability)


  attr(new_true_effect, "class") <- "TrueEffect"

  # Return instance
  return(new_true_effect)

} # new_TrueEffect
#========================================================================================================


#========================================================================================================
# validate_TrueEffect - error checking for raw TrueEffect instance
#========================================================================================================
validate_TrueEffect <- function(true_effect)
{

  if ((true_effect$effect_size < 0) || (true_effect$effect_size > 1)){
    stop("Effect size value must be between 0 and 1, inclusive")
  }


  if ((true_effect$effect_size_probability < 0) || (true_effect$effect_size_probability > 1)){
    stop("Effect size probability values must be between 0 and 1, inclusive")
  }

} # end validate_TrueEffect
#========================================================================================================


#========================================================================================================
# TrueEffect - Public ctor
#========================================================================================================
#-----------------
# Roxygen comments

#' TrueEffect
#'
#' Class TrueEffect holds effect size (TrueEffect$effect_size) and probablity (TrueEffect$effect_size_probability),
#' information needed to describe the true state of the world.
#'
#' An effect size of 0 implies that the null hypothesis (H0) is true.
#
#' An experimental scenario may have one or more true effects, each with a probability of occurrence,
#' (more generally called a "base rate"). The most typical situations are:
#'
#'    1) a single effect with probability 1, meaning we consider only those cases where the effect is present
#'
#'    2) a single effect with probability (base rate) < 1. The researcher expects the effect to be present in this
#'    proportion of studies across an infinite number of studies in the area. The researcher expects there to be no effect (i.e. H0 to be true) in the remaining 1- base rate
#'    proportion of studies.
#'
#' \code{TrueEffect}  is the public-facing TrueEffect constructor. It creates an unvalidated class instance with
#'  new_TrueEffect, and passes it to validate_TrueEffect. An exception is thrown if either parameter value
#'  is illegal (outside range {0..1}). It returns the instance if validation does not fail.
#'
#' @param effect_size A \code{numeric} that describes the effect size on a scale from 0 to 1 (e.g. Cohen's d or eta-squared)
#' @param effect_size_probability A \code{numeric} between 0 and 1; the base rate of this effect size. Set at 1 if only interested in cases where
#'   the effect is known to exist. Set at less than 1 to allow probability of Ho true.
#' @return Validated TrueEffect instance
#' @export
TrueEffect <- function(effect_size = 0, effect_size_probability = 1)
{

  # Create instance
  instance = new_TrueEffect(effect_size, effect_size_probability)

  # Validate instance. Execution halts here if things go wrong
  validate_TrueEffect(instance)

  # Return validated instance
  return(instance)

} # end TrueEffect ctor
#========================================================================================================


#########################################################################################################
# Public methods
#########################################################################################################
#========================================================================================================
# print.TrueEffect: System will dispatch to this method on print(instance)
#========================================================================================================
#' @export
print.TrueEffect <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}
