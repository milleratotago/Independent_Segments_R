# Class SegmentedResearcher: Holds collection of summary values to describe a specific research scenario



#########################################################################################################
# Class properties
#########################################################################################################
# Passed in to constructor

#   max_n_segments: Maximum number of segments to test.
#   n_per_segment:  SampleN per segment (total N across all samples if there are multiple samples per segment, e.g., 2-sample t-test).
#   alpha_total:    Total probability of Type I error across all segments.
#   alpha_strong:   Stop sampling and reject Ho if p<AlphaStrong in any segment. Defaults to 0 if passed null
#   alpha_weak:     Stop sampling Maximum p value for observed results in each segment. Computed if not provided


#########################################################################################################
# Creation methods - following Wickham
#########################################################################################################
#========================================================================================================
# new_SegmentedResearcher: Default values determined emperically. Users are advised to override.
#========================================================================================================
# Low-level constructor
new_SegmentedResearcher <- function(max_n_segments = 3,
                                    n_per_segment = 50,
                                    alpha_total = 0.05,
                                    alpha_strong = 0,
                                    alpha_weak = NA)
{

  new_seg_res <- list(max_n_segments = max_n_segments,
                      n_per_segment = n_per_segment,
                      alpha_total = alpha_total,
                      alpha_strong = alpha_strong,
                      alpha_weak = alpha_weak)


  attr(new_seg_res, "class") <- "SegmentedResearcher"

  # Return instance
  return(new_seg_res)

} # new_SegmentedResearcher
#========================================================================================================

#========================================================================================================
# validate_SegmentedResearcher
#========================================================================================================
# Checks additional properties of elements, for example, that values are in allowed ranges -- business rules
validate_SegmentedResearcher <- function(segmented_researcher)
{
  # Pull field values for clarity
  max_n_segments <- segmented_researcher$max_n_segments
  alpha_total <- segmented_researcher$alpha_total
  alpha_strong <- segmented_researcher$alpha_strong
  alpha_weak <- segmented_researcher$alpha_weak

  # Make sure that alpha weak is initialised correctly.
  # If a value was supplied (e.g. to improve efficiency in simulations), make sure that it is
  # mathematically correct.
  # If no value was supplied, compute it here.

  # If supplied, check the maths (pers. comm. Miller, 2020)
  tolerance <- 0.0001
  if (!is.na(alpha_weak)) {

    alpha_dif <- alpha_weak - alpha_strong
    pred_alpha_total <- alpha_strong*(1 - alpha_dif^(max_n_segments-1)) / (1- alpha_dif) + alpha_weak * alpha_dif^(max_n_segments-1)

    if (abs(pred_alpha_total - alpha_total) > tolerance) {
      stop('Cannot proceed because provided values of alpha_strong, and alpha_weak do not produce alpha_total')
      }
    } else {

      # Computes appropriate alpha weak value to obtain desired alpha total. See method for computation.
      alpha_weak <- alpha_weak_computation(segmented_researcher, max_n_segments, alpha_total, alpha_strong)

      # Assign to instance. Objects are passed by reference in this construction, so this works.
      segmented_researcher$alpha_weak <- alpha_weak
    } # end alpha_weak passed as default NA
  # Return updated object instance
  return(segmented_researcher)

} # end validate_SegmentedResearcher
#========================================================================================================


#========================================================================================================
# SegmentedResearcher - public-facing ctor
#========================================================================================================
#-----------------
# Roxygen comments

#' SegmentedResearcher Cosntructor
#'
#' A SegmentedResearcher instance holds a collection of values which
#' define a specific research scenario.
#'
#'
#' \code{SegmentedResearcher}  is the public-facing SegmentedResearcher constructor. It
#' provides typical default values for its parameters, but in most cases the user
#' will supply their own values to define the study scenario they wish to explore.
#'
#' @param max_n_segments An integer specifying the maximum number of segments to
#'   be run prior to termination
#' @param n_per_segment A \code{numeric} specifying the number of subjects to be
#'   tested in each segment. Note that this value should comprise all subjects
#'   in a segment, summed across any groups.
#' @param alpha_total Required Type 1 error probability across the entire study.
#'   Defaults to 0.05
#' @param alpha_strong For segmented hypothesis testing decision rule: If
#'   p-observed for any segment is less than alpha_strong, reject H0 and stop
#'   testing.
#' @param alpha_weak For decision rule: If p-observed for any segment is
#'   greater than alpha_weak, fail to reject H0 and stop testing. NB: The correct value
#'   of alpha_weak is determined entirely by the values of alpha_total and
#'   alpha_strong and must be derived numerically from those parameters. In most
#'   cases, a value SHOULD NOT be supplied for this argument. When no value is
#'   supplied, the correct value is automatically computed in succeeding
#'   computations.
#' @return Validated SegmentedResearcher instance
#'
#' @export
SegmentedResearcher <- function(max_n_segments = 3,
                                n_per_segment = 50,
                                alpha_total = 0.05,
                                alpha_strong = 0,
                                alpha_weak = NA)
{

    # Create instance
    instance <- new_SegmentedResearcher(max_n_segments, n_per_segment, alpha_total, alpha_strong, alpha_weak)

    # Validate instance. Initialises alpha weak if necessary.
    valid_instance <- validate_SegmentedResearcher(instance)

    # Return validated instance
    return(valid_instance)

} # end SegmentedResearcher
#========================================================================================================


#########################################################################################################
# Public methods
#########################################################################################################
#========================================================================================================
# print.SegmentedResearcher: System will dispatch to this method on print(instance)
#========================================================================================================
#' @export
print.SegmentedResearcher <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}


#########################################################################################################
# Internal methods
#########################################################################################################
#========================================================================================================
# alpha_weaK_computation - For given values of max_n_segments and alpha_strong, compute the value of alpha_weak
# that will give the desired overall value of alpha_total (total Type 1 error rate). See Miller & Ulrich, 2020
# for detailed discussion.
#========================================================================================================
alpha_weak_computation <- function(segmented_researcher, max_n_segments, alpha_total, alpha_strong, tolerance = 1e-8)
{
  # Logically, alpha_strong can not be larger than alpha_total
  if (alpha_strong > alpha_total)
  {
    stop("Alpha strong must be less than or equal to the overall alpha")
  }

  # If alpha_strong has been specified as zero (boundary case) alpha weak has an analytical solution. Otherwise
  # use numerical search function fzero to minimise function which computes the difference between alpha_total obtained
  # with a target value of alpha weak, and the desired alpha_total.

  if (alpha_strong == 0) {
    final_alpha_weak <- alpha_total^(1/max_n_segments)
  } else {

    # Error function to minimise. Computation by definition. Cf. Miller & Ulrich, 2020.
    compute_error <- function(try_alpha_weak)
    {
      diff <- try_alpha_weak - alpha_strong
      diff_accumulated  <- diff^(max_n_segments - 1)
      err <- alpha_strong * (1 - diff_accumulated )/(1 - diff )+try_alpha_weak * diff_accumulated - alpha_total
      return(err)
    }

    # Define sensible range for numeric search
    aw_range <- c(alpha_strong , alpha_total^(1/max_n_segments ))

    # Perform numeric search to minimise function
    search_result <- pracma::fzero(compute_error, aw_range )

    # Pull value of x that minimises error function
    final_alpha_weak <- search_result$x
  } # end if alpha_strong != 0


  return (final_alpha_weak)
}

