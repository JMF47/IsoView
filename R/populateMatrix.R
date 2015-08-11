#' Helper function for main visualize function
#'
#'@param query Query
#'@param subject Subject
#'@export
#'@examples
#'populateMatrix()
#'
populateMatrix = function(query, subject){
      info = findOverlaps(query, subject)
      ind = 1:length(subject)
      return((ind %in% subjectHits(info)))
}