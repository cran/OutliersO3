O3plotColours <- function(colours = c("khaki", "yellow", "red", "lightgreen",
                           "lightblue", "red", "slategray1", "slategray2",
                           "slategray3", "slategray4", "orange", "red"), colors) {
#Allow either colours or colors
    if(!missing(colors)) {
    if(!missing(colours))
      stop("Please specify colours or colors but not both.")
    else
      colours <- colors
  }
  
  list(colours=colours)
}