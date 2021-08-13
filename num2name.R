num2name <- function(groups) {
  names = list()
  for(i in 1:length(groups)) {
    names[[i]] = aal2.120$name[groups[[i]]]
  }
  return(names)
}