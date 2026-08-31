#' @title findIQTREE
#'
#' @description Locates an IQ-TREE executable and reports its major version.
#'   IQ-TREE 2 installs its binary as "iqtree2" and IQ-TREE 3 installs it as
#'   "iqtree", so a hard-coded name works on one machine and fails on the next.
#'   Both names are tried, "iqtree2" first, and the version is read back from
#'   the executable rather than inferred from its name.
#'
#' @param iqtree.path path to an IQ-TREE executable, or to the directory
#'   containing it; use NULL to search the system PATH
#'
#' @param quiet if TRUE suppresses the message naming the executable found
#'
#' @return a list with "path" (the executable), "version" (the version string
#'   reported by IQ-TREE, or NA if it could not be read) and "major" (the major
#'   version as an integer, or NA)
#'
#' @examples
#'
#' iqtree = findIQTREE()
#' iqtree$path
#' iqtree$major
#'
#' @export

findIQTREE = function(iqtree.path = NULL,
                      quiet = FALSE) {

  #IQ-TREE 2 ships as iqtree2, IQ-TREE 3 ships as iqtree
  exe.names = c("iqtree2", "iqtree")

  if (is.null(iqtree.path)){
    found = unname(Sys.which(exe.names))
    found = found[nzchar(found)]
    if (length(found) == 0){ stop("No iqtree2 or iqtree executable was found on the PATH. Provide iqtree.path.") }
    iqtree.exe = found[1]
  } else if (dir.exists(iqtree.path)){
    candidates = file.path(iqtree.path, exe.names)
    candidates = candidates[file.exists(candidates)]
    if (length(candidates) == 0){ stop("No iqtree2 or iqtree executable was found in ", iqtree.path) }
    iqtree.exe = candidates[1]
  } else {
    iqtree.exe = iqtree.path
    if (file.exists(iqtree.exe) == FALSE){ stop("The IQ-TREE executable could not be found at ", iqtree.exe) }
  }

  #Asks the executable what it is, rather than trusting its file name. A version
  #that cannot be read is reported as NA rather than being fatal: a wrapper
  #script or an unusual build should still be allowed to run.
  version.text = suppressWarnings(tryCatch(system2(iqtree.exe, "--version", stdout = TRUE, stderr = TRUE),
                                           error = function(e) character(0)))
  version.line = grep("IQ-TREE", version.text, value = TRUE)
  version.number = NA_character_
  major.version = NA_integer_

  if (length(version.line) > 0){
    parsed = regmatches(version.line[1], regexpr("[0-9]+\\.[0-9]+(\\.[0-9]+)?", version.line[1]))
    if (length(parsed) > 0){
      version.number = parsed
      major.version = as.integer(sub("\\..*$", "", parsed))
    }
  }

  if (quiet == FALSE){
    if (is.na(version.number)){
      cat(paste0("Using IQ-TREE at ", iqtree.exe, " (version could not be read)\n"))
    } else {
      cat(paste0("Using IQ-TREE ", version.number, " at ", iqtree.exe, "\n"))
    }
  }

  return(list(path = iqtree.exe, version = version.number, major = major.version))

}#end function
