##' Class "AstralPlane"
##' This class stores astral phylogenetic data
##'
##' @name AstralPlane-class
##' @docType class
##' @slot fileName path to the ASTRAL output tree file
##' @slot samples character vector of sample (tip) names in the tree
##' @slot phylo rooted phylo object from ape
##' @slot nodeData data.frame of per-node ASTRAL statistics (q1-q3, f1-f3, pp1-pp3, QC, EN)
##' @slot edgeData data.frame mapping edges to node pairs with branch lengths
##' @slot concordanceFactorData data.frame of IQ-TREE concordance factor data (empty if unused)
##' @exportClass AstralPlane
##' @keywords classes

#Sets the class for the phylogenetic object
setOldClass(Classes = "phylo")

setClass("AstralPlane", slots=list(fileName = "character",
                                   samples="character",
                                   phylo = "phylo",
                                   nodeData="data.frame",
                                   edgeData="data.frame",
                                   concordanceFactorData = "data.frame"))

