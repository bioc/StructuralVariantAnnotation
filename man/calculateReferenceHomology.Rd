% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BreakpointGRanges.R
\name{calculateReferenceHomology}
\alias{calculateReferenceHomology}
\title{Calculates the length of inexact homology between the breakpoint sequence
and the reference}
\usage{
calculateReferenceHomology(
  gr,
  ref,
  anchorLength = 300,
  margin = 5,
  match = 2,
  mismatch = -6,
  gapOpening = 5,
  gapExtension = 3
)
}
\arguments{
\item{gr}{reakpoint GRanges}

\item{ref}{reference BSgenome}

\item{anchorLength}{Number of bases to consider for homology}

\item{margin}{Number of additional reference bases include. This allows
for inexact homology to be detected even in the presence of indels.}

\item{match}{see pwalign::pairwiseAlignment}

\item{mismatch}{see pwalign::pairwiseAlignment}

\item{gapOpening}{see pwalign::pairwiseAlignment}

\item{gapExtension}{see pwalign::pairwiseAlignment}
}
\value{
A dataframe containing the length of inexact homology between the 
breakpoint sequence and the reference.
}
\description{
Calculates the length of inexact homology between the breakpoint sequence
and the reference
}
