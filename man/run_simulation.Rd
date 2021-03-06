% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{run_simulation}
\alias{run_simulation}
\title{run_simulation}
\usage{
run_simulation(nreps, niter, mc.cores, simtype = "permutation",
  J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = 0.2,
  optimizer = "bobyqa", re.form = NULL, permtype = "within",
  boottype = "parametric")
}
\arguments{
\item{nreps}{Number of repetitions using new data.}

\item{niter}{For each simulated data set, number of permutation or bootstrap
iterations.}

\item{mc.cores}{Number of cores to use.}

\item{simtype}{Permutations or bootstrap.}

\item{permtype}{This can be used to specify either the permutation type
('within', 'between', or 'between_within'), or if simtype is 'bootstrap',
and boottype is 'cases', it is used to specify how cases are resampled,
('groups', 'within', 'both').}

\item{boottype}{}
}
\description{
run_simulation
}
