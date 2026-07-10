# Constrain survival probabilities according to hazards in a lifetable Recalculated constrained survival probabilities (by week) as the lower of the original unadjusted survival probability and the survival implied by the given lifetable (assumed indexed as years).

Constrain survival probabilities according to hazards in a lifetable
Recalculated constrained survival probabilities (by week) as the lower
of the original unadjusted survival probability and the survival implied
by the given lifetable (assumed indexed as years).

## Usage

``` r
constrain_survprob(
  survprob1,
  survprob2 = NA,
  lifetable = NA,
  timevec = 0:(length(survprob1) - 1)
)
```

## Arguments

- survprob1:

  (Unconstrained) survival probability value or vector

- survprob2:

  Optional survival probability value or vector to constrain on (default
  = NA)

- lifetable:

  Lifetable (default = NA)

- timevec:

  Vector of times corresponding with survival probabilities above

## Value

Vector of constrained survival probabilities

## Examples

``` r
ltable <- tibble::tibble(lttime=0:20, lx=c(1,0.08,0.05,0.03,0.01,rep(0,16)))
survprob <- c(1,0.5,0.4,0.2,0)
constrain_survprob(survprob, lifetable=ltable)
#> [1] 1.0 0.5 0.4 0.2 0.0
timevec <- 100*(0:4)
constrain_survprob(survprob, lifetable=ltable, timevec=timevec)
#> [1] 1.00000000 0.05200139 0.01201389 0.00000000 0.00000000
survprob2 <- c(1,0.45,0.35,0.15,0)
constrain_survprob(survprob, survprob2)
#> [1] 1.00 0.45 0.35 0.15 0.00
```
