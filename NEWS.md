news:
 one_page: false
 cran_dates: false
 releases:
 - text: "usethis 1.3.0"
   href: https://www.tidyverse.org/articles/2018/02/usethis-1-3-0/
 - text: "usethis 1.0.0 (and 1.1.0)"
   href: https://www.tidyverse.org/articles/2017/11/usethis-1.0.0/

# 15 Jan 2024 - Version 0.2

This version provides additional functionality to the calculation of restricted mean durations in `calc_allrmds`. These estimates may now be constrained by a lifetable (see `calc_ltsurv`) and discounting may now be applied. A vignette describing how to use this functionality is provided ("background-mortality").

There are also some changes to _pkgdown.yml to avoid the theming changes imposed in the 0.6.0 release of bslib 0.6.0, with thanks to @nanxstats.

# 5 Jan 2024 - Version 0.1.1

## Added PSM analysis functions (experimental)

I’ve merged some experimental functions into the main branch following the version 0.1 release package. These functions provide analyses of the constraints on mortality hazards and therefore survival implied by a PSM. They are: `calc_haz_psm()`, `calc_surv_psmpps()`, `pickout_psmhaz`, `graph_psm_hazards()`, and `graph_psm_survs()`.

# 1 Jan 2024 - Version 0.1

This is the initial release of the package, rather belatedly. The code dates to October 2023.

