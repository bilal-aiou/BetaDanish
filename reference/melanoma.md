# Malignant Melanoma Survival After Surgery

Survival times for 205 patients with malignant melanoma after surgery.
This rich clinical dataset includes multiple covariates and heavy
censoring.

## Usage

``` r
melanoma
```

## Format

A data frame with 205 rows and 6 columns:

- time:

  Survival time in days

- status:

  Event indicator (1 = died from melanoma, 0 = alive, 2 = other death)

- thickness:

  Tumor thickness in mm

- sex:

  Patient sex (1 = male, 0 = female)

- age:

  Patient age in years

- ulcer:

  Ulceration indicator (1 = present, 0 = absent)

- year:

  Year of operation

## Source

Andersen, P. K., Borgan, O., Gill, R. D., & Keiding, N. (1993).
Statistical Models Based on Counting Processes. Springer.

## Examples

``` r
data(melanoma)
# \donttest{
# Treat status 1 as event, others as censored
melanoma$event <- ifelse(melanoma$status == 1, 1, 0)
fit <- fit_betadanish(survival::Surv(time, event) ~ age + thickness, data = melanoma)
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=1.82804e-322): z=3.83395e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=2.17572e-319): z=4.69193e-318, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=8.01538e-319): z=1.72854e-317, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=1.10962e-318): z=2.39293e-317, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=7.40529e-317): z=1.59699e-315, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=1.02087e-316): z=2.20156e-315, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0655, b=1.65454e-10, x=1, y=1.65148e-316): z=3.56149e-315, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.035, b=7.41774e-05, x=1, y=4.05134e-322): z=8.7252e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.1909, b=0.000334684, x=1, y=1.03754e-322): z=2.04543e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2131, b=0.000111587, x=1, y=6.91692e-323): z=1.43279e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.236, b=0.000168256, x=1, y=1.08694e-322): z=2.25294e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2135, b=0.000111982, x=1, y=2.66795e-322): z=5.52859e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2127, b=0.000110472, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2138, b=0.000112071, x=1, y=3.85371e-322): z=7.9841e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2128, b=0.00011055, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.214, b=0.000112153, x=1, y=5.03947e-322): z=1.04396e-320, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2129, b=0.000110627, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2143, b=0.000112233, x=1, y=6.52167e-322): z=1.35078e-320, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.213, b=0.000110703, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2127, b=0.0001104, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2444, b=0.000231612, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2333, b=0.000196461, x=1, y=4.79244e-322): z=9.83191e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3151, b=0.000269307, x=1, y=9.38725e-323): z=2.05531e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3476, b=0.000274228, x=1, y=4.44659e-323): z=8.2509e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.4569, b=0.000288208, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2407, b=0.000195453, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2629, b=0.000313182, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3236, b=0.000303598, x=1, y=1.4822e-323): z=4.10074e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2415, b=0.00020378, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3536, b=0.000300654, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2419, b=0.000203894, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3884, b=0.000302965, x=1, y=6.91692e-323): z=1.44267e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3884, b=0.000302965, x=1, y=8.89318e-323): z=1.85769e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3884, b=0.000302965, x=1, y=8.89318e-323): z=1.85769e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2423, b=0.000204034, x=1, y=9.88131e-324): z=2.02567e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.4222, b=0.000304559, x=1, y=1.23516e-322): z=2.48021e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2428, b=0.000204171, x=1, y=1.4822e-323): z=4.10074e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.461, b=0.000314021, x=1, y=4.94066e-323): z=1.03754e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.461, b=0.000314021, x=1, y=7.90505e-323): z=1.65512e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2433, b=0.000204352, x=1, y=3.45846e-323): z=8.20149e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.242, b=0.000203786, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2435, b=0.000204262, x=1, y=2.47033e-323): z=4.10074e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2421, b=0.000203791, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2462, b=0.000199559, x=1, y=2.47033e-323): z=4.10074e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2419, b=0.000203664, x=1, y=9.88131e-324): z=2.02567e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2418, b=0.00020369, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2419, b=0.00020401, x=1, y=9.88131e-324): z=2.02567e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2418, b=0.000203758, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2954, b=0.000304586, x=1, y=3.45846e-323): z=8.20149e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2422, b=0.000204365, x=1, y=3.95253e-323): z=8.20149e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2419, b=0.000203839, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2418, b=0.000203734, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2418, b=0.000203713, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3427, b=0.000213867, x=1, y=4.94066e-323): z=1.02766e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2457, b=0.000204106, x=1, y=3.80431e-322): z=7.79142e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2426, b=0.000203789, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.242, b=0.000203725, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2419, b=0.000203712, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.3472, b=0.000235721, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2426, b=0.000203948, x=1, y=9.88131e-324): z=2.02567e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.242, b=0.000203757, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2419, b=0.000203719, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.2418, b=0.000203712, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9315, b=0.000184449, x=1, y=1.97626e-323): z=4.24896e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.8023, b=0.000407754, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.8023, b=0.000407754, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9323, b=0.000184264, x=1, y=9.88131e-324): z=2.12448e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.0654, b=0.000208488, x=1, y=4.94066e-323): z=1.06718e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9324, b=0.000184242, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9409, b=0.000185646, x=1, y=1.33398e-322): z=2.96439e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9351, b=0.000184784, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.4917, b=0.000270942, x=1, y=1.97626e-323): z=4.34778e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.4917, b=0.000270942, x=1, y=1.97626e-323): z=4.34778e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9379, b=0.00018518, x=1, y=1.97626e-323): z=4.24896e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9348, b=0.000184725, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.5001, b=0.000256782, x=1, y=1.72923e-322): z=3.913e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9381, b=0.000185123, x=1, y=1.4822e-323): z=4.24896e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9349, b=0.000184732, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=22.6335, b=0.000280164, x=1, y=1.63042e-322): z=3.49798e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9391, b=0.000185271, x=1, y=2.96439e-323): z=6.37345e-322, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9353, b=0.000184777, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9427, b=0.000185737, x=1, y=1.82804e-322): z=3.81419e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9361, b=0.00018489, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9348, b=0.000184721, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9346, b=0.000184687, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9356, b=0.000184476, x=1, y=8.89318e-323): z=1.90709e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9347, b=0.00018464, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9345, b=0.000184672, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9345, b=0.000184679, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.9345, b=0.00018468, x=1, y=4.94066e-324): z=0, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.0971, b=0.000245906, x=1, y=9.88131e-323): z=2.03555e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.0971, b=0.000245906, x=1, y=9.88131e-323): z=2.03555e-321, b*z == 0 underflow, hence inaccurate pbeta()
# }
```
