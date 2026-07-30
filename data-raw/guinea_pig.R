## Guinea pig survival data (Bjerkedal 1960).
##
## Survival times in days of 72 guinea pigs injected with virulent tubercle
## bacilli, from the principal regimen of the study. This is a complete sample:
## every animal was observed to death, so status is 1 throughout.
##
## Included because it is the one dataset in the thesis on which the
## four-parameter Beta-Danish model attains a genuine finite interior optimum,
## with b-hat = 3.64 (SE 1.20), placing b about 2.2 standard errors clear of the
## b = 1 identifiability ridge. Every other application sits on the flat
## (a, c) direction with a-hat < 1.
##
## Re-run with:  source("data-raw/guinea_pig.R")

gp <- c( 12,  15,  22,  24,  24,  32,  32,  33,  34,  38,
         38,  43,  44,  48,  52,  53,  54,  54,  55,  56,
         57,  58,  58,  59,  60,  60,  60,  60,  61,  62,
         63,  65,  65,  67,  68,  70,  70,  72,  73,  75,
         76,  76,  81,  83,  84,  85,  87,  91,  95,  96,
         98,  99, 109, 110, 121, 127, 129, 131, 143, 146,
        146, 175, 175, 211, 233, 258, 258, 263, 297, 341,
        341, 376)

stopifnot(length(gp) == 72L, !anyNA(gp), all(gp > 0),
          identical(gp, sort(gp)))

guinea_pig <- data.frame(time = as.numeric(gp), status = 1L)

save(guinea_pig, file = "data/guinea_pig.rda", version = 2, compress = "xz")
