library(SoilTaxonomy)
library(purrr)

data('ST_unique_list', package = 'SoilTaxonomy')


# formative element parsing requires taxa one level higher
x <- OrderFormativeElements(ST_unique_list$suborder)
any(is.na(x$defs$order))
any(is.na(x$char.index))

# formative element parsing requires taxa one level higher
x <- OrderFormativeElements(ST_unique_list$greatgroup)
any(is.na(x$defs$order))

# formative element parsing requires taxa one level higher
x <- OrderFormativeElements(ST_unique_list$subgroup)
any(is.na(x$defs$order))


##

# formative element parsing requires taxa one level higher
x <- SubOrderFormativeElements(ST_unique_list$greatgroup)
# OK
any(is.na(x$defs$order))
# some bugs
any(is.na(x$char.index))

# "fibristels" "folistels"  "hemistels"  "glacistels" "sapristels"
idx <- which(is.na(x$char.index))
ST_unique_list$greatgroup[idx]

# fixed
x <- 'calciargids'
SubOrderFormativeElements(x)

# correct parsing of formative elements, but not character position
x <- 'folistels'
SubOrderFormativeElements(x)

# missing subgroup formative elements: petrosalic
x <- 'petrosalic anhyturbels'
OrderFormativeElements(x)
SubOrderFormativeElements(x)
GreatGroupFormativeElements(x)
SubGroupFormativeElements(x)


# hmm
x <- c('argixerolls', 'folistels', 'calciargids', 'fartyblogfish')
SubOrderFormativeElements(x)



# works !
x <- SubOrderFormativeElements(ST_unique_list$subgroup)
any(is.na(x$defs$order))
any(is.na(x$char.index))

# find errors
s <- ST_unique_list$subgroup
names(s) <- s
s.test <- map(s, safely(SubOrderFormativeElements))
s.test <- transpose(s.test)

# find errors
idx <- ! sapply(s.test$error, is.null)
names(s.test$error[idx])


# how about some problematic subgroups
x <- c('argixerolls', 'acrustoxic kanhaplustults')
SubOrderFormativeElements(x)




##
## look for missing subgroup formative elements
##

data('ST_formative_elements', package = 'SoilTaxonomy')

# all subgroups
s <- ST_unique_list$subgroup

# tokenize via whitespace
tok <- SoilTaxonomy:::.tokenizeST(s)

# extract all but last token
tok.sg <- sapply(tok, function(i) {
  len <- length(i)
  idx <- seq(from=1, to=len - 1)
  return(i[idx])
})

# unique vector of subgroup formative elements
tok.sg <- unique(unlist(tok.sg))

# find those missing in the dictionary
missing.in.dictionary <- sort(setdiff(tok.sg, ST_formative_elements$subgroup$element))

# find those missing in actual use
missing.in.ST <- sort(setdiff(ST_formative_elements$subgroup$element, tok.sg))

# save for later... need to define these
cat(missing.in.dictionary, sep='\n')


