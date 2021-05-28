# I worked up a quick set of validations to check the consistency of taxonomy reported in SSURGO components.
# 
# Part 0: Identifying issues in the SC database
# -	1203 inactive series
# -	830 series do not have modern taxonomic subgroup
# -	788 of 830 are inactive series: 42 series are active without a modern code (most are established; attached as CSV)
# 
# Part 1: Validating Component:Series:Taxonomic Subgroup correlations in SSURGO against the SC database
# 
# This validation seeks to identify components that have not been updated as their corresponding series have been updated. “Fixing” these things in practice may not be as simple as replacing old with new, but this identifies the scope of the mismatches for further analysis.
# 
# -	I read the “component” table out of the gSSURGO .gdb and filter to remove components that do not at least have subgroup level taxonomy. Then I calculated the 4 or 5 letter subgroup-level taxonomy letter code for taxsubgrp in SSURGO. This analysis can be repeated (without cross-check to SC) to check higher taxa.
# 
# -	A total of 60900 components do not get a subgroup  code assigned using the 12th edition Keys to Soil Taxonomy lookup table. Of those 60900, 54802 are a “series”-level component based on compkind; we attempt to create an update table based on more “recent” entries in SC database.
# 
# -	54802 series-level components are missing a modern subgroup level taxonomy; in this set there are 4947 unique component names. 
# 
# -	4668 of 4947 unique component names (with compkind series and obsolete or missing taxonomy) appear to [at least partially] match a series defined in the SC database.
# 
# -	Create a lookup table that helps relate the component names that are missing modern taxonomy to their [partial] matched series names; then join the SC database table to the component table using calculated series index.
# 
# -	Calculate the subgroup-level taxonomy letter code for the subgroup derived from the more recent entry in the SC database. In most cases, this results in a more modern but “equivalent” taxonomic class relative to what is populated in SSURGO.
# 
# -	53358 components could be updated using table of current series taxonomy derived from partial matching on compname against the SC database (attached as CSV)
# 
# -	1238 components appear to have compkind=”Series” but compname not a named series defined in the SC database(attached as CSV; compname such as "Soils calcareous throughout", "Slopes more than 15 percent" , "Very deep silt loam soils"); these can be “fixed” by updating their compname, compkind, phase etc. but probably require at least minimal local investigation. Some of these may actually be the "same" as dominant component in MU, but differ in slope etc, fixing those semi-automatically could be possible if mapunits are flagged as such                                     

library(sf)
library(soilDB)
library(SoilTaxonomy)

# get the SC database from local NASIS
scdb <- soilDB::get_soilseries_from_NASIS()
scdb$newcode <- SoilTaxonomy::taxon_to_taxon_code(scdb$taxsubgrp)

# missing order
scdb$scdb_taxorder <- SoilTaxonomy::taxon_to_taxon_code(scdb$taxorder)
sum(is.na(scdb$scdb_taxorder) & scdb$soilseriesstatus != "inactive")

# missing suborder
scdb$scdb_taxsuborder <- SoilTaxonomy::taxon_to_taxon_code(scdb$taxsuborder)
sum(is.na(scdb$scdb_taxsuborder) & scdb$soilseriesstatus != "inactive")

# missing great group
scdb$scdb_taxgrtgroup <- SoilTaxonomy::taxon_to_taxon_code(scdb$taxgrtgroup)
sum(is.na(scdb$scdb_taxgrtgroup) & scdb$soilseriesstatus != "inactive")

# missing subgroup
scdb$scdb_taxsubgroup <- SoilTaxonomy::taxon_to_taxon_code(scdb$taxsubgrp)
sum(is.na(scdb$scdb_taxsubgroup) & scdb$soilseriesstatus != "inactive")

decompcodes <- data.table::rbindlist(lapply(decompose_taxon_code(scdb$newcode), 
                                            function(x) data.table::data.table(t(x))), fill = T)
decompcodes[] <- lapply(decompcodes, as.character)
colnames(decompcodes) <- c("calc_taxorder","calc_taxsuborder","calc_taxgrtgroup","calc_taxsubgrp")

scdb <- cbind(scdb, decompcodes)

# 1203 inactive series in SC db
sum(scdb$soilseriesstatus == "inactive")

# 830 series without modern taxonomic subgroup codes in SC db
sum(is.na(scdb$newcode))

# 788 of 830 that are not "modern" are inactive series in SC db
sum(is.na(scdb$newcode) & scdb$soilseriesstatus == "inactive")

# therefore, 49 series are active without modern order to subgroup code (and 42 are establishes)
fix_series <- scdb[is.na(scdb$newcode) & scdb$soilseriesstatus != "inactive" |
                     is.na(scdb$scdb_taxgrtgroup) & scdb$soilseriesstatus != "inactive" |
                      is.na(scdb$scdb_taxsuborder) & scdb$soilseriesstatus != "inactive",
                   c("soilseriesstatus","soilseriesname","taxclname","soiltaxclasslastupdated",
                     "taxorder","calc_taxorder",
                     "taxsuborder","calc_taxsuborder",
                     "taxgrtgroup","calc_taxgrtgroup",
                     "taxsubgrp","calc_taxsubgrp")]
nrow(fix_series)
table(fix_series$soilseriesstatus)

# read gNATSGO fGDB
gnatsgo <- st_read("E:/Geodata/soils/gNATSGO_CONUS.gdb", layer = "component")

# get components that have subgroup level taxonomy
gnatsgo_subgrp <- subset(gnatsgo, !is.na(taxsubgrp))

# use SoilTaxonomy package to look up subgroup letter codes
gnatsgo_subgrp$codes <- SoilTaxonomy::taxon_to_taxon_code(gnatsgo_subgrp$taxsubgrp)

# first we focus on the subgroup codes we cannot identify for components called "series"
na.ldx <- is.na(gnatsgo_subgrp$codes)
isseries.ldx <- gnatsgo_subgrp$compkind == "Series"

## Find components with obsolete, incomplete or incorrect taxonomy
# we will attempt to create an update table based on the [hopefully] more recent entries in SC database

# 60900 components have subgroup level taxon letter codes that do not exist in 12th edition look up table
sum(na.ldx)

# of those 60900, 54802 have a series-level component (based on compkind)
sum(na.ldx & isseries.ldx, na.rm = TRUE)

# in those 54802 series-level components, there are 4947 unique component names
badcode <- subset(gnatsgo_subgrp, na.ldx & isseries.ldx)
badcodeseries <- badcode$compname
uniquebadseries <- unique(badcodeseries)

# TODO: use data.table
test <- sapply(trimws(scdb$soilseriesname), 
               function(y) grepl(y, trimws(toupper(uniquebadseries))))
test2 <- apply(test, 2, which)

# 4668 of 4947 components with compkind series (and obsolete or missing taxonomy)
#  appear to [at least partially] match a series defined in the SC database
ers <- pmatch(trimws(toupper(scdb$soilseriesname)), 
              trimws(toupper(uniquebadseries)))
sum(!is.na(ers))

# inspect problematic names
unmatched.compnames <- uniquebadseries[!(1:length(uniquebadseries) %in% ers)]

"Geiger" %in% badcode$compname

# these need 
fix_components <- badcode[badcode$compname %in% unmatched.compnames,]

# these probably just have some garbage in the compname
update_components <- badcode[!badcode$compname %in% unmatched.compnames,]

# all failed to have a letter code assigned based on 12th ed taxonomy
all(is.na(fix_components$codes))
all(is.na(update_components$codes))

# obsolete taxonomy
unique(fix_components$taxsubgrp)
unique(update_components$taxsubgrp)

# several different vintages
unique(fix_components$soiltaxedition)
unique(update_components$soiltaxedition)

comp_series_lut <- trimws(toupper(unique(update_components$compname)))
lut_series_index <- pmatch(comp_series_lut, trimws(toupper(scdb$soilseriesname)))
names(lut_series_index) <- comp_series_lut
update_components$series_index <- lut_series_index[trimws(toupper(update_components$compname))]


# assign an integer index for position in NASIS SC db
scdb$series_index <- 1:nrow(scdb)
update_components <- data.table::merge.data.table(suffixes = c('.SSURGO','.SCDB'), 
                                             all.x = TRUE,
                                             update_components, 
                                             scdb, 
                                             by = "series_index")

# calculate a new subgroup code based on the subgroup from SC database (instead of component)
update_components$newcode <- SoilTaxonomy::taxon_to_taxon_code(update_components$taxsubgrp.SCDB)

# separate out series that appear to need update in SC db
update_series <- subset(update_components, is.na(newcode) | is.na(series_index))
update_series$compname

# 63 unique compnames have obsolete taxonomy in the SC database
unique(update_series$compname)

# keep those that are related to a "modern" series family class
update_components <- subset(update_components, !is.na(newcode))

# these components have an erroneous "twelfth edition" assigned to older taxonomy
subset(update_components, soiltaxedition == "twelfth edition")

# 53356 components could in theory be semi-auto-updated (using table of current series taxonomy)
sum(complete.cases(update_components[,c('taxclname.SSURGO','taxclname.SCDB', 
                                        "taxorder.SCDB", "taxsuborder.SCDB",
                                        "taxgrtgroup.SCDB", "taxsubgrp.SCDB",
                                        'newcode')]))

res0 <- fix_series
write.csv(as.data.frame(res0), file = "fix_series_SC_db.csv")

# these are not-properly-named components with compkind == 'Series' and old taxonomy
#  (cannot be related to a particular series that exists in SC)
res1 <- fix_components[,c('mukey','cokey','compname','compkind','codes','taxclname')]
write.csv(res1, "fix_SSURGO_components.csv")

# convert the decomposed codes back to taxonomy
update_components[["calc_taxorder"]] <- taxon_code_to_taxon(update_components[["calc_taxorder"]])
update_components[["calc_taxsuborder"]] <- taxon_code_to_taxon(update_components[["calc_taxsuborder"]])
update_components[["calc_taxgrtgroup"]] <- taxon_code_to_taxon(update_components[["calc_taxgrtgroup"]])
update_components[["calc_taxsubgrp"]] <- taxon_code_to_taxon(update_components[["calc_taxsubgrp"]])

# these can [in theory] be updated based on records in the SC DB
#  show the tax levels in SC DB alongside the "calculated" (based on decomposing the subgroup code)
res2 <- update_components[,c('mukey','cokey','compname','compkind','codes',
                'soilseriesname','taxclname.SSURGO','taxclname.SCDB',
                "taxorder.SCDB", "calc_taxorder",
                "taxsuborder.SCDB", "calc_taxsuborder",
                "taxgrtgroup.SCDB", "calc_taxgrtgroup",
                "taxsubgrp.SCDB", "calc_taxsubgrp")]
write.csv(res2, "update_SSURGO_components_from_SC.csv")

# compare higher taxa (no subgroup or lower level taxa defined)
gnatsgo$order_code <- taxon_to_taxon_code(gnatsgo$taxorder)
gnatsgo$suborder_code <- taxon_to_taxon_code(gnatsgo$taxsuborder)
gnatsgo$grtgroup_code <- taxon_to_taxon_code(gnatsgo$taxgrtgroup)
gnatsgo$subgroup_code <- taxon_to_taxon_code(gnatsgo$taxsubgrp)

ordermatch_suborder <- gnatsgo$order_code == substr(gnatsgo$suborder_code,1,1)
table(ordermatch_suborder, useNA = "ifany")

ordermatch_grtgroup <- gnatsgo$order_code == substr(gnatsgo$grtgroup_code,1,1)
table(ordermatch_grtgroup, useNA = "ifany")

subordermatch_grtgroup <- gnatsgo$suborder_code == substr(gnatsgo$grtgroup_code,1,2)
table(subordermatch_grtgroup, useNA = "ifany")

ordermatch_subgroup <- gnatsgo$order_code == substr(gnatsgo$subgroup_code,1,1)
table(ordermatch_subgroup, useNA = "ifany")
subordermatch_subgroup <- gnatsgo$suborder_code == substr(gnatsgo$subgroup_code,1,2)
table(subordermatch_subgroup, useNA = "ifany")
greatgroupmatch_subgroup <- gnatsgo$grtgroup_code == substr(gnatsgo$subgroup_code,1,3)
table(greatgroupmatch_subgroup, useNA = "ifany")

# TODO: compare the modified {SoilTaxonomy} lookup table against NASIS domains
#   - there are a few taxa defined in domains that do NOT exist in 12th edition


