#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- setup ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(tidyverse, taxize)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- load data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bd.raw <- vroom::vroom('Data y ambiente R diciembre 2023/BASE DE DATOS FLORA_MIRA MATAJE_2021.csv')

col.taxonomy <- 
  vroom::vroom('Data y ambiente R diciembre 2023/taxonomy_plants/col/Taxon.tsv')

wcvp.taxonomy <- 
  vroom::vroom('Data y ambiente R diciembre 2023/taxonomy_plants/wcvp_actualizado/wcvp_names.csv')

wcvp.distribution <- 
  vroom::vroom('Data y ambiente R diciembre 2023/taxonomy_plants/wcvp_actualizado/wcvp_distribution.csv')

wfo.taxonomy <- 
  vroom::vroom("C:/Users/rober/folder_mega/archivos_varios/taxonomy_plants/WFO_Backbone/classification.csv")

bd.raw.01 <- 
bd.raw %>% 
  filter(str_detect(str_to_lower(basisOfRecord), 'specimen' )) %>% 
  mutate(across(everything(), 
                ~str_replace(.x, '<Null>|(s|S)in (d|D)atos|(U|u)nknown|^s\\.n\\.$', NA_character_)),
         scientificName = 
           str_replace_all(scientificName, '\\s', ' ') %>% 
           str_replace('fo\\.', 'f\\.') %>% 
           str_squish() %>% 
           str_to_sentence(),
         recordedBy = 
           recordedBy %>% 
           stringi::stri_trans_general("Latin-ASCII") %>% 
           str_replace_all(',|;', '\\|') %>% 
           str_remove_all('\\.|,|;|\\)|\\(||\\?|\\*') %>% 
           str_squish(),
         locality =
           str_to_lower(locality) %>% 
           stringi::stri_trans_general("Latin-ASCII") %>% 
           str_remove_all('\\.|,|;|\\)|\\(') %>% 
           str_squish()) %>% 
  nest(data = -scientificName) %>% 
  mutate(id.scientificName = 1:n()) %>% 
  unnest(data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- harmonize scientificNames ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## WCVP names
old.names <- 
bd.raw.01 %>% 
  distinct(scientificName, id.scientificName) %>% 
  rename(old.scientificName = scientificName ) %>% 
  filter(!is.na(old.scientificName)) %>% 
  mutate(taxonomy.checked = NA_character_)

## Left join with wcvp database to obtain accepted names
new.names.wcvp <- 
old.names %>% 
  left_join(wcvp.taxonomy, by = c('old.scientificName' = 'taxon_name'))

## Names without change 
names.withoutChange.wcvp <- 
new.names.wcvp %>% 
  filter(taxon_status %in% c('Accepted')) %>% 
  select(id.scientificName, 
         family, 
         scientificName = old.scientificName, 
         scientificNameAuthorship = taxon_authors,
         taxonRank = taxon_rank,
         genus, 
         specificEpithet = species,
         infraspecificEpithet = infraspecies) %>% 
  mutate(taxonomy.checked = 
           'wcvp')

## count number of accepted names in each taxon_status 
new.names.wcvp %>%
  count(taxon_status, is.na(accepted_plant_name_id))

## priorize taxon_status synonym
new.names.wcvp.01 <- 
new.names.wcvp %>% 
  filter(!taxon_status %in% c('Accepted', 'Unplaced'), !is.na(taxon_status)) %>% 
  mutate(
    prior.taxon_status = 
      as.factor(taxon_status) %>% 
      fct_relevel('Synonym', 'Illegitimate', 'Misapplied', 'Invalid')
  )
    
new.names.wcvp.02 <- 
  new.names.wcvp.01 %>% 
  arrange(prior.taxon_status) %>% 
  distinct(id.scientificName, .keep_all = T) %>% 
  select(id.scientificName, old.scientificName,
         accepted_plant_name_id) %>% 
  left_join(wcvp.taxonomy, by = c('accepted_plant_name_id' = 'plant_name_id')) %>% 
  # relocate(taxon_name, .before = old.scientificName) %>% 
  select(id.scientificName, 
         family, 
         scientificName = taxon_name, 
         scientificNameAuthorship = taxon_authors,
         taxonRank = taxon_rank,
         genus, 
         specificEpithet = species,
         infraspecificEpithet = infraspecies) %>% 
  mutate(taxonomy.checked = 'wcvp')
  
new.names.wcvp.03 <- 
names.withoutChange.wcvp %>% 
  bind_rows(new.names.wcvp.02) %>% 
  distinct(id.scientificName, .keep_all = T)

## update names
updated.names <- 
old.names %>% 
  select(id.scientificName) %>% 
  left_join(new.names.wcvp.03,
            by = 'id.scientificName')
  
## Names not available in wcvp database
old.names.01 <- 
updated.names %>% 
  filter(is.na(taxonomy.checked)) %>% 
  select(id.scientificName) %>% 
  left_join(old.names)

## WFO names
new.names.wfo <- 
old.names.01 %>% 
  filter(!str_detect(old.scientificName, ' (sect|subg)\\.? ')) %>% 
  left_join(wfo.taxonomy, by = c('old.scientificName' = 'scientificName'))

## names without change
names.withoutChange.wfo <- 
  new.names.wfo %>% 
  filter(taxonomicStatus %in% c('Accepted', 'Unchecked')) %>% 
  select(id.scientificName, 
         family, 
         scientificName = old.scientificName, 
         scientificNameAuthorship,
         taxonRank,
         genus, 
         specificEpithet,
         infraspecificEpithet) %>% 
  mutate(taxonomy.checked = 
           'wfo')

## names synonyms
new.names.wfo.01 <- 
new.names.wfo %>% 
  filter(taxonomicStatus %in% c('Synonym')) %>% 
  select(id.scientificName, acceptedNameUsageID) %>% 
  left_join(wfo.taxonomy, by = c('acceptedNameUsageID' = 'taxonID')) %>% 
  select(id.scientificName, 
         family, 
         scientificName, 
         scientificNameAuthorship,
         taxonRank,
         genus, 
         specificEpithet,
         infraspecificEpithet) %>% 
  mutate(taxonomy.checked = 
           'wfo')

## unite accepted names
new.names.wfo.02 <- 
names.withoutChange.wfo %>% 
  bind_rows(new.names.wfo.01) %>% 
  distinct(id.scientificName, .keep_all = T)

## update names with wfo base
updated.names.01 <- 
updated.names %>% 
  rows_update(new.names.wfo.02, by = 'id.scientificName')

## Tropicos names
old.names.02 <- 
  updated.names.01 %>% 
  filter(is.na(taxonomy.checked)) %>% 
  select(id.scientificName) %>% 
  left_join(old.names) %>% 
  filter(!str_detect(old.scientificName, ' (sect|subg)\\.? '))

key.tp <- '74b997b6-c313-4158-b64d-7b02450047a8'

new.names.tp <- 
old.names.02 %>% 
  mutate(id.tropicos = 
           map(
             old.scientificName,
             ~tryCatch(
               {tp_search(.x, key = key.tp)}
             )
           ))

new.names.tp.01 <- 
new.names.tp %>% 
  unnest() %>% 
  filter(is.na(error)) %>% 
  mutate(id.prior = 
          as.factor(nomenclaturestatusname) %>% 
           fct_relevel('Legitimate', 'No opinion', 'Illegitimate')) %>% 
  arrange(id.prior) %>% 
  distinct(id.scientificName, .keep_all = T)

new.names.tp.02 <- 
new.names.tp.01 %>% 
  mutate(accnames = 
           map(
             nameid,
             ~tryCatch(
               {tp_accnames(.x, key = key.tp)}
             )
           ))

names.withoutChange.tp <- 
new.names.tp.02 %>% 
  unnest() %>%
  distinct(id.scientificName, .keep_all = T) %>% 
  filter(str_detect(accnames, 'No accepted names')) %>% 
  select(id.scientificName, nameid)

new.names.tp.03  <- 
  new.names.tp.02 %>% 
  unnest() %>%
  distinct(id.scientificName, .keep_all = T) %>% 
  filter(!str_detect(accnames, 'No accepted names')) %>% 
  select(id.scientificName, accnames) %>% 
  unnest(accnames) %>% 
  select(id.scientificName, nameid)

new.names.tp.04 <- 
  bind_rows(names.withoutChange.tp, new.names.tp.03) %>% 
  distinct(.keep_all = T)

## search taxonomic information in tropicos
dummy.tibble <- 
  tibble(author = NA_character_,
         scientificname = NA_character_,
         family= NA_character_,
         taxonRank= NA_character_,
         genus= NA_character_,
         specificEpithet= NA_character_,
         infraspecificEpithet = NA_character_)

new.names.tp.05 <- 
new.names.tp.04 %>% 
  mutate(taxon.info =
           map(
             nameid,
             ~tryCatch(
               {tp_summary(.x, key = key.tp) %>% 
                   bind_rows(dummy.tibble[0,])
                   }
             )
           ))

new.names.tp.06 <- 
new.names.tp.05 %>% 
  unnest() %>% 
  select(id.scientificName,
         scientificName = scientificname,
         scientificNameAuthorship = author,
         family,
         taxonRank = rank,
         genus,
         specificEpithet = speciesepithet,
         infraspecificEpithet = otherepithet
  ) %>% 
  mutate(taxonomy.checked = 'tropicos')

## update names with tropicos
updated.names.02 <- 
  updated.names.01 %>% 
  rows_update(new.names.tp.06, by = 'id.scientificName')

## Eliminate sub-taxons (section, var, subsp) and harmonize name of
## higher level taxon
old.names.03 <-
  updated.names.02 %>% 
  filter(is.na(taxonomy.checked)) %>% 
  select(id.scientificName) %>% 
  left_join(old.names) %>% 
  mutate(old.scientificName =
           str_remove(old.scientificName, ' [A-z]+\\. [A-z-]+$'))

new.names.residual <- 
old.names.03 %>% 
  left_join(wcvp.taxonomy, by = c('old.scientificName' = 'taxon_name'))

names.withoutChange.residual <- 
new.names.residual %>% 
  filter(taxon_status %in% c('Accepted')) %>% 
  select(id.scientificName, 
         family, 
         scientificName = old.scientificName, 
         scientificNameAuthorship = taxon_authors,
         taxonRank = taxon_rank,
         genus, 
         specificEpithet = species,
         infraspecificEpithet = infraspecies) %>% 
  mutate(taxonomy.checked = 
           'wcvp')

new.names.residual.01 <- 
new.names.residual %>% 
  filter(!taxon_status %in% c('Accepted'),
         !is.na(plant_name_id)) %>% 
  select(id.scientificName, accepted_plant_name_id) %>% 
  left_join(wcvp.taxonomy, by = c('accepted_plant_name_id' = 'plant_name_id')) %>% 
  select(id.scientificName, 
         family, 
         scientificName = taxon_name, 
         scientificNameAuthorship = taxon_authors,
         taxonRank = taxon_rank,
         genus, 
         specificEpithet = species,
         infraspecificEpithet = infraspecies) %>% 
  mutate(taxonomy.checked = 
           'wcvp')

new.names.residual.02 <- 
  bind_rows(
    names.withoutChange.residual,
    new.names.residual.01
  ) %>% 
  distinct(id.scientificName, .keep_all = T)

## Last update of residual names
updated.names.03 <- 
  updated.names.02 %>% 
  rows_update(new.names.residual.02, by = 'id.scientificName') %>% 
  mutate(taxonRank = str_to_lower(taxonRank))

## export to visualize new names
# updated.names.03 %>% 
#   write_csv('Data y ambiente R diciembre 2023/intermediate_files/updated.names.03.csv', na='')

## Update database with new scientificnames
bd.raw.02 <- 
bd.raw.01 %>% 
  rows_update(updated.names.03 %>% 
                select(-taxonomy.checked),
              by = 'id.scientificName') %>% 
  nest(data=-recordedBy) %>% 
  mutate(id.recordedBy = 1:n()) %>% 
  unnest() %>% 
  mutate(taxonRank = str_to_lower(taxonRank))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- harmonize collector names  ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Add a key for each name
collector.names.01 <- 
bd.raw.02 %>% 
  mutate(recordedBy.01 = 
           case_when(
             str_detect(recordedBy, '\\||&') ~
               str_extract(recordedBy, '^.+?(?>(\\||&))') %>% 
               str_remove('\\|') %>% 
               str_squish(),
             TRUE ~ recordedBy
           )
           ) %>% 
  count(id.recordedBy, recordedBy, recordedBy.01) %>% 
  mutate(n.space = 
           str_count(recordedBy.01, ' '),
         extracted =
           case_when(
             n.space == 0 ~
               recordedBy.01,
             
             n.space == 1 ~ 
               str_extract(recordedBy.01, '[A-z-]+$'),
             
             n.space == 2 ~ 
               str_extract(recordedBy.01, '[A-z-]+$'),
             
             n.space == 3 ~ 
               str_extract(recordedBy.01, '(?<=^[A-z-]{1,100} [A-z-]{1,100} )[A-z-]+'),
             
             TRUE ~ recordedBy.01),
         n.letters = str_length(extracted))

## obtain surnames of irregular names
collector.names.02 <- 
collector.names.01 %>% 
  mutate(extracted =
           case_when(
             str_detect(recordedBy.01, '(v|V)an (d)er (w|W)erff|HHvd Werff|Hvd Werff') ~ 
               'van der Werff',
             
             str_detect(recordedBy.01, 'Olga S de Benavides') ~ 
               'Benavides',
             
             str_detect(recordedBy.01, 'H GarcA-a Barriga|H Garcia Barriga') ~ 
               'Garcia-Barriga',
             
             str_detect(recordedBy.01, 'Jansen B') ~ 
               'Jansen',
             
             str_detect(recordedBy.01, 'Alice D A Fay') ~ 
               'Fay',
             
             str_detect(recordedBy.01, 'Lorena Endara A') ~ 
               'Endara',
             
             str_detect(recordedBy.01, 'E Asplund Coll Delessert 1948') ~ 
               'Coll',
             
             str_detect(recordedBy.01, 'Alina Freire Fierro|A Freire Fierro|A Freire-Fierro|Alina Freire') ~ 
               'Freire-Fierro',
             
             str_detect(recordedBy.01, 'Mats H G Gustafsson') ~ 
               'Gustafsson',
             
             str_detect(recordedBy.01, 'Edgar Gudino Jara') ~ 
               'Gudino',
             
             str_detect(recordedBy.01, 'Ynes Enriquetta Julietta Mexia') ~ 
               'Mexia',
             
             str_detect(recordedBy.01, 'Homero Vargas L') ~ 
               'Vargas',
             
             str_detect(recordedBy.01, 'Susana Leon|S Leon|S Leon-Yanez|Susana Leon Yanez') ~ 
               'Leon-Yanez',
             
             str_detect(recordedBy.01, 'Linda K Albert de Escobar') ~ 
               'Escobar',
             
             str_detect(recordedBy.01, 'Moscol Olivera M') ~ 
               'Olivera',
             
             str_detect(recordedBy.01, 'Maria del Carmen Ulloa Ulloa|Carmen Ulloa U') ~ 
               'Ulloa',
             
             str_detect(recordedBy.01, 'A Grijalva P') ~ 
               'Grijalva',
             
             str_detect(recordedBy.01, 'Xavier Cornejo S') ~ 
               'Cornejo',
             
             str_detect(recordedBy.01, 'M Acosta Solis') ~ 
               'Acosta-Solis',
             
             str_detect(recordedBy.01, 'P Mena V') ~ 
               'Mena',
             
             str_detect(recordedBy.01, 'Kenneth A Wilson') ~ 
               'Wilson',
             
             str_detect(recordedBy.01, 'Hoover Ws') ~ 
               'Hoover',
             
             str_detect(recordedBy.01, 'Jorgesen|Jorgenson') ~ 
               'Jorgensen',
             
             str_detect(recordedBy.01, 'Paola Pedraza') ~ 
               'Pedraza-Penalosa',
             
             str_detect(recordedBy.01, 'Gavilanes') ~ 
               'Gavilanez',
             
             TRUE ~ extracted
           ),
         extracted = str_to_title(extracted)) %>% 
  filter(!is.na(extracted))


collector.names.03 <- 
collector.names.02 %>% 
  select(id.recordedBy, recordedBy.extracted = extracted)

## Update database with extracted collector surnames
bd.raw.03 <- 
bd.raw.02 %>% 
  left_join(collector.names.03, by = 'id.recordedBy') %>% 
  mutate(decimalLongitude =
           if_else(!str_detect(decimalLongitude, '\\.'),
                   str_replace(decimalLongitude, '(\\d{2})(\\d*)', '\\1\\.\\2'),
                   decimalLongitude
                 )
         )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- filter duplicates ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# potential filters: eventDate, locality, decimalLongitude, decimalLatitude
# minimumElevationInMeters, recordedBy, recordNumber

# bd.raw.03 %>% 
#   filter(is.na(recordedBy)) %>% 
#   select(family, scientificName, eventDate, locality, decimalLongitude, decimalLatitude,
#          minimumElevationInMeters, recordedBy, recordNumber) %>% 
#   view()

## 42,446 with recordNumber and recordedBy
bd.raw.04 <- 
bd.raw.03 %>% 
  filter(!is.na(recordedBy)&!is.na(recordNumber)) %>% 
  mutate(id.duplicate.family = 
           paste(family, recordedBy.extracted, recordNumber),
         
         id.duplicate.decimalLongitude = 
         paste(round(parse_number(decimalLongitude), digits = 2), 
               recordedBy.extracted, 
               recordNumber),
         
         id.duplicate =
           paste(recordedBy.extracted, recordNumber),
         
         duplicate.family = duplicated(id.duplicate.family),
         duplicate.decimalLongitude = duplicated(id.duplicate.decimalLongitude),
         duplicate = duplicated(id.duplicate),
         
         identifiedBy =
           str_replace(identifiedBy, '\\*', NA_character_)
         )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- contrast recordedBy and recordNumber in different context 
# (family and geographic)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## False negatives rate
bd.raw.04 %>% 
  count(duplicate, duplicate.decimalLongitude)


bd.raw.04 %>% 
  count(duplicate, duplicate.family)

## using decimalLongitude as context
contrast.collectorVsdecimalLongitudeContext <- 
bd.raw.04 %>% 
  filter(duplicate != duplicate.decimalLongitude) 

duplicates.not.detected.withDecimalLongitude <- 
bd.raw.04 %>% 
  select(occid,
         recordedBy, recordedBy.extracted, recordNumber, 
         duplicate, duplicate.family, duplicate.decimalLongitude, 
         family, 
         scientificName, 
         decimalLongitude, 
         locality,
         institutionCode) %>% 
  arrange(recordedBy.extracted, recordNumber) %>% 
  semi_join(contrast.collectorVsdecimalLongitudeContext, 
            by = c('recordedBy.extracted', 'recordNumber'))

duplicates.not.detected.withDecimalLongitude %>% 
  count(duplicate, duplicate.decimalLongitude)

## using family as context
contrast.collectorVsfamilyContext <- 
  bd.raw.04 %>% 
  filter(duplicate != duplicate.family) 

duplicates.not.detected.withFamily <- 
bd.raw.04 %>% 
  select(occid,
         recordedBy, recordedBy.extracted, recordNumber, 
         duplicate, duplicate.family, duplicate.decimalLongitude, 
         family, 
         scientificName, 
         decimalLongitude, 
         locality,
         institutionCode) %>% 
  arrange(recordedBy.extracted, recordNumber) %>% 
  semi_join(contrast.collectorVsfamilyContext, 
            by = c('recordedBy.extracted', 'recordNumber')) 

duplicates.not.detected.withFamily %>% 
  count(duplicate, duplicate.family)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- prioritize duplicates with better taxonomic treatment ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## order by:
## identififiedBy 
## Herbarium
## taxonRank : species, genus, family
bd.raw.05 <- 
bd.raw.04 %>% 
  mutate(
    prior.identifiedBy =
      if_else(!is.na(identifiedBy), 1, 2),
    
    prior.Herbarium =
      as.factor(institutionCode) %>% 
      fct_relevel(
        'MO', 'NY', 'US', 'F', 
        'INABIOEC',
        'QCA', 'AAU'),
    
    prior.taxonRank =
      as.factor(taxonRank) %>% 
      fct_relevel(
        'subspecies', 'variety', 'species', 'genus', 'family'
      )
  )

# remove duplicates
db.without.duplicates <- 
bd.raw.05 %>% 
  arrange(prior.identifiedBy, prior.Herbarium, prior.taxonRank) %>% 
  distinct(id.duplicate, .keep_all = T)


# db.without.duplicates %>% 
#   write_csv('export/db_flora_without_duplicates.csv', na = '')


db.without.duplicates %>% 
  # count(taxonRank)
  filter(taxonRank == 'species') %>% 
  distinct(scientificName)











