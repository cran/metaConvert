## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE, message=FALSE, warning=FALSE, results='hide'-----------------
library(metaConvert)
library(DT)

## ----eval = FALSE-------------------------------------------------------------
#  see_input_data(measure = "or")

## ----echo=FALSE---------------------------------------------------------------
DT::datatable(see_input_data(measure = "or", extension="data.frame"), options = list(  
    scrollX = TRUE,
    dom = c('t'),
    scrollY = "300px", 
    pageLength = 500,
    ordering = FALSE,
    rownames = FALSE,
    columnDefs = list(
                  list(width = '130px',
                       targets = "_all"),
                  list(className = 'dt-center', 
                                     targets = "_all"))))

## ----eval = FALSE-------------------------------------------------------------
#  data_extraction_sheet(measure = "or")

## ----echo=FALSE---------------------------------------------------------------
DT::datatable(data_extraction_sheet(measure = "or", extension="data.frame"), options = list(  
    scrollX = TRUE,
    dom = c('t'),
   ordering = FALSE,
    rownames = FALSE,
    columnDefs = list(
                  list(width = '130px',
                       targets = "_all"),
                  list(className = 'dt-center', 
                                     targets = "_all"))))

## ----message=FALSE, fig.width= 11---------------------------------------------
compare_df(
    df_extractor_1 = df.compare1,
    df_extractor_2 = df.compare2,
    output = "html")

## ----message=FALSE------------------------------------------------------------
res = convert_df(x = df.short, measure = "g")

## ----eval= FALSE--------------------------------------------------------------
#  summary(res)

## ----echo=FALSE---------------------------------------------------------------
DT::datatable(summary(res), options = list(  
    scrollX = TRUE,
    dom = c('t'),
    ordering = FALSE,
    rownames = FALSE,
    scrollY = "300px", 
    pageLength = 100,
    columnDefs = list(
                  list(width = '130px',
                       targets = "_all"),
                  list(className = 'dt-center', 
                                     targets = "_all"))))

