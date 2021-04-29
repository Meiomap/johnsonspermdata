## code to prepare `DATASET` dataset goes here
library(GEOquery)
library(tidyverse)

get_and_folder_in <- function(x,download=TRUE)
{
  metadata=data.frame()
  serverdata=getGEOSuppFiles(x,fetch_files = FALSE,makeDirectory=FALSE,filter_regex='idat.gz')
  if (is.null(serverdata))
  {
    Sys.sleep(5)
    serverdata=getGEOSuppFiles(x,fetch_files = FALSE,makeDirectory=FALSE,filter_regex='idat.gz')
  }
  if (!is.null(serverdata))
  {
    metadata= serverdata %>%
      separate(fname,into = c('sampleid','arrayid','position','channel'),sep='_',remove=FALSE) %>%
      mutate(outfile=str_extract(fname,'[0-9]+_R[0-9]+C[0-9]+_(Grn|Red).idat'))

    if (download==TRUE)
    {
      metadata %>%
        select(arrayid) %>% as.matrix() %>%
        map(function(y) dir.create(path=as.character(y),showWarnings=FALSE))
      metadata %>%
        select(arrayid,fname,sampleid)  %>%
        pmap(function(arrayid,fname,sampleid) getGEOSuppFiles(sampleid,fetch_files = TRUE,makeDirectory = FALSE,baseDir=arrayid,filter_regex=fname))
      #map2(fnamefunction(x)getGEOSuppFiles(x,fetch_files = TRUE,makeDirectory=TRUE,filter_regex=z))


      metadata %>%
        select(sampleid,arrayid,fname,outfile) %>%
        pmap(function(sampleid,arrayid,fname,outfile) gunzip(file.path(arrayid,fname),destname = file.path(arrayid,outfile),overwrite=TRUE,remove=FALSE))
    }
  }
  return(metadata)
}


write_samplesheet <- function(filename,df,gtcpath)
{
  header=c('[Header],,,,',
           'INVESTIGATOR NAME,Ivan,,,',
           'PROJECT NAME,Test1,,,',
           'EXPERIMENT NAME,SingleCellReference,,,',
           'DATE,19/10/2016 00:00,,,',
           '[Manifests],,,,',
           'A,HumanKaryomap-12v1_A,,,',
           '[Data],,,,')

  writeLines(header, filename)


  df_1 = df %>%
    select(sampleid,arrayid,position)  %>%
    distinct() %>%
    mutate(Path=gtcpath,
           Aux=0)
  colnames(df_1)=c('Sample_ID','SentrixBarcode_A','SentrixPosition_A','Path','Aux')
  df_1 %>%
    write.table(file=filename,append=TRUE,row.names=FALSE,quote=FALSE,sep=',')
}


gse <- getGEO("GSE19247",GSEMatrix = TRUE)
samplelist_sperm=as.data.frame(gse$`GSE19247-GPL6985_series_matrix.txt.gz`) %>%
  filter((str_detect(cell.type.ch1,'sperm')) & str_detect(cell.amplication.ch1,'MDA')) %>%
  rownames()


metadata=as.data.frame(gse$`GSE19247-GPL6985_series_matrix.txt.gz`) %>%
  filter((str_detect(cell.type.ch1,'sperm')) & str_detect(cell.amplication.ch1,'MDA'))%>%
  select(title,source_name_ch1,cell.amplication.ch1) %>%
  mutate(familyid=gsub('.*family\ ([0-9]+)$','\\1',source_name_ch1),
         embryoid=gsub('.*embryo\ ([0-9]+)\ from\ family [0-9]+','\\1',source_name_ch1),
         celltype=gsub('single\ ([a-z ]+)\ amplified by MDA','\\1',cell.amplication.ch1),
         id=rownames(.) %>% str_to_lower()) %>%
  mutate(embryoid=gsub('.*family\ ([0-9]+)$','\\1',embryoid),)  %>%
  select(id,title,familyid,embryoid,celltype)

metadata[metadata$celltype=='sperm cell',]$embryoid=metadata[metadata$celltype=='sperm cell',]$title
metadf_sperm=map_if(samplelist_sperm,function(x) !is.null(x),function(x) get_and_folder_in(x,download=FALSE))
metadf_sperm_merged = Reduce(function(...) merge(..., all=T), metadf_sperm)
write_samplesheet('data/johnsonsperm.csv',metadf_sperm_merged,gtcpath = './' )

samplesheet=file.path('johnsonsperm.csv')

usethis::use_data(metadata, overwrite = TRUE)
usethis::use_data(samplesheet,overwrite=TRUE)
