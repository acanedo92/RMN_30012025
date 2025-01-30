#######################################################
#
# Explorando matriz de abundancias generada por RSEM 
# para realizar preparar el análisis de expresión génica con
# datos provenientes de RNA-Seq
# 
#######################################################
library(dplyr)
#######################################################
data <-  read.table("Databases/AbundanceMatrix.isoform.counts.matrix" , sep = "\t", header = T)
names(data)[1] <-"GeneID"
names(data)[2:4] <-c("Control", "Tratamiento.1", "Tratamiento.2" )
names(data)
class(data)
head(data)
ncol(data)
dim(data)
#######################################################
#  Remover  unigenes con 0 cuentas en todas las columnas (condiciones)
#######################################################
data.clean <-  data %>%
  filter(!if_all(-GeneID, ~ . == 0))   

#######################################################
# 2.1 Obtener promedio por réplicas para Flores
#######################################################
averageFlores <-  data.clean %>%
  head() %>%
  rowwise() %>% # Trabajar fila por fila
  mutate(averageFlores = mean(c_across(2:5))) 
dim(averageFlores)
#######################################################
# 2.2 Obtener promedio por réplicas de cada tejido
#######################################################
data.clean <-   data.clean %>%
  rowwise() %>% # Trabajar fila por fila
  mutate(averageFlores = mean(c_across(2:4))) %>%
  mutate(averageHojas = mean(c_across(5:7))) %>%
  mutate(averageTallos = mean(c_across(8:10))) %>%
  mutate(averageVainas = mean(c_across(11:13))) %>%
  ungroup() # Liberar el modo rowwise

dim(data.clean)
#######################################################
# Crear un dataframe unicamente con las columnas de promedios entre réplicas
#######################################################
data.mean <- data.clean %>% 
  select(GeneID, averageFlores, averageHojas, averageTallos, averageVainas)
dim(data.mean)
summary(data.mean)
#######################################################
# Filtrar unigenes con al menos 10 cuentas mapeadas en cada experimento
#######################################################
c10.F <-  data.mean %>% filter(averageFlores > 10) 
c10.H <-  data.mean %>% filter(averageHojas > 10)    
c10.T <-  data.mean %>% filter(averageTallos > 10) 
c10.V <-  data.mean %>% filter(averageVainas > 10)

c10.number <- data.frame(
  "Flores" = nrow(value.F), 
  "Hojas" = nrow(value.H),
  "Tallos" = nrow(value.T),
  "Vainas" = nrow(value.V))


#######################################################
# Obtener genes con conteos mayoritarios para cada tejido(condición)
#######################################################
cm.F <-  data.mean %>%
  rowwise() %>%
  filter(all(averageFlores > c(averageHojas & averageTallos& averageVainas ))) %>% #verifica si averageFlores es mayor que todas las otras columnas en la fila
  ungroup() 

cm.H <-  data.mean %>%
  rowwise() %>%
  filter(all(averageHojas > c(averageFlores & averageTallos& averageVainas ))) %>%
  ungroup()

cm.T <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTallos > c(averageHojas & averageFlores& averageVainas ))) %>%
  ungroup()
 

cm.V <-  data.mean %>%
  rowwise() %>%  # Evalúa fila por fila
  filter( all(averageVainas > c(averageHojas, averageTallos, averageFlores)))%>%  
  ungroup()

# Crear un vector para almacenar los resultados
cm.Genes <- data.frame(
  Vainas = nrow(cm.V),
  Tallos = nrow(cm.T),
  Hojas = nrow(cm.H),
  Flores = nrow(cm.F) )
 


