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
averageTx1 <-  data.clean %>%
  head() %>%
  rowwise() %>% # Trabajar fila por fila
  mutate(averageTx1 = mean(c_across(2:5))) 
dim(averageTx1)
#######################################################
# 2.2 Obtener promedio por réplicas de cada tejido
#######################################################
data.clean <-   data.clean %>%
  rowwise() %>% # Trabajar fila por fila
  mutate(averageTx1 = mean(c_across(2:4))) %>%
  mutate(averageTx2 = mean(c_across(5:7))) %>%
  mutate(averageTx3 = mean(c_across(8:10))) %>%
  mutate(averageControl = mean(c_across(11:13))) %>%
  ungroup() # Liberar el modo rowwise
names(data.clean)

dim(data.clean)
data.clean$
#######################################################
# Crear un dataframe unicamente con las columnas de promedios entre réplicas
#######################################################
data.mean <- data.clean %>% 
  select(GeneID, averageTx1, averageTx2, averageTx3, averageControl) %>% 
names(data.mean)
dim(data.mean)
summary(data.mean)

#######################################################
# Filtrar unigenes con al menos 10 cuentas mapeadas en cada experimento
#######################################################
c10.Tx1 <-  data.mean %>% filter(averageTx1 > 10) 
c10.Tx2 <-  data.mean %>% filter(averageTx2 > 10)    
c10.Tx3 <-  data.mean %>% filter(averageTx3 > 10) 
c10.C <-  data.mean %>% filter(averageControl > 10)

c10.number <- data.frame(
  "Tx1" = nrow(c10.Tx1), 
  "Tx2" = nrow(c10.Tx2),
  "Tx3" = nrow(c10.Tx3),
  "Control" = nrow(c10.C))


#######################################################
# Obtener genes con conteos mayoritarios para cada tejido(condición)
#######################################################
cm.Tx1 <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTx1 > c(averageTx2 & averageTx3& averageControl ))) %>% #verifica si averageTx1 es mayor que todas las otras columnas en la fila
  ungroup() 

cm.Tx2 <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTx2 > c(averageTx1 & averageTx3& averageControl ))) %>%
  ungroup()

cm.Tx3 <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTx3 > c(averageTx2 & averageTx1& averageControl ))) %>%
  ungroup()
 

cm.C <-  data.mean %>%
  rowwise() %>%  # Evalúa fila por fila
  filter( all(averageControl > c(averageTx2, averageTx3, averageTx1)))%>%  
  ungroup()


# Crear un vector para almacenar los resultados
cm.Genes <- data.frame(
  Control = nrow(cm.C),
  Tx1 = nrow(cm.Tx1),
  Tx2 = nrow(cm.Tx2),
  Tx3 = nrow(cm.Tx3) )

###############################
# rename()
##############################
Renombrar columna en el dataframe
data.mean <- rename(data.mean, PPM3 = averageTx1 ) # renombrer una columna
names(data.mean)
data.mean <- data.mean %>% # renombrer varias columnas, asignando el tipo o concentración del tratamiento experimental
  rename("PPM5" = "averageTx2",
         "PPM10" = "averageTx3")

############################################################
# Obtener foldchange Control vs. Tratamiento1
############################################################
# FC = cuentas gen1 (Tx1) / cuentas gen1 (Control)
############################################################
data.mean %>% 
mutate(averageControl = ifelse(is.infinite(averageControl), min(data.mean$averageControl), averageControl))%>%
  head()

data.mean <- data.mean %>% 
  rowwise() %>% # Trabajar fila por fila
  mutate(FC = averageTx1 /averageControl) %>% 
  select(GeneID,averageTx1, averageControl, FC) 
dim(data.mean)
boxplot(log10(data.mean$averageControl), log10(data.mean$averageTx1))

FCmax <-  data.mean %>%
   filter(FC > 2) %>%
  mutate(FC = ifelse(is.infinite(FC), NA, FC)) %>%
  filter(!is.na(FC))

FCmin <-  data.mean %>%
  filter(FC < 2)  %>%
  mutate(FC = ifelse(is.infinite(FC), NA, FC)) %>%
  filter(!is.na(FC))
  

