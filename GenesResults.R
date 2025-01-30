#######################################################
#
# Explorando dataframe para valores de expresión génica provenientes de RNA-Seq
# 
#######################################################
library(dplyr)
#######################################################
data <-  read.table("Databases/Tx1.rep1.genes.results" , sep = "\t", header = T)
names(data)
class(data) # clase de objeto
dim(data) #obtener las dimensiones del objeto 
ncol(data) # numero de columnas
nrow(data) # numero de genes en el data frame

---
  # Conocer cual es el valor mínimo y máximo de una variable:
  
min(data$expected_count)
max(data$expected_count)
---
  # Identificar el gen de mayor longitud:
  
max.Length <-row.names(data[which.max(data$length), ] )
x.Length <-row.names(data)[data$length=="615"]
y.Length <- (data)[data$length > 10000, ] 
---
  # Función subsets:
  
  y.Length <- (data)[data$length > 10000, ]
z.Length <- subset(data, data$length > 10000 ) 
max.Counts <- subset(data, data$expected_count== max(data$expected_count))

---
# Obtener un subset de genes con un número de conteos determinado:
  
threshold.1 <- subset(data, data$expected_count > 1000)
length(threshold.1$transcript_id.s.)

---
  # Es tu turno:
  # Obten genes que con una longitud mayor a 1000 pb y que tegan más de 50 conteos asignados por RSEM.
  
  ---
############################################
#  dplyr::filter():
#############################################
library(dplyr)

data %>% 
  filter(expected_count > 1000) %>%
  count()

---
############################################
#  dplyr::summarize()
#############################################
#  ¿Cómo se distribuyen los valores de TPM?
data %>%
  summarise(TOTAL = n() , MEAN= mean(TPM) , SD= sd(TPM), VAR= var(TPM))
%>% datawizard::data_rotate(colnames = T, rownames = "TOTAL")


---
  #  Tu turno:
  # Crea un nuevo objeto donde puedas obtener genes con valores de FPKM > 15 considerados como "sobrerepresetnados"; y valores de FPKM < 15 considerados como "reprimidos".
  
  OverExpression <-   
  
  DownExpression <-
  ---
  ############################################
#  dplyr::mutate()
#############################################
# Transformación logaritmica de los conteos observados para cada gen:
data <- data %>% 
  mutate(LOG = log10(expected_count))
data$LOG

############################################
#  dplyr::select()
#############################################

data %>% 
  select(transcript_id.s.,length, expected_count)%>% 
  head()
---
  # ¿Podemos combinar funciones dplyr?
  #  filter() %>% select() %>% mutate() 
  
  ---
  # Obtener unigenes con valores de FPKM mayores a 5, y con una longitud mayor a 1000 pb 
  # Abundancia Relativa= Valor en cada celda/ Suma total de la columna
  
  AbRel <-  OverExpression %>% 
  select(FPKM, length) %>% 
  filter(FPKM > 5 & length >1000) %>% 
  mutate(AbunRel= FPKM/sum(data$FPKM))

