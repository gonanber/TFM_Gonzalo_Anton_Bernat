![](https://capsule-render.vercel.app/api?type=waving&height=300&color=gradient&&customColorList=3&text=TFM%20Gonzalo%20Antón%20Bernat&fontSize=60&textBg=false&desc=Identificación%20de%20genes%20específicos%20del%20sexo%20en%20el%20Trastorno%20del%20Espectro%20Autista&descAlignY=60&fontAlignY=30&reversal=false)

## Índice

- [1_Microarrays](#1_Microarrays)
- [2_RNA-seq](#2_RNA-seq)
- [3_Genera_todos_fit](#3_Genera_todos_fit)
- [4_Metaanalisis](#4_Metaanalisis)
- [5_Analisis_funcional](#5_Analisis_funcional)
- [6_Material_suplementario](#6_Material_suplementario)


## 0_Descarga_datos

- En esta carpeta se encuentra el código necesario para descargar los estudios de microarrays y RNA-seq
- Debe ser lo primero en ser ejecutado


## 1_Microarrays

- En esta carpeta se encuentra el código necesario para ejecutar todo el análisis de los estudios de microarrays
- Debe ser lo segundo en ser ejecutado


## 2_RNA-seq

- En esta carpeta se encuentra el código necesario para ejecutar todo el análisis de los estudios de RNA-seq
- Debe ser lo tercero en ser ejecutado


## 3_Genera_todos_fit

- En esta carpeta se almacenan los datos de los fits de cada estudio
- Es necesario elegir si se quiere generar este objeto para todos los estudios, para los de sangre o para los de cerebro
- Debe ser lo cuarto en ser ejecutado


## 4_Metaanalisis

- En esta carpeta se encuentra el código encesario para realizar el metaanálisis
- También se encuentra el código necesario para obtener las listas de genes que se cargan en STRING
- Debe ser lo quinto en ser ejecutado


## 5_Analisis_funcional

- En esta carpeta se encuentra el código encesario para realizar el análisis funcional
- Debe ser lo sexto en ser ejecutado

## 6_Material_suplementario

### 6.1_Resultados_Revision_Sistematica

Se trata de un Excel que contiene la información de todo el proceso realizado en la Revisión Sistemática. Se pueden encontrar las siguientes pestañas:

- **Revisión Microarrays:** esta pestaña contiene la revisión sistemática de los estudios de microarrays
- **Revisión RNA-seq:** esta pestaña contiene la revisión sistemática de los estudios de RNA-seq
- **Selección final:** esta pestaña muestra la información de los estudios seleccioandos que se utilizan en este trabajo
- **Muestras estudios:** esta pestaña muestra el número de muestras que tiene cada grupo experimental de cada estudio
- **Resultado imputación:** esta pestaña muestra el porcentaje de acierto de la imputación de las muestras de los estudios de los microarrays
- **Muestras eliminadas:** esta pestaña muestra el número de muestras eliminadas y especifica las que se han eliminado para cada estudio
- **Expresión diferencial:** esta pestaña muestra el número de genes diferencialmente expresados para cada estudio (individual)
- **Metaanálisis:** esta pestaña muestra el número de genes diferencialmente expresado para cada metaanálisis y contraste.
- **OpenT_MA_todos:** esta pestaña muestra la informacióne estadistica de los genes diferencialmente expresados en el metaanálisis y que están presentes en la base de datos de Open Targets
  
### 6.2_Resultados_MA

- Se muestran las carpetas de los metaanálisis realizados
- Dentro de cada carpeta hay subcarpetas de los contrastes que tienen resultados significativos
- Dentro de cada subcarpeta se encuentra el informe html que permite ver los forest y funnel plots, así como el análisis de influencia para cada gen diferencialmente expresado en el metaanálisis.
