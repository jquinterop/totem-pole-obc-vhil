# Cargador Bidireccional Totem-Pole de Una Etapa — Modelos y Código de Soporte

Modelos de MATLAB/Simulink y scripts de soporte desarrollados para el trabajo de grado *"Diseño y Validación en Tiempo Real del Controlador de un Cargador Bidireccional de Vehículo Eléctrico de Una Etapa"* (Juan Pablo Quintero Pachón, Universidad de los Andes, Departamento de Ingeniería Eléctrica y Electrónica, 2026).

El proyecto adapta la topología totem-pole interleaved con aislamiento de alta frecuencia propuesta por Belkamel, Kim y Choi (2021) a las especificaciones de 120 V<sub>rms</sub>/60 Hz y 1 kW de la competencia IFEC 2026, y valida el sistema de control resultante mediante Virtual Hardware-in-the-Loop sobre un simulador OPAL-RT OP4500.

## Contenido del repositorio

| Archivo | Descripción |
|---|---|
| `OBC_Belkamel2021.slx` | Modelo de réplica del artículo original (220 V<sub>rms</sub>/3,3 kW), empleado para la validación presentada en la Sección 7.1 del documento. |
| `OBC_Tesis_IELE3002.slx` | Modelo offline de la adaptación (120 V<sub>rms</sub>/1 kW), con la planta de potencia en Simscape Electrical y el controlador discreto, empleado en las etapas de diseño de los Capítulos 4 y 5. |
| `OBC_Tesis_IELE3002_2023.slx` | Modelo de la adaptación preparado para ejecución en tiempo real sobre la plataforma OPAL-RT (Capítulo 6), empleado en la validación de los modos G2V y V2G presentada en las Secciones 7.2 y 7.3. |
| `script.m` | Script de MATLAB que calcula y exporta la estructura de parámetros correspondiente al caso de réplica (220 V/3,3 kW), consumida por `OBC_Belkamel2021.slx`. |
| `Script_Tesis.m` | Script de MATLAB que calcula y exporta la estructura de parámetros `belkamel_params` para el caso adaptado (120 V/1 kW): dimensionamiento de los componentes pasivos (Sección 4.7), sintonía de los tres lazos de control (Capítulo 5) y cálculo del punto de operación nominal, consumida por ambos modelos `OBC_Tesis_IELE3002*.slx`. |
| `LICENSE` | Licencia del repositorio. |
| `.gitignore` | Exclusión de archivos temporales y generados de MATLAB/Simulink (autosaves, `slprj/`, `codegen/`, etc.). |

## Requisitos de software

- **MATLAB/Simulink R2023b**, con las librerías **Simscape**, **Simscape Electrical** y **Simscape Battery**, para abrir y ejecutar `OBC_Belkamel2021.slx` y `OBC_Tesis_IELE3002.slx` en modo offline.
- **RT-LAB** (licencia 2024) junto con un simulador **OPAL-RT OP4500** con FPGA Kintex-7, para desplegar `OBC_Tesis_IELE3002_2023.slx` en tiempo real. Este modelo se guardó específicamente en formato R2023b por requisito de compatibilidad de la licencia de RT-LAB disponible durante el desarrollo del proyecto; **abrirlo con una versión de MATLAB distinta a R2023b puede impedir su despliegue correcto sobre RT-LAB**, aunque no debería afectar su apertura o edición en modo offline.

## Cómo reproducir los resultados

1. **Réplica (220 V/3,3 kW, Sección 7.1):** ejecutar `script.m` en MATLAB para generar la estructura de parámetros en el workspace, y luego abrir y simular `OBC_Belkamel2021.slx`.
2. **Adaptación, modo offline (Capítulos 4–5):** ejecutar `Script_Tesis.m`, luego abrir y simular `OBC_Tesis_IELE3002.slx`.
3. **Adaptación, modo tiempo real (Secciones 7.2–7.3):** ejecutar `Script_Tesis.m` en MATLAB R2023b, abrir `OBC_Tesis_IELE3002_2023.slx`, y seguir el procedimiento estándar de RT-LAB (*Build → Load → Execute*) para desplegarlo sobre el OPAL-RT OP4500, ajustando la consigna de tensión de batería desde la consola de monitoreo (`SC_eHS`) según el modo G2V o V2G deseado.

## Versión citada en el documento de tesis

La versión correspondiente exactamente a los resultados reportados en el Capítulo 7 del documento está fijada mediante la etiqueta (*tag*) `v1.0-tesis` de este repositorio.
