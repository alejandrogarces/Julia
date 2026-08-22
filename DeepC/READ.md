Un ejemplo simple e incompleto de DeepC.
- El modelo funciona bien, aunque requiere una etapa inicial de entrenamiento muy larga
- No resulta mejor a un simple control proporcional
- El modelo de optimización se resuelve con JuMP pero se puede encontrar fácilmente una forma analítica, ya que es un problema cuadrático con restricciones de igualdad
- ¿Hay alguna ventaja en usar este método, cuando al final hay que almacenar mas datos?  Nótese que con el modelo, solo debemos almacenar A,B,C.

