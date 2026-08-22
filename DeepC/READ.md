Un ejemplo simple e incompleto de DeepC.
- El modelo funciona bien, aunque requiere una etapa inicial de entrenamiento
- Problema:  un simple control proporcional funciona mejor
- El modelo de optimización se resuelve con JuMP. Se puede encontrar fácilmente una forma cerrada ya que es un problema cuadrático con restricciones de igualdad
- ¿Hay realmente alguna ventaja en usar este método?, Finalmente hay que almacenar muchos datos.  Nótese que con el modelo, solo debemos almacenar A,B,C.

