
const SEIIRRRSTRANSITIONMATRIX = [
    # S   E   I1  I2  R1  R2  R3  cumulative cases 
     -1   1   0   0   0   0   0   1  # infection 
      0  -1   1   0   0   0   0   0  # transition E -> I1
      0   0  -1   1   0   0   0   0  # transition I1 -> I2
      0   0   0  -1   1   0   0   0  # transition I2 -> R1
      0   0   0   0  -1   1   0   0  # transition R1 -> R2
      0   0   0   0   0  -1   1   0  # transition R2 -> R3
      1   0   0   0   0   0  -1   0  # transition R3 -> S
      0   0   0   0   1  -1   0   0  # boosting plus vaccination from R2 
      0   0   0   0   1   0  -1   0  # boosting plus vaccination from R3 
     -1   0   0   0   1   0   0   0  # vaccination from S
]

WXYYZSEIIRRRSTRANSITIONMATRIX = [
    # W   X   Y1  Y2  Zr  Za  S   E   I1  I2  I′  R1  R2  R3  ccp cch dia
     -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from W and admission to X 
     -1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from W and admission to Y1 
     -1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from W and admission to Y2 
     -1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0  # discharge from W and admission to Za 
      1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from X and admission to W 
      0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from X and admission to Y1 plus transition X -> Y1
      0  -1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from X and admission to Y2 
      0  -1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0  # discharge from X and admission to Za 
      1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y1 and admission to W 
      0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y1 and admission to X 
      0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y1 and admission to Y2 plus transition Y1 -> Y2
      0   0  -1   0   0   1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y1 and admission to Za 
      1   0   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y2 and admission to W 
      0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y2 and admission to X 
      0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y2 and admission to Y1 
      0   0   0  -1   0   1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Y2 and admission to Za 
      1   0   0   0  -1   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Zr and admission to W 
      0   1   0   0  -1   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Zr and admission to X 
      0   0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Zr and admission to Y1 
      0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0  # discharge from Zr and admission to Y2 
      0   0   0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Zr and admission to Za 
      1   0   0   0   0  -1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Za and admission to W 
      0   1   0   0   0  -1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Za and admission to X 
      0   0   1   0   0  -1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Za and admission to Y1 
      0   0   0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0  # discharge from Za and admission to Y2 
     -1   1   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0  # patient infection 
      0   0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0  # transition Y2 -> Zr
      0   0   0   0   0   0  -1   1   0   0   0   0   0   0   0   1   0  # healthcare worker infection 
      0   0   0   0   0   0   0  -1   1   0   0   0   0   0   0   0   0  # transition E -> I1
      0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0   0   0  # transition I1 -> I2
      0   0   0   0   0   0   0   0  -1   0   1   0   0   0   0   0   1  # diagnosis from I1
      0   0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0   1  # diagnosis from I2
      0   0   0   0   0   0   0   0   0  -1   0   1   0   0   0   0   0  # transition I2 -> R1
      0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0  # transition I′ -> R1
      0   0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0   0  # transition R1 -> R2
      0   0   0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0  # transition R2 -> R3
      0   0   0   0   0   0   1   0   0   0   0   0   0  -1   0   0   0  # transition R3 -> S
      0   0   0   0   0   0   0   0   0   0   0   1  -1   0   0   0   0  # boosting from R2 plus vaccination from R2
      0   0   0   0   0   0   0   0   0   0   0   1   0  -1   0   0   0  # boosting from R3 plus vaccination from R3
      0   0   0   0   0   0  -1   0   0   0   0   1   0   0   0   0   0  # vaccination from S
]
