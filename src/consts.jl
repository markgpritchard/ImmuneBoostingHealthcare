
const SEIIRRRSTRANSITIONMATRIX = [
    # S   E   I1  I2  R1  R2  R3  cumulative cases 
     -1   1   0   0   0   0   0   1  # infection 
      0  -1   1   0   0   0   0   0  # transition E -> I1
      0   0  -1   1   0   0   0   0  # transition I1 -> I2
      0   0   0  -1   1   0   0   0  # transition I2 -> R1
      0   0   0   0  -1   1   0   0  # transition R1 -> R2
      0   0   0   0   0  -1   1   0  # transition R2 -> R3
      1   0   0   0   0   0  -1   0  # transition R3 -> S
      0   0   0   0   1  -1   0   0  # boosting from R2 
      0   0   0   0   1   0  -1   0  # boosting from R3 
     -1   0   0   0   1   0   0   0  # vaccination from S
      0   0   0   0   1  -1   0   0  # vaccination from R2
      0   0   0   0   1   0  -1   0  # vaccination from R3
]

WXYYZSEIIRRRSTRANSITIONMATRIX = [
    # W   X   Y1  Y2  Zr  Za  S   E   I1  I2  R1  R2  R3  ccp cch
     -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from W and admission to X 
     -1   0   1   0   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from W and admission to Y1 
     -1   0   0   1   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from W and admission to Y2 
     -1   0   0   0   0   1   0   0   0   0   0   0   0   0   0  # next discharge is from W and admission to Za 
      1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from X and admission to W 
      0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from X and admission to Y1 
      0  -1   0   1   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from X and admission to Y2 
      0  -1   0   0   0   1   0   0   0   0   0   0   0   0   0  # next discharge is from X and admission to Za 
      1   0  -1   0   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from Y1 and admission to W 
      0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from Y1 and admission to X 
      0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from Y1 and admission to Y2 
      0   0  -1   0   0   1   0   0   0   0   0   0   0   0   0  # next discharge is from Y1 and admission to Za 
      1   0   0  -1   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from Y2 and admission to W 
      0   1   0  -1   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from Y2 and admission to X 
      0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0  # next discharge is from Y2 and admission to Y1 
      0   0   0  -1   0   1   0   0   0   0   0   0   0   0   0  # next discharge is from Y2 and admission to Za 
      1   0   0   0  -1   0   0   0   0   0   0   0   0   0   0  # next discharge is from Zr and admission to W 
      0   1   0   0  -1   0   0   0   0   0   0   0   0   0   0  # next discharge is from Zr and admission to X 
      0   0   1   0  -1   0   0   0   0   0   0   0   0   0   0  # next discharge is from Zr and admission to Y1 
      0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0  # next discharge is from Zr and admission to Y2 
      0   0   0   0  -1   1   0   0   0   0   0   0   0   0   0  # next discharge is from Zr and admission to Za 
      1   0   0   0   0  -1   0   0   0   0   0   0   0   0   0  # next discharge is from Za and admission to W 
      0   1   0   0   0  -1   0   0   0   0   0   0   0   0   0  # next discharge is from Za and admission to X 
      0   0   1   0   0  -1   0   0   0   0   0   0   0   0   0  # next discharge is from Za and admission to Y1 
      0   0   0   1   0  -1   0   0   0   0   0   0   0   0   0  # next discharge is from Za and admission to Y2 
     -1   1   0   0   0   0   0   0   0   0   0   0   0   1   0  # patient infection 
      0  -1   1   0   0   0   0   0   0   0   0   0   0   0   0  # transition X -> Y1
      0   0  -1   1   0   0   0   0   0   0   0   0   0   0   0  # transition Y1 -> Y2
      0   0   0  -1   1   0   0   0   0   0   0   0   0   0   0  # transition Y2 -> Zr
      0   0   0   0   0   0  -1   1   0   0   0   0   0   0   1  # healthcare worker infection 
      0   0   0   0   0   0   0  -1   1   0   0   0   0   0   0  # transition E -> I1
      0   0   0   0   0   0   0   0  -1   1   0   0   0   0   0  # transition I1 -> I2
      0   0   0   0   0   0   0   0   0  -1   1   0   0   0   0  # transition I2 -> R1
      0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0  # transition R1 -> R2
      0   0   0   0   0   0   0   0   0   0   0  -1   1   0   0  # transition R2 -> R3
      0   0   0   0   0   0   1   0   0   0   0   0  -1   0   0  # transition R3 -> S
      0   0   0   0   0   0   0   0   0   0   1  -1   0   0   0  # boosting from R2 
      0   0   0   0   0   0   0   0   0   0   1   0  -1   0   0  # boosting from R3 
      0   0   0   0   0   0  -1   0   0   0   1   0   0   0   0  # vaccination from S
      0   0   0   0   0   0   0   0   0   0   1  -1   0   0   0  # vaccination from R2
      0   0   0   0   0   0   0   0   0   0   1   0  -1   0   0  # vaccination from R3
]
