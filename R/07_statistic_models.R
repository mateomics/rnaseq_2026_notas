## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

colnames(mat)

## Hacemos la regresión lineal y obseramos valores ---> Y ~ X1 + X2
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))
