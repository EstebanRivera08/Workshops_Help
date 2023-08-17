import numpy as np
import numpy.linalg as la

A = np.array([[-2,3,9],[11,1,2],[4,6,-1]])


b = np.array([[45,4,10]])
b = b.T

# Gauss-Seidel

print("GAUSS-SEIDEL")

error = 10
it = 0
x = np.ones((3,1))

while (error > 1e-5) or (it < 10):
    x0 = x

    x1 = (45 - 3*x[1] - 9*x[2] )/-2
    x2 = (4 - 11*x1 - 2*x[2] )/1
    x3 = (10 - 6*x2 - 4*x1 )/-1

    x = np.array([x1,x2,x3])

    it = it + 1

    error = la.norm(x-x0, 2)


print("El número de iteraciones es: ",it)
print("La solución con Gauss-Seidel es:")
print(x)

