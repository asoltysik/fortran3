from mylib import mylib
import numpy as np

print("Matrix multiplication sequential:")

a = np.matrix('1 2; 3 4')
b = np.matrix('5 6; 7 8')

c = mylib.matmul(a, b)
print(c)

print("\n\nGauss elimination:")

a = np.array([[1.0, -2.753, 1.0], [2.0, -0.3, 3.0], [-1.5, 1.3, -4.2]], order='F')
x = np.array(np.random.rand(3), order='F')

print(f"a = {a}")
print(f"x = {x}")

ret = np.linalg.solve(a, x)
print("\nsolved by numpy: ")
print(ret)

mylib.gauss(a, x, 3)

print("\nSolved my my fortran function (not even close): ")
print(a)
print(x)
